#include "DistributeGadget.h"
#include "GadgetStreamInterface.h"
#include "gadgetron_xml.h"
#include "CloudBus.h"
#include <stdint.h>

namespace Gadgetron{

  DistributionConnector::DistributionConnector(DistributeGadget* g)
  : distribute_gadget_(g)
  {

  }

  int DistributionConnector::process(size_t messageid, ACE_Message_Block* mb) {
    return distribute_gadget_->collector_putq(mb);
  }


  DistributeGadget::DistributeGadget()
  : BasicPropertyGadget()
  , mtx_("distribution_mtx")
  , prev_connector_(0)
  {
  }


  const char* DistributeGadget::get_node_xml_config()
  {
    return node_xml_config_.c_str();
  }

  int DistributeGadget::collector_putq(ACE_Message_Block* m)
  {
    if (!collect_gadget_) {
      GERROR("Collector gadget not set\n");
      m->release();
      return GADGET_FAIL;
    }

    if (collect_gadget_->putq(m) == -1) {
      m->release();
      GERROR("DistributeGadget::collector_putq, passing data on to next gadget\n");
      return GADGET_FAIL;
    }
    return GADGET_OK;
  }

  int DistributeGadget::process(ACE_Message_Block* m)
  {
    int node_index = this->node_index(m);

    if (single_package_mode.value()) {
      node_index = ++started_nodes_;
    }

    if (node_index < 0) {
      GERROR("Negative node index received");
      return GADGET_FAIL;
    }

    //If we are not supposed to use this node for compute, add one to make sure we are not on node 0
    //if (!use_this_node_for_compute.value()) {
    //  node_index = node_index+1;
    //}

    // instead of sending down the stream, processing is done by making connections
    //if (node_index == 0) { //process locally
    //  if (this->next()->putq(m) == -1) {
    //    m->release();
    //    GERROR("DistributeGadget::process, passing data on to next gadget\n");
    //    return GADGET_FAIL;
    //  }
    //  return GADGET_OK;
    //}

    //At this point, the node index is positive, so we need to find a suitable connector.
    mtx_.acquire();
    auto n = node_map_.find(node_index);
    mtx_.release();
    GadgetronConnector* con = 0;
    if (n != node_map_.end()) { //We have a suitable connection already.
      con = n->second;
    } else {
      std::vector<GadgetronNodeInfo> nl;
      CloudBus::instance()->get_node_info(nl);

      GDEBUG("Number of network nodes found: %d\n", nl.size());

      GadgetronNodeInfo me;
      me.address = "127.0.0.1";//We may have to update this
      me.port = CloudBus::instance()->port();
      me.uuid = CloudBus::instance()->uuid();
      me.active_reconstructions = CloudBus::instance()->active_reconstructions();

      //This would give the current node the lowest possible priority
      if (!use_this_node_for_compute.value()) {
        me.active_reconstructions = UINT32_MAX;
      }

      for (auto it = nl.begin(); it != nl.end(); it++) {
        if (it->active_reconstructions < me.active_reconstructions) {
          me = *it;
        }

        //Is this a free node
        if (me.active_reconstructions == 0) break;
      }

      // first job, send to current node if required
      if (use_this_node_for_compute.value() && node_index==0)
      {
        size_t num_of_ip = local_address_.size();

          for (auto it = nl.begin(); it != nl.end(); it++)
          {
              for (size_t ii=0; ii<num_of_ip; ii++)
              {
                  if (it->address == local_address_[ii])
                  {
                      me = *it;
                  }
              }
          }

          GDEBUG_STREAM("Send first job to current node : " << me.address);
      }

      con = new DistributionConnector(this);

      GadgetronXML::GadgetStreamConfiguration cfg;
      try {
        deserialize(node_xml_config_.c_str(), cfg);
      }  catch (const std::runtime_error& e) {
        GERROR("Failed to parse Node Gadget Stream Configuration: %s\n", e.what());
        return GADGET_FAIL;
      }

      //Configuration of readers
      for (auto i = cfg.reader.begin(); i != cfg.reader.end(); ++i) {
        GadgetMessageReader* r =
        controller_->load_dll_component<GadgetMessageReader>(i->dll.c_str(),
        i->classname.c_str());
        if (!r) {
          GERROR("Failed to load GadgetMessageReader from DLL\n");
          return GADGET_FAIL;
        }
        con->register_reader(i->slot, r);
      }

      for (auto i = cfg.writer.begin(); i != cfg.writer.end(); ++i) {
        GadgetMessageWriter* w =
        controller_->load_dll_component<GadgetMessageWriter>(i->dll.c_str(),
        i->classname.c_str());
        if (!w) {
          GERROR("Failed to load GadgetMessageWriter from DLL\n");
          return GADGET_FAIL;
        }
        con->register_writer(i->slot, w);
      }


      char buffer[10];
      sprintf(buffer,"%d",me.port);
      if (con->open(me.address,std::string(buffer)) != 0) {
        GERROR("Failed to open connection to node %s : %d\n", me.address.c_str(), me.port);
        return GADGET_FAIL;
      }

      if (con->send_gadgetron_configuration_script(node_xml_config_) != 0) {
        GERROR("Failed to send XML configuration to compute node\n");
        return GADGET_FAIL;
      }

      if (con->send_gadgetron_parameters(node_parameters_) != 0) {
        GERROR("Failed to send XML parameters to compute node\n");
        return GADGET_FAIL;
      }
      
      mtx_.acquire();
      node_map_[node_index] = con;
      mtx_.release();
    }


    if (!con) {
      //Zero pointer for the connection means that either a) connection creation failed or b) using local chain.
      //Either way, we will send it down the chain.

      if (!use_this_node_for_compute.value()) {
        GERROR("This node cannot be used for computing and no other node is available\n");
        m->release();
        return GADGET_FAIL;
      }

      if (this->next()->putq(m) == -1) {
        m->release();
        GERROR("DistributeGadget::process, passing data on to next gadget\n");
        return GADGET_FAIL;
      } else {
        return GADGET_OK;
      }

    } else {
     
      //Let's make sure that we did not send a close message to this connector already
      auto c = std::find(closed_connectors_.begin(),closed_connectors_.end(),con);
      if (c != closed_connectors_.end()) {
	//This is a bad situation, we need to bail out. 
	m->release();
	GERROR("The valid connection for incoming data has already been closed. Distribute Gadget is not configured properly for this type of data\n");
	return GADGET_FAIL;
      }

      //If nodes receive their data sequentially (default), we should see if we should be closing the previos connection
      if (nodes_used_sequentially.value() && !single_package_mode.value()) {
	//Is this a new connection, if so, send previous one a close
	if (prev_connector_ && prev_connector_ != con) {
	  GDEBUG("Sending close to previous connector, not expecting any more data for this one\n");
	  auto mc = new GadgetContainerMessage<GadgetMessageIdentifier>();
	  mc->getObjectPtr()->id = GADGET_MESSAGE_CLOSE;
	  
	  if (prev_connector_->putq(mc) == -1) {
	    GERROR("Unable to put CLOSE package on queue of previous connection\n");
	    return -1;
	  }
	  closed_connectors_.push_back(prev_connector_);
	}
      }
     
      //Update previous connection
      prev_connector_ = con;

      //We have a valid connector
      auto m1 = new GadgetContainerMessage<GadgetMessageIdentifier>();
      m1->getObjectPtr()->id = message_id(m);

      m1->cont(m);

      if (con->putq(m1) == -1) {
        GERROR("Unable to put package on connector queue\n");
        m1->release();
        return GADGET_FAIL;
      }

      if (single_package_mode.value()) {
        auto m2 = new GadgetContainerMessage<GadgetMessageIdentifier>();
        m2->getObjectPtr()->id = GADGET_MESSAGE_CLOSE;

        if (con->putq(m2) == -1) {
          GERROR("Unable to put CLOSE package on queue\n");
          return -1;
        }
	closed_connectors_.push_back(con);
      }
    }


    return 0;
  }

  int DistributeGadget::process_config(ACE_Message_Block* m)
  {

    started_nodes_ = 0;
    node_parameters_ = std::string(m->rd_ptr());

    //Grab the original XML conifguration
    std::string xml = controller_->get_xml_configuration();

    GadgetronXML::GadgetStreamConfiguration cfg;
    GadgetronXML::deserialize(xml.c_str(),cfg);

    //Delete Gadgets up to this Gadget
    std::vector<GadgetronXML::Gadget>::iterator it = cfg.gadget.begin();
    while ((it->name != std::string(this->module()->name())) && (it != cfg.gadget.end())) it++; it++;
    cfg.gadget.erase(cfg.gadget.begin(),it);

    //Delete Gadgets after collector
    it = cfg.gadget.begin();
    while ((it->name != collector.value()) && (it != cfg.gadget.end())) it++; it++;
    cfg.gadget.erase(it,cfg.gadget.end());

    std::stringstream o;
    GadgetronXML::serialize(cfg,o);

    node_xml_config_ = o.str();

    Gadget* tmp = this;
    while (tmp->next()) {
      if (std::string(tmp->module()->name()) == collector.value()) break;
      tmp = dynamic_cast<Gadget*>(tmp->next());
    }

    collect_gadget_ = tmp;

    if (!collect_gadget_) {
      GERROR("Failed to locate collector Gadget with name %s\n", collector.value().c_str());
      return GADGET_FAIL;
    } else {
      collect_gadget_->set_parameter("pass_through_mode","true");
    }

    // get current node ip addresses
    ACE_INET_Addr* the_addr_array = NULL;
    size_t num_of_ip = 0;

    int rc = ACE::get_ip_interfaces (num_of_ip, the_addr_array);
    if (rc != 0)
    {
        GERROR_STREAM("Retreive local ip addresses failed ... ");
        num_of_ip = 0;
    }

    if (the_addr_array!=NULL ) delete [] the_addr_array;

    for (size_t ii=0; ii<num_of_ip; ii++)
    {
        std::string ip = std::string(the_addr_array[ii].get_host_addr());
        local_address_.push_back(ip);
        GDEBUG_STREAM("--> Local address  : " << ip);
    }

    return GADGET_OK;
  }

  int DistributeGadget::close(unsigned long flags)
  {
    int ret = Gadget::close(flags);
    if (flags) {
      mtx_.acquire();

      for (auto n = node_map_.begin(); n != node_map_.end(); n++) {
        if (n->second) {
          auto m1 = new GadgetContainerMessage<GadgetMessageIdentifier>();
          m1->getObjectPtr()->id = GADGET_MESSAGE_CLOSE;

          if (n->second->putq(m1) == -1) {
            GERROR("Unable to put CLOSE package on queue\n");
            return -1;
          }
	}
      }

      auto it = node_map_.begin();
      while (it != node_map_.end()) {
        if (it->second) {
          it->second->wait();
          delete it->second;
        }
        node_map_.erase(it);
        it = node_map_.begin();
      }

      mtx_.release();
      GDEBUG("All connectors closed. Waiting for Gadget to close\n");
    }
    return ret;
  }

  int DistributeGadget::node_index(ACE_Message_Block* m)
  {
    return 0;
  }

  GADGET_FACTORY_DECLARE(DistributeGadget)
}
