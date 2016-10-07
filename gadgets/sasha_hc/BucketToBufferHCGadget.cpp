#include "BucketToBufferHCGadget.h"
namespace Gadgetron {

    BucketToBufferHCGadget::BucketToBufferHCGadget()
    {
    }

    BucketToBufferHCGadget::~BucketToBufferHCGadget()
    {
        //The buckets array should be empty but just in case, let's make sure all the stuff is released.
    }

    int BucketToBufferHCGadget
        ::process( GadgetContainerMessage<IsmrmrdAcquisitionBucket>* m1 )
    {

        size_t key;
        std::map<size_t, GadgetContainerMessage<IsmrmrdReconData>* > recon_data_buffers;

        //GDEBUG("BucketToBufferGadget::process\n");

        //Some information about the bucket
        //GDEBUG_STREAM("The Reference part: " << m1->getObjectPtr()->refstats_.size() << std::endl);
        //GDEBUG_STREAM("   nslices: " << m1->getObjectPtr()->refstats_[0].slice.size() << std::endl);
        //for (int e=0; e<m1->getObjectPtr()->refstats_.size() ; e++) {
        //    for (std::set<uint16_t>::iterator it = m1->getObjectPtr()->refstats_[e].kspace_encode_step_1.begin();
        //         it != m1->getObjectPtr()->refstats_[e].kspace_encode_step_1.end(); ++it) {
        //        GDEBUG_STREAM("   K1: " <<  *it << std::endl);
        //    }
        //}
        //GDEBUG_STREAM("The data part: " << m1->getObjectPtr()->datastats_.size() << std::endl);
        //GDEBUG_STREAM("   nslices: " << m1->getObjectPtr()->datastats_[0].slice.size() << std::endl);
        //for (int e=0; e<m1->getObjectPtr()->datastats_.size() ; e++) {
        //    for (std::set<uint16_t>::iterator it = m1->getObjectPtr()->datastats_[e].kspace_encode_step_1.begin();
        //         it != m1->getObjectPtr()->datastats_[e].kspace_encode_step_1.end(); ++it) {
        //        GDEBUG_STREAM("   K1: " <<  *it << std::endl);
        //    }
        //}

        //Iterate over the reference data of the bucket
        IsmrmrdDataBuffered* pCurrDataBuffer = NULL;
        for (std::vector<IsmrmrdAcquisitionData>::iterator it = m1->getObjectPtr()->ref_.begin();
        it != m1->getObjectPtr()->ref_.end(); ++it)
        {
            //Get a reference to the header for this acquisition
            ISMRMRD::AcquisitionHeader & acqhdr = *it->head_->getObjectPtr();

            //Generate the key to the corresponding ReconData buffer
            key = getKey( acqhdr.idx );

            //The storage is based on the encoding space
            uint16_t espace = acqhdr.encoding_space_ref;

            //GDEBUG_STREAM("espace: " << acqhdr.encoding_space_ref << std::endl);
            //GDEBUG_STREAM("slice: " << acqhdr.idx.slice << std::endl);
            //GDEBUG_STREAM("rep: " << acqhdr.idx.repetition << std::endl);
            //GDEBUG_STREAM("k1: " << acqhdr.idx.kspace_encode_step_1 << std::endl);
            //GDEBUG_STREAM("k2: " << acqhdr.idx.kspace_encode_step_2 << std::endl);
            //GDEBUG_STREAM("seg: " << acqhdr.idx.segment << std::endl);
            //GDEBUG_STREAM("key: " << key << std::endl);

            //Get some references to simplify the notation
            //the reconstruction bit corresponding to this ReconDataBuffer and encoding space
            IsmrmrdReconBit & rbit = getRBit( recon_data_buffers, key, espace );
            //and the corresponding data buffer for the reference data
            if (!rbit.ref_)
                rbit.ref_ = IsmrmrdDataBuffered();
            IsmrmrdDataBuffered & dataBuffer = *rbit.ref_;
            //this encoding space's xml header info
            ISMRMRD::Encoding & encoding = hdr_.encoding[espace];
            //this bucket's reference stats
            IsmrmrdAcquisitionBucketStats & stats = m1->getObjectPtr()->refstats_[espace];

            //Fill the sampling description for this data buffer, only need to fill the sampling_ once per recon bit
            if (&dataBuffer != pCurrDataBuffer)
            {
                fillSamplingDescription( dataBuffer.sampling_, encoding, stats, acqhdr, true );
                pCurrDataBuffer = &dataBuffer;
            }

            //Make sure that the data storage for this data buffer has been allocated
            //TODO should this check the limits, or should that be done in the stuff function?
            allocateDataArrays( dataBuffer, acqhdr, encoding, stats, true );

            // Stuff the data, header and trajectory into this data buffer
            stuff( it, dataBuffer, encoding, true );
        }


        // 2016-10-03 Kelvin: Find the limits of the high-contrast lines
        int16_t iHCLineMin = (int16_t)6536;
        int16_t iHCLineMax = 0;

        //Iterate over the imaging data of the bucket
        // this is exactly the same code as for the reference data except for
        // the chunk of the data buffer.
        pCurrDataBuffer = NULL;
        for (std::vector<IsmrmrdAcquisitionData>::iterator it = m1->getObjectPtr()->data_.begin();
        it != m1->getObjectPtr()->data_.end(); ++it)
        {
            //Get a reference to the header for this acquisition
            ISMRMRD::AcquisitionHeader & acqhdr = *it->head_->getObjectPtr();

            //Generate the key to the corresponding ReconData buffer
            key = getKey( acqhdr.idx );

            //The storage is based on the encoding space
            uint16_t espace = acqhdr.encoding_space_ref;

            //GDEBUG_STREAM("espace: " << acqhdr.encoding_space_ref << std::endl);
            //GDEBUG_STREAM("slice: " << acqhdr.idx.slice << std::endl);
            //GDEBUG_STREAM("rep: " << acqhdr.idx.repetition << std::endl);
            //GDEBUG_STREAM("k1: " << acqhdr.idx.kspace_encode_step_1 << std::endl);
            //GDEBUG_STREAM("k2: " << acqhdr.idx.kspace_encode_step_2 << std::endl);
            //GDEBUG_STREAM("seg: " << acqhdr.idx.segment << std::endl);
            //GDEBUG_STREAM("key: " << key << std::endl);

            //Get some references to simplify the notation
            //the reconstruction bit corresponding to this ReconDataBuffer and encoding space
            IsmrmrdReconBit & rbit = getRBit( recon_data_buffers, key, espace );
            //and the corresponding data buffer for the imaging data
            IsmrmrdDataBuffered & dataBuffer = rbit.data_;
            //this encoding space's xml header info
            ISMRMRD::Encoding & encoding = hdr_.encoding[espace];
            //this bucket's imaging data stats
            IsmrmrdAcquisitionBucketStats & stats = m1->getObjectPtr()->datastats_[espace];

            //Fill the sampling description for this data buffer, only need to fill sampling_ once per recon bit
            if (&dataBuffer != pCurrDataBuffer)
            {
                fillSamplingDescription( dataBuffer.sampling_, encoding, stats, acqhdr, false );
                pCurrDataBuffer = &dataBuffer;
            }

            //Make sure that the data storage for this data buffer has been allocated
            //TODO should this check the limits, or should that be done in the stuff function?
            allocateDataArrays( dataBuffer, acqhdr, encoding, stats, false );

            // 2016-10-03 Kelvin: Update limits for HC images
            if (acqhdr.idx.contrast)
            {
                if ((int16_t)acqhdr.idx.kspace_encode_step_1 < iHCLineMin)
                {
                    iHCLineMin = (int16_t)acqhdr.idx.kspace_encode_step_1;
                    //				GDEBUG_STREAM("Changed iHCLineMin: set: " << (int16_t)acqhdr.idx.set << " lin: " << (int16_t)acqhdr.idx.kspace_encode_step_1 << std::endl);
                }

                if ((int16_t)acqhdr.idx.kspace_encode_step_1 > iHCLineMax)
                {
                    iHCLineMax = (int16_t)acqhdr.idx.kspace_encode_step_1;
                }
            }

            // Stuff the data, header and trajectory into this data buffer
            stuff( it, dataBuffer, encoding, false );

            // 2016-09-28 Kelvin: Duplicate regular images for HC
        /*
            if (((int16_t)acqhdr.idx.user[1] == 0)
              && (((int16_t)acqhdr.idx.kspace_encode_step_1 > 65)))
        //	  && (((int16_t)acqhdr.idx.kspace_encode_step_1 < 23) || ((int16_t)acqhdr.idx.kspace_encode_step_1 > 65)))
            {
                acqhdr.idx.contrast = 1;
                stuff(it, dataBuffer, encoding, false);
                acqhdr.idx.contrast = 0;
            }
        */
        }

        // 2016-10-03 Kelvin: Re-stuff shared lines
        GDEBUG_STREAM( "iHCLineMin: " << iHCLineMin << " iHCLineMax: " << iHCLineMax );
        pCurrDataBuffer = NULL;
        for (std::vector<IsmrmrdAcquisitionData>::iterator it = m1->getObjectPtr()->data_.begin();
        it != m1->getObjectPtr()->data_.end(); ++it)
        {
            //Get a reference to the header for this acquisition
            ISMRMRD::AcquisitionHeader & acqhdr = *it->head_->getObjectPtr();

            //Generate the key to the corresponding ReconData buffer
            key = getKey( acqhdr.idx );

            //The storage is based on the encoding space
            uint16_t espace = acqhdr.encoding_space_ref;

            //GDEBUG_STREAM("espace: " << acqhdr.encoding_space_ref << std::endl);
            //GDEBUG_STREAM("slice: " << acqhdr.idx.slice << std::endl);
            //GDEBUG_STREAM("rep: " << acqhdr.idx.repetition << std::endl);
            //GDEBUG_STREAM("k1: " << acqhdr.idx.kspace_encode_step_1 << std::endl);
            //GDEBUG_STREAM("k2: " << acqhdr.idx.kspace_encode_step_2 << std::endl);
            //GDEBUG_STREAM("seg: " << acqhdr.idx.segment << std::endl);
            //GDEBUG_STREAM("key: " << key << std::endl);

            //Get some references to simplify the notation
            //the reconstruction bit corresponding to this ReconDataBuffer and encoding space
            IsmrmrdReconBit & rbit = getRBit( recon_data_buffers, key, espace );
            //and the corresponding data buffer for the imaging data
            IsmrmrdDataBuffered & dataBuffer = rbit.data_;
            //this encoding space's xml header info
            ISMRMRD::Encoding & encoding = hdr_.encoding[espace];
            //this bucket's imaging data stats
            IsmrmrdAcquisitionBucketStats & stats = m1->getObjectPtr()->datastats_[espace];

            //Fill the sampling description for this data buffer, only need to fill sampling_ once per recon bit
            if (&dataBuffer != pCurrDataBuffer)
            {
                fillSamplingDescription( dataBuffer.sampling_, encoding, stats, acqhdr, false );
                pCurrDataBuffer = &dataBuffer;
            }

            //Make sure that the data storage for this data buffer has been allocated
            //TODO should this check the limits, or should that be done in the stuff function?
            allocateDataArrays( dataBuffer, acqhdr, encoding, stats, false );

            // 2016-10-03 Kelvin: Share lines
            /*
            if (acqhdr.idx.contrast == 0)
            {
                if (((int16_t)acqhdr.idx.kspace_encode_step_1 < iHCLineMin) || ((int16_t)acqhdr.idx.kspace_encode_step_1 > iHCLineMax))
                {
                    acqhdr.idx.contrast = 1;
                    stuff(it, dataBuffer, encoding, false);
                    acqhdr.idx.contrast = 0;

                    GDEBUG_STREAM("     sharing set: " << (int16_t)acqhdr.idx.set << " lin: " << (int16_t)acqhdr.idx.kspace_encode_step_1);
                } else {
                    GDEBUG_STREAM(" Not sharing set: " << (int16_t)acqhdr.idx.set << " lin: " << (int16_t)acqhdr.idx.kspace_encode_step_1);
                }
            } else {
                GDEBUG_STREAM("*Not sharing set: " << (int16_t)acqhdr.idx.set << " lin: " << (int16_t)acqhdr.idx.kspace_encode_step_1);
            }
            */
        }

        //Send all the ReconData messages
        GDEBUG( "End of bucket reached, sending out %d ReconData buffers\n", recon_data_buffers.size() );
        for (std::map<size_t, GadgetContainerMessage<IsmrmrdReconData>* >::iterator it = recon_data_buffers.begin(); it != recon_data_buffers.end(); it++)
        {
            //GDEBUG_STREAM("Sending: " << it->first << std::endl);
            if (it->second) {
                if (this->next()->putq( it->second ) == -1) {
                    it->second->release();
                    throw std::runtime_error( "Failed to pass bucket down the chain\n" );
                }
            }
        }

        //Clear the recondata buffer map
        recon_data_buffers.clear();  // is this necessary?

        //We can release the incoming bucket now. This will release all of the data it contains.
        m1->release();

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE( BucketToBufferHCGadget )

}
