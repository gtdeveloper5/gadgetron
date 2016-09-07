#include "fatwater.h"

#include "hoMatrix.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "hoArmadillo.h"

#include <math.h>
#include <iostream>
#include <boost/config.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/timer/timer.hpp>

// Curve fitting includes (from Hui's example)
//#include "hoNDHarrWavelet.h"
//#include "hoNDRedundantWavelet.h"
#include "hoNDArray_math.h"
#include "simplexLagariaSolver.h"
#include "twoParaExpDecayOperator.h"
#include "curveFittingCostFunction.h"
// #include <gtest/gtest.h>
#include <boost/random.hpp>

#define GAMMABAR 42.576 // MHz/T
#define PI 3.141592

using namespace boost;
 
typedef int EdgeWeightType;
 
typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
typedef adjacency_list < vecS, vecS, directedS,
			 property < vertex_name_t, std::string,
				    property < vertex_index_t, long,
			 property < vertex_color_t, boost::default_color_type,
			 property < vertex_distance_t, long,
			 property < vertex_predecessor_t, Traits::edge_descriptor > > > > >,
			 property < edge_capacity_t, EdgeWeightType,
			 property < edge_residual_capacity_t, EdgeWeightType,
			 property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;

/*
Traits::edge_descriptor AddEdge(Traits::vertex_descriptor &v1,
				Traits::vertex_descriptor &v2,
				property_map < Graph, edge_reverse_t >::type &rev,
				const double capacity,
				Graph &g);
*/

void AddEdge(Traits::vertex_descriptor &v1, Traits::vertex_descriptor &v2, property_map < Graph, edge_reverse_t >::type &rev, const double capacity, Graph &g)
{
  Traits::edge_descriptor e1 = add_edge(v1, v2, g).first;
  Traits::edge_descriptor e2 = add_edge(v2, v1, g).first;
  put(edge_capacity, g, e1, capacity);
  put(edge_capacity, g, e2, 0*capacity);
  
  rev[e1] = e2;
  rev[e2] = e1;
}


namespace Gadgetron {
  hoNDArray< std::complex<float> > fatwater_separation(hoNDArray< std::complex<float> >& data, FatWaterParameters p, FatWaterAlgorithm a)
  {
    


    // Start the timer
    boost::timer::cpu_timer mytimer;
    //    mytimer.start();


    //Get some data parameters
    //7D, fixed order [X, Y, Z, CHA, N, S, LOC]
    uint16_t X = data.get_size(0);
    uint16_t Y = data.get_size(1);
    uint16_t Z = data.get_size(2);
    uint16_t CHA = data.get_size(3);
    uint16_t N = data.get_size(4);
    uint16_t S = data.get_size(5);
    uint16_t LOC = data.get_size(6);
    
    hoNDArray< std::complex<float> > out(X,Y,Z,CHA,N,2,LOC); // S dimension gets replaced by water/fat stuff
    
    float fieldStrength = p.fieldStrengthT_;
    std::vector<float> echoTimes = p.echoTimes_;
    bool precessionIsClockwise = p.precessionIsClockwise_;
    for (auto& te: echoTimes) {
      te = te*0.001; // Echo times in seconds rather than milliseconds
    }
    
    //Get or set some algorithm parameters
    //Gadgetron::ChemicalSpecies w = a.species_[0];
    //Gadgetron::ChemicalSpecies f = a.species_[1];
    
    //	GDEBUG("In toolbox - Fat peaks: %f  \n", f.ampFreq_[0].first);
    //	GDEBUG("In toolbox - Fat peaks 2: %f  \n", f.ampFreq_[0].second);
    
    // Set some initial parameters so we can get going
    // These will have to be specified in the XML file eventually
    std::pair<float,float> range_r2star = std::make_pair(0.0,100.0);
    uint16_t num_r2star = 11;
    std::pair<float,float> range_fm = std::make_pair(-200.0,200.0);
    uint16_t num_fm = 201;
    uint16_t size_clique = 1;
    uint16_t num_iterations = 20;
    uint16_t subsample = 1;
    float lmap_power = 2.0;
    float lambda = 0.02;
    float lambda_extra = 0.01;
    
    //Check that we have reasonable data for fat-water separation
    
    
    //Calculate residual
    //
    float relAmp, freq_hz;
    uint16_t npeaks;
    uint16_t nspecies = a.species_.size();
    uint16_t nte = echoTimes.size();
    
    hoMatrix< std::complex<float> > phiMatrix(nte,nspecies);
    for( int k1=0;k1<nte;k1++) {
      for( int k2=0;k2<nspecies;k2++) {
	phiMatrix(k1,k2) = 0.0;
	npeaks = a.species_[k2].ampFreq_.size();
	for( int k3=0;k3<npeaks;k3++) {
	  relAmp = a.species_[k2].ampFreq_[k3].first;
	  freq_hz = fieldStrength*GAMMABAR*a.species_[k2].ampFreq_[k3].second;
	  phiMatrix(k1,k2) += relAmp*std::complex<float>(cos(2*PI*echoTimes[k1]*freq_hz),sin(2*PI*echoTimes[k1]*freq_hz));
	}
      }
    }
    //auto a_phiMatrix = as_arma_matrix(&phiMatrix);
    //auto mymat2 = mymat.t()*mymat;
    
    hoMatrix< std::complex<float> > IdentMat(nte,nte);
    for( int k1=0;k1<nte;k1++) {
      for( int k2=0;k2<nte;k2++) {
	if( k1==k2 ) {
	  IdentMat(k1,k2) = std::complex<float>(1.0,0.0);
	} else {
	  IdentMat(k1,k2) = std::complex<float>(0.0,0.0);
	}
      }
    }
    //	auto a_phiMatrix = as_arma_matrix(&IdentMat);
    
    float fm;
    std::vector<float> fms(num_fm);
    fms[0] = range_fm.first;
    for(int k1=1;k1<num_fm;k1++) {
      fms[k1] = range_fm.first + k1*(range_fm.second-range_fm.first)/(num_fm-1);
    }
    
    float r2star;
    std::vector<float> r2stars(num_r2star);
    r2stars[0] = range_r2star.first;
    for(int k2=1;k2<num_r2star;k2++) {
      r2stars[k2] = range_r2star.first + k2*(range_r2star.second-range_r2star.first)/(num_r2star-1);
    }
    
    
    std::complex<float> curModulation;
    hoMatrix< std::complex<float> > tempM1(nspecies,nspecies);
    hoMatrix< std::complex<float> > tempM2(nspecies,nte);
    hoMatrix< std::complex<float> > psiMatrix(nte,nspecies);
    hoNDArray< std::complex<float> > Ps(nte,nte,num_fm,num_r2star);
    hoMatrix< std::complex<float> > P1(nte,nte);
    hoMatrix< std::complex<float> > P(nte,nte);
    
    for(int k3=0;k3<num_fm;k3++) {
      fm = fms[k3];
      for(int k4=0;k4<num_r2star;k4++) {
	r2star = r2stars[k4];
	
	
	for( int k1=0;k1<nte;k1++) {
	  curModulation = exp(-r2star*echoTimes[k1])*std::complex<float>(cos(2*PI*echoTimes[k1]*fm),sin(2*PI*echoTimes[k1]*fm));
	  for( int k2=0;k2<nspecies;k2++) {
	    psiMatrix(k1,k2) = phiMatrix(k1,k2)*curModulation;
	  }
	}
	
	herk( tempM1, psiMatrix, 'L', true );
	//	    tempM1.copyLowerTriToUpper();
	for (int ka=0;ka<tempM1.get_size(0);ka++ ) {
	  for (int kb=ka+1;kb<tempM1.get_size(1);kb++ ) {
	    tempM1(ka,kb) = conj(tempM1(kb,ka));
	  }
	}
	
	potri(tempM1);
	for (int ka=0;ka<tempM1.get_size(0);ka++ ) {
	  for (int kb=ka+1;kb<tempM1.get_size(1);kb++ ) {
	    tempM1(ka,kb) = conj(tempM1(kb,ka));
	  }
	}
	
	
	//GDEBUG(" (%d,%d) = (%d,%d) X (%d,%d) \n", tempM2.get_size(0),tempM2.get_size(1),tempM1.get_size(0),tempM1.get_size(1),psiMatrix.get_size(1),psiMatrix.get_size(0));
	gemm( tempM2, tempM1, false, psiMatrix, true );
	
	
	//GDEBUG(" (%d,%d) = (%d,%d) X (%d,%d) \n", P1.get_size(0),P1.get_size(1),psiMatrix.get_size(0),psiMatrix.get_size(1),tempM2.get_size(0),tempM2.get_size(1));
	gemm( P1, psiMatrix, false, tempM2, false );
	
	subtract(IdentMat,P1,P);
	
	// Keep all projector matrices together
	for( int k1=0;k1<nte;k1++) {
	  for( int k2=0;k2<nte;k2++) {
	    Ps(k1,k2,k3,k4) = P(k1,k2);
	  }
	}
      }
    }
    
    
    // Need to check that S = nte
    // N should be the number of contrasts (eg: for PSIR)
    hoMatrix< std::complex<float> > tempResVector(S,N);
    hoMatrix< std::complex<float> > tempSignal(S,N);
    hoNDArray<float> residual(num_fm,X,Y,Z);
    hoNDArray<uint16_t> r2starIndex(num_fm,X,Y,Z);
    hoNDArray<uint16_t> fmIndex(X,Y,Z);
    float curResidual, minResidual, minResidual2;
    for( int kx=0;kx<X;kx++) {
      for( int ky=0;ky<Y;ky++) {
	for( int kz=0;kz<Z;kz++) {
	  // Get current signal
	  for( int kn=0;kn<N;kn++) {
	    for( int ks=0;ks<S;ks++) {
	      tempSignal(ks,kn) = data(kx,ky,kz,0,kn,ks,0);
	    }
	  }
	  
	  minResidual2 = 1.0 + pow(nrm2(&tempSignal),2.0);
	  
	  for(int kfm=0;kfm<num_fm;kfm++) {
	    
	    minResidual = 1.0 + pow(nrm2(&tempSignal),2.0);
	    
	    for(int kr=0;kr<num_r2star;kr++) {
	      // Get current projector matrix
	      for( int kt1=0;kt1<nte;kt1++) {
		for( int kt2=0;kt2<nte;kt2++) {
		  P(kt1,kt2) = Ps(kt1,kt2,kfm,kr);
		}
	      }
	      
	      // Apply projector
	      gemm( tempResVector, P, false, tempSignal, false );
	      
	      curResidual = pow(nrm2(&tempResVector),2);
	      
	      if (curResidual < minResidual) {
		minResidual = curResidual;
		r2starIndex(kfm,kx,ky,kz) = kr;
	      }
	    }
	    residual(kfm,kx,ky,kz) = minResidual;
	    
	    if (minResidual < minResidual2) {
	      minResidual2 = minResidual;
	      fmIndex(kx,ky,kz) = kfm;
	    }
	  }
	}
      }
    }

    
    float residual_max = Gadgetron::max(&residual);
    float BIG_NUMBER = 100000;
    for(int kfm=0;kfm<num_fm;kfm++) {
      for( int kx=0;kx<X;kx++ ) {
	for( int ky=0;ky<Y;ky++ ) {
	  for( int kz=0;kz<Z;kz++ ) {
	    residual(kfm,kx,ky,kz) = (BIG_NUMBER/residual_max)*residual(kfm,kx,ky,kz);
	  }
	}
      }
    }
    


  
 

    // Initialize field map (indexed on the quantized fm values) to be in the middle of the field map range (ie: all zeroes if symmetric range)
    hoNDArray< uint16_t > cur_ind(X,Y,Z); // field map index
    cur_ind.fill((int)(num_fm/2));
    // If we want to consider only one field map value, no need for any fancy stuff
    if( num_fm>1 ) { // Otherwise, do graph cut iterations


      float delta_fm = fms[2]-fms[1];
      hoNDArray< float > lmap(X,Y,Z); // regularization parameter map
      int fm_min_index;
      int fm_min;
      // Set regularization parameter map for spatially-varying regularization
      for( int kx=0;kx<X;kx++ ) {
	for( int ky=0;ky<Y;ky++ ) {
	  for( int kz=0;kz<Z;kz++ ) {
	  
	    // Find minimum residual at this (kx,ky,kz) pixel
	    fm_min = residual(1,kx,ky,kz);
	    fm_min_index = 1;
	    for( int kfm=1;kfm<num_fm-1;kfm++ ) {
	      if( residual(kfm,kx,ky,kz) < fm_min ) {
		fm_min = residual(kfm,kx,ky,kz);
		fm_min_index = kfm;
	      }
	    }
	    lmap(kx,ky,kz) = abs((residual(fm_min_index+1,kx,ky,kz) + residual(fm_min_index-1,kx,ky,kz) - 2*residual(fm_min_index,kx,ky,kz))/(delta_fm*delta_fm));
	    lmap(kx,ky,kz) = lambda*pow(lmap(kx,ky,kz),lmap_power/2.0);
	  
	  }
	}
      }


      float lmap_mean = Gadgetron::mean(&lmap);
      // Set regularization parameter map for spatially-varying regularization
      for( int kx=0;kx<X;kx++ ) {
	for( int ky=0;ky<Y;ky++ ) {
	  for( int kz=0;kz<Z;kz++ ) {
	    lmap(kx,ky,kz) += lmap_mean*lambda_extra;
	  }
	}
      }    


      for(int ks=0;ks<S;ks++) {
	std::cout << "Voxel 90, 100: signal = " << data(90,100,0,0,0,ks,0) <<  std::endl;
      }
	
	
      for(int kfm=0;kfm<num_fm;kfm++) {
	std::cout << "Voxel 90, 100: fm = " << fms[kfm] << ", residual = " << residual(kfm,90,100,0) << std::endl;
      }


    
      // Form the graph
      uint32_t num_nodes = X*Y*Z; // One node per voxel, num_nodes excludes source and sink
      uint32_t num_edges = num_nodes*(2 + (size_clique+1)^2); // Number of edges, including data and regularization terms


      // Add some graph stuff
     
      using namespace boost;
 
      typedef int EdgeWeightType;
 
      typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
      typedef adjacency_list < vecS, vecS, directedS,
			       property < vertex_name_t, std::string,
					  property < vertex_index_t, long,
						     property < vertex_color_t, boost::default_color_type,
								property < vertex_distance_t, long,
									   property < vertex_predecessor_t, Traits::edge_descriptor > > > > >,
 
			       property < edge_capacity_t, EdgeWeightType,
					  property < edge_residual_capacity_t, EdgeWeightType,
						     property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;
 
      int num_big_jumps = 5;
      for(int kiter=0;kiter<num_iterations;kiter++) {
	
	
	// Graphcut: big jumps
	
	// Find the next candidate index at each voxel (ie: next local minimum)
	hoNDArray< uint16_t > next_ind(X,Y,Z); // field map index
	float cur_sign = pow(-1,kiter);
	uint16_t next_ind_voxel;
	bool found_local_min;
	for( int kx=0;kx<X;kx++ ) {
	  for( int ky=0;ky<Y;ky++ ) {
	    for( int kz=0;kz<Z;kz++ ) {
	      
	      if(kiter<num_big_jumps) {
		// Find next local minimizer of residual at this (kx,ky,kz) pixel
		fm_min = residual(1,kx,ky,kz);
		next_ind_voxel = cur_ind(kx,ky,kz) + (int)cur_sign;
		found_local_min = false;
		while (next_ind_voxel < num_fm-1 && next_ind_voxel > 0 && !found_local_min) {
		  if( residual(next_ind_voxel,kx,ky,kz) < residual(next_ind_voxel-1,kx,ky,kz) &&  residual(next_ind_voxel,kx,ky,kz) < residual(next_ind_voxel+1,kx,ky,kz) ) {
		    found_local_min = true;
		  } else { 
		    next_ind_voxel = next_ind_voxel + (int)cur_sign;
		  }
		}
		
		if( !found_local_min ) {
		  if(cur_sign>0) {
		    next_ind_voxel = num_fm-1;
		  } else {
		    next_ind_voxel = 0;
		  }
		}	
		
		
	      } else {
		
		next_ind_voxel = cur_ind(kx,ky,kz) + (int)(cur_sign); // DH* Need to add bigger jumps here

		if(next_ind_voxel>=num_fm) 
		  next_ind_voxel=num_fm-1;

		if(next_ind_voxel<=0) 
		  next_ind_voxel=0;

	      }
	      next_ind(kx,ky,kz) = next_ind_voxel;

	    }
	  }
	}
      
	std::cout << "Voxel 90, 100: _____ cur ind = " << cur_ind(90,100,0) << ", fm = " << fms[cur_ind(90,100,0)] << ", next ind = " << next_ind(90,100,0)  << ", fm = " << fms[next_ind(90,100,0)] << std::endl;
	
	
	Graph g; //a graph with 0 vertices
	
	property_map < Graph, edge_reverse_t >::type rev = get(edge_reverse, g);
	
	//add a source and sink node, and store them in s and t, respectively
	Traits::vertex_descriptor s = add_vertex(g);
	hoNDArray<Traits::vertex_descriptor> v(X,Y,Z);
	for(int kx=0;kx<X;kx++) {
	  for(int ky=0;ky<Y;ky++) {
	    for(int kz=0;kz<Z;kz++) {
	      v(kx,ky,kz) = add_vertex(g);
	    }
	  }
	}
	Traits::vertex_descriptor t = add_vertex(g);
	


	int min_edge = 10000;
	int max_edge = -10000;

	float dist;
	float a,b,c,d;
	float curlmap;
	for(int kx=0;kx<X;kx++) {
	  for(int ky=0;ky<Y;ky++) {
	    for(int kz=0;kz<Z;kz++) {


	      //float val_sv = std::max(float(0.0),residual(next_ind(kx,ky,kz),kx,ky,kz)-residual(cur_ind(kx,ky,kz),kx,ky,kz));
	      //float val_vt = std::max(float(0.0),residual(cur_ind(kx,ky,kz),kx,ky,kz)-residual(next_ind(kx,ky,kz),kx,ky,kz));
	      //AddEdge(s, v(kx,ky,kz), rev, (int)round(val_sv), g);
	      //AddEdge(v(kx,ky,kz), t, rev, (int)round(val_vt), g);


	      float resNext = residual(next_ind(kx,ky,kz),kx,ky,kz);
	      float resCur = residual(cur_ind(kx,ky,kz),kx,ky,kz);

	      double val_sv = (double)(std::max(float(0.0),resNext-resCur));
	      double val_vt = (double)(std::max(float(0.0),resCur-resNext));

	      AddEdge(s, v(kx,ky,kz), rev, val_sv, g);
	      AddEdge(v(kx,ky,kz), t, rev, val_vt, g);
	      
	      if(kx==90 && ky==100)
		std::cout << " Val_sv = " << val_sv << ", val_vt = " << val_vt << std::endl;

	      if (val_sv > max_edge)
		max_edge = val_sv;
	      if (val_vt > max_edge)
		max_edge = val_vt;

	      if (val_sv < min_edge)
		min_edge = val_sv;
	      if (val_vt < min_edge)
		min_edge = val_vt;
	      



	      //	      if(kx==90 && ky==100)
	      //		std::cout << " Wanted.... Val_sv = " << (int)round(val_sv) << ", val_vt = " << (int)round(val_vt) << std::endl;

	      
	      for(int dx=-size_clique;dx<=size_clique;dx++) {
		for(int dy=-size_clique;dy<=size_clique;dy++) {
		  for(int dz=-size_clique;dz<=size_clique;dz++) {
		    
		    dist = pow(dx*dx+dy*dy+dz*dz,0.5);
		    
		    if(kx+dx>=0 && kx+dx<X && ky+dy>=0 && ky+dy<Y && kz+dz>=0 && kz+dz<Z && dist>0) {
		      
		      curlmap = std::min(lmap(kx,ky,kz),lmap(kx+dx,ky+dy,kz+dz));
		      
		      a = curlmap/dist*pow(cur_ind(kx,ky,kz) - cur_ind(kx+dx,ky+dy,kz+dz),2); 
		      b = curlmap/dist*pow(cur_ind(kx,ky,kz) - next_ind(kx+dx,ky+dy,kz+dz),2); 
		      c = curlmap/dist*pow(next_ind(kx,ky,kz) - cur_ind(kx+dx,ky+dy,kz+dz),2); 
		      d = curlmap/dist*pow(next_ind(kx,ky,kz) - next_ind(kx+dx,ky+dy,kz+dz),2); 
		      
		      // if(kx==109 && ky==149)
		      //   std::cout << " Val_bcad = " << (b+c-a-d) << ", a = " << a << ", b = " << b << ", c = " << c << ", d = " << d << std::endl;

		      AddEdge(v(kx,ky,kz), v(kx+dx,ky+dy,kz+dz), rev, (int)(b+c-a-d), g);
		      AddEdge(s, v(kx,ky,kz), rev, (int)(std::max(float(0.0),c-a)), g);
		      AddEdge(v(kx,ky,kz), t, rev, (int)(std::max(float(0.0),a-c)), g);
		      AddEdge(s, v(kx+dx,ky+dy,kz+dz), rev, (int)(std::max(float(0.0),d-c)), g);
		      AddEdge(v(kx+dx,ky+dy,kz+dz), t, rev, (int)(std::max(float(0.0),c-d)), g);
		      
		    }
		    
		  }
		}
	      }

	    }
	  }
	}


	std::cout << " MIN EDGE = " << min_edge << ", MAX_EDGE = " << max_edge << std::endl;
		
		
	//EdgeWeightType flow = edmonds_karp_max_flow(g, s, t); // a list of sources will be returned in s, and a list of sinks will be returned in t
	//EdgeWeightType flow = push_relabel_max_flow(g, s, t); // a list of sources will be returned in s, and a list of sinks will be returned in t
	EdgeWeightType flow = boykov_kolmogorov_max_flow(g, s,t); // a list of sources will be returned in s, and a list of sinks will be returned in t



	std::cout << "Max flow is: " << flow << std::endl;

	property_map<Graph, edge_capacity_t>::type
	  capacity = get(edge_capacity, g);

	property_map<Graph, edge_residual_capacity_t>::type
	  residual_capacity = get(edge_residual_capacity, g);

	property_map<Graph, vertex_color_t>::type
	  colormap = get(vertex_color, g);


	//	std::cout << "c flow values:" << std::endl;
	//	graph_traits<Graph>::vertex_iterator u_iter, u_end;
	//	graph_traits<Graph>::out_edge_iterator ei, e_end;
	//	for (tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
	//	  for (tie(ei, e_end) = out_edges(*u_iter, g); ei != e_end; ++ei)
	//	    if (capacity[*ei] > 0)
	//	      std::cout << "f " << *u_iter << " " << target(*ei, g) << " " << (capacity[*ei]) << " " << (residual_capacity[*ei]) << " " << (capacity[*ei] - residual_capacity[*ei]) << std::endl;
				


	for( int kx=0;kx<X;kx++) {
	  for( int ky=0;ky<Y;ky++) {
	    for( int kz=0;kz<Z;kz++) {

	      //	      if(colormap[1 + kx + ky*X + kz*X*Y] > 0)
	      //		std::cout << " ColorMap = " << colormap[1 + kx + ky*X + kz*X*Y] << ", x = " << kx << ", y = " << ky << std::endl;

	      if(kx==90 && ky==100) {
		std::cout << " ColorMap = " << colormap[1 + ky + kx*Y + kz*X*Y] << ", Colormap0 = " << colormap[0] << ", ColormapEnd = " << colormap[X*Y*Z+1] << std::endl;
		std::cout << " posVoxel = " << 1 + ky + kx*Y + kz*X*Y  << ", posEnd = " << X*Y*Z+1 << std::endl;
	      }


	      if(colormap[1 + ky + kx*Y + kz*X*Y]!=colormap[0])
		cur_ind(kx,ky,kz) = next_ind(kx,ky,kz);
	    }
	  }
	}
    
	// Graphcut: small jumps
      }
      
      
    } // End else (all the graph cut iteration portion)


    std::cout << " Final estimate at voxel [90,100]... fm = " << fms[cur_ind(90,100,0)] << std::endl; 
    
    //Do final calculations once the field map is done
    hoMatrix< std::complex<float> > curWaterFat(2,N);
    hoMatrix< std::complex<float> > AhA(2,2);
    // Do fat-water separation with current field map and R2* estimates
    for( int kx=0;kx<X;kx++) {
      for( int ky=0;ky<Y;ky++) {
	for( int kz=0;kz<Z;kz++) {
	  // Get current signal
	  for( int kn=0;kn<N;kn++) {
	    for( int ks=0;ks<S;ks++) {
	      tempSignal(ks,kn) = data(kx,ky,kz,0,kn,ks,0);
	    }
	  }
	  // Get current Psi matrix
	  fm = fms[cur_ind(kx,ky,kz)];
	  r2star =r2stars[r2starIndex(cur_ind(kx,ky,kz),kx,ky,kz)];
	  for( int kt=0;kt<nte;kt++) {
	    curModulation = exp(-r2star*echoTimes[kt])*std::complex<float>(cos(2*PI*echoTimes[kt]*fm),sin(2*PI*echoTimes[kt]*fm));
	    for( int ksp=0;ksp<nspecies;ksp++) {
	      psiMatrix(kt,ksp) = phiMatrix(kt,ksp)*curModulation;
	    }
	  }
	  
	  // Solve for water and fat
	  gemm( curWaterFat, psiMatrix, true, tempSignal, false );
	  herk( AhA, psiMatrix, 'L', true );
	  //	    AhA.copyLowerTriToUpper();
	  for (int ka=0;ka<AhA.get_size(0);ka++ ) {
	    for (int kb=ka+1;kb<AhA.get_size(1);kb++ ) {
	      AhA(ka,kb) = conj(AhA(kb,ka));
	    }
	  }
	  
	  hesv(AhA,curWaterFat);
	  for ( int kn=0;kn<N;kn++ ) {
	    for ( int ks=0;ks<nspecies;ks++ ) { // 2 elements for water and fat currently
	      out(kx,ky,kz,0,kn,ks,0) = curWaterFat(ks,kn);
	    }
	  }
	  

	  //	  out(kx,ky,kz,0,0,0,0) = fm;
	  //	  out(kx,ky,kz,0,0,1,0) = r2star;


	}
      }
    }    
    std::cout << mytimer.format() << '\n';  
    
    //Clean up as needed
    
    
    return out;
  }
}



