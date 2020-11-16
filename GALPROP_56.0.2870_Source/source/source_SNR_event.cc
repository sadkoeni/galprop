//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * source_SNR_event.cc *                         galprop package * 10/12/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include"galprop_classes.h"
#include"galprop_internal.h"
#include"constants.h"  //For Rsun

#include"SourceClass_Compatibility.h"

#include <ErrorLogger.h>

#include <cstring>

using namespace std;//AWS20050624

void SourceClass_Compatibility::source_SNR_event(Particle& particle, double t) const {

  INFO("Entry");

  //Never called unless 3D
  double p1, p2;
  
  //To allow for different source distribution of electrons and nuclei
  const std::vector<double>* parameters;
  const std::vector<double>* values;
  const std::vector<double>* radius;
  int model;
  
  const string priElecStr = "primary_electrons";

  if (priElecStr != particle.name) {
    
    model = source_model;
    parameters = &source_parameters;
    values = &source_values;
    radius = &source_radius;
    
  } else {
    
    model = source_model_electron;
    parameters = &source_parameters_electron;
    values = &source_values_electron;
    radius = &source_radius_electron;
    
  }

  //This cannot be done.  We must subtract the steady state distribution from this class.
  //particle.primary_source_function=0.0;
  
  //Use the convenience function to get the spectral shape
  auto spec_shape = createSpecShape(particle);

  const double cell_volume=particle.dx * particle.dy * particle.dz;
  const double factor= SNR_interval/SNR_livetime*
     source_distribution(Rsun,0,0, model,particle.n_spatial_dimensions,*parameters,*values,*radius)/cell_volume;
  
  int n_SNR=0;
  
  if (removeSteadyState) {
     for (int ix=0;ix<particle.n_xgrid;ix++){
        for (int iy=0;iy<particle.n_ygrid;iy++){
           for (int iz=0;iz<particle.n_zgrid;iz++){

              const double sd = source_distribution(particle.x[ix],particle.y[iy],particle.z[iz], model,particle.n_spatial_dimensions,*parameters,*values,*radius);
              for (int ip = 0; ip < particle.n_pgrid; ++ip) {

                 double spec_dg_ratio;                             //AWS20010410

                 if ("primary_electrons" == particle.name) { 

                    spec_dg_ratio=
                       pow(particle.rigidity[ip]/SNR_electron_dgpivot, 1.*SNR_electron_dg.d3[ix][iy][iz].s[0]);

                    /*
                       if(galdef.verbose==-501)// selectable debug
                       cout<<"SNR_electron_dg="<<galaxy.SNR_electron_dg.d3[ix][iy][iz].s[0]
                       <<" rigidity="<<particle.rigidity[ip]<<"  spec_dg_ratio="<<spec_dg_ratio<<endl;
                       */
                 }

                 if ("primary_electrons" != particle.name) { 
                    spec_dg_ratio=
                       pow(particle.rigidity[ip]/SNR_nuc_dgpivot, 1.*SNR_nuc_dg.d3[ix][iy][iz].s[0]);

                    /*
                       if(galdef.verbose==-501)// selectable debug
                       cout<<"SNR_nuc_dg="<<galaxy.SNR_nuc_dg.d3[ix][iy][iz].s[0] 
                       <<" rigidity="<<particle.rigidity[ip]<<"  spec_dg_ratio="<<spec_dg_ratio<<endl;
                       */
                 }

                 particle.primary_source_function.d3[ix][iy][iz].s[ip]-=sd*spec_shape[ip]*spec_dg_ratio;   ;

              }
           }
        }
     }
     removeSteadyState = false;
  }

  for (int ix=0;ix<particle.n_xgrid;ix++){
    for (int iy=0;iy<particle.n_ygrid;iy++){
      for (int iz=0;iz<particle.n_zgrid;iz++){
		
	p1 = fmod(t/SNR_cell_time.d3[ix][iy][iz].s[0], 1.);
	p2 = fmod((t + SNR_livetime)/SNR_cell_time.d3[ix][iy][iz].s[0], 1.);
		
	if (p1 < SNR_cell_phase.d3[ix][iy][iz].s[0] && 
	    p2 > SNR_cell_phase.d3[ix][iy][iz].s[0]) {
	  	  
          for (int ip = 0; ip < particle.n_pgrid; ++ip) {
	    
	    double spec_dg_ratio;                             //AWS20010410
	    
	    if ("primary_electrons" == particle.name) { 
	    
	      spec_dg_ratio=
		pow(particle.rigidity[ip]/SNR_electron_dgpivot, 1.*SNR_electron_dg.d3[ix][iy][iz].s[0]);

              /*
	      if(galdef.verbose==-501)// selectable debug
		cout<<"SNR_electron_dg="<<galaxy.SNR_electron_dg.d3[ix][iy][iz].s[0]
		    <<" rigidity="<<particle.rigidity[ip]<<"  spec_dg_ratio="<<spec_dg_ratio<<endl;
                    */
	    }

	    if ("primary_electrons" != particle.name) { 
	      spec_dg_ratio=
		pow(particle.rigidity[ip]/SNR_nuc_dgpivot, 1.*SNR_nuc_dg.d3[ix][iy][iz].s[0]);

              /*
	      if(galdef.verbose==-501)// selectable debug
		cout<<"SNR_nuc_dg="<<galaxy.SNR_nuc_dg.d3[ix][iy][iz].s[0] 
		    <<" rigidity="<<particle.rigidity[ip]<<"  spec_dg_ratio="<<spec_dg_ratio<<endl;
                    */
	    }
	    
	    particle.primary_source_function.d3[ix][iy][iz].s[ip]+=factor*spec_shape[ip]*spec_dg_ratio;   ;

	  }//ip

	  n_SNR += 1;

	  //cout<<" source_SNR_event at t="<<t<<" x y z= "<<galaxy.x[ix]<<" "<<galaxy.y[iy]<<" "<<galaxy.z[iz]<<" SNR_cell_time "<<galaxy.SNR_cell_time .d3[ix][iy][iz].s[0]<<" p1 SNR_cell_phase p2 "<<p1<<" "<<galaxy.SNR_cell_phase.d3[ix][iy][iz].s[0]<<" "<<p2<<" primary_source_function="<<particle.primary_source_function.d3[ix][iy][iz].s[0]<<endl; //AWS20001121
     
	}
	
      }//iz
   
    }//iy
   
  }//ix

  std::ostringstream os;
  os<<"source_SNR_event: number of live SNR at this time = "<<n_SNR<<endl;
  INFO(os.str());
  
  /*
  if(galdef.verbose>=1){
    
    cout<<"source_SNR_event: primary source function for particle "<<particle.name<<endl;
    particle.primary_source_function .print();
    
  }
  */

  INFO("Exit");

}
