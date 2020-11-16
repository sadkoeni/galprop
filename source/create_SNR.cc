
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * create_SNR.cc *                               galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

using namespace std;//AWS20050624
#include<cstdlib>
#include"galprop_classes.h"
#include"galprop_internal.h"
#include"constants.h" //for Rsun
#include"SourceClass_Compatibility.h"


void SourceClass_Compatibility::create_SNR(const Particle &particle) const
{
   INFO("Entry");

   if (particle.n_spatial_dimensions != 3 && SNR_events != 1)
      return;

   if (! createSNRDistributions ) {
      //Make sure the dimensions haven't changed
      if ( particle.n_xgrid == SNR_cell_time.n_xgrid &&
            particle.n_ygrid == SNR_cell_time.n_ygrid &&
            particle.n_zgrid == SNR_cell_time.n_zgrid )
         return;
   }

   SNR_cell_time.init(particle.n_xgrid, particle.n_ygrid, particle.n_zgrid, 1);
   SNR_cell_phase.init(particle.n_xgrid, particle.n_ygrid, particle.n_zgrid, 1);

   SNR_electron_dg.init(particle.n_xgrid, particle.n_ygrid, particle.n_zgrid, 1); //AWS20010410
   SNR_nuc_dg.init(particle.n_xgrid, particle.n_ygrid, particle.n_zgrid, 1); //AWS20010410

   const double cell_volume=particle.dx * particle.dy * particle.dz;

   const unsigned seed=1234;
   srand(seed); // eventually use a parameter

   const double localDensity = source_distribution(Rsun, 0, 0, source_model, particle.n_spatial_dimensions, source_parameters, source_values, source_radius);

   for(int ix=0;ix<particle.n_xgrid;ix++){
      for(int iy=0;iy<particle.n_ygrid;iy++){
         for(int iz=0;iz<particle.n_zgrid;iz++){
  
  
            SNR_cell_time .d3[ix][iy][iz].s[0] = SNR_interval * localDensity/
               (source_distribution(particle.x[ix], particle.y[iy], particle.z[iz], source_model, particle.n_spatial_dimensions, source_parameters, source_values, source_radius)+1.e-30)/cell_volume;

     
            SNR_cell_phase.d3[ix][iy][iz].s[0]=float(rand())/RAND_MAX;


            // Gaussian distributed source index delta                                       AWS20010410
            SNR_electron_dg.d3[ix][iy][iz].s[0]=gauss(0.0,SNR_electron_sdg); //AWS20010410
            SNR_nuc_dg     .d3[ix][iy][iz].s[0]=gauss(0.0,SNR_nuc_sdg);      //AWS20010410

         }
      }
   }
 
   createSNRDistributions = false;
}
