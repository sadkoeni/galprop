
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Galaxy.cc *                                   galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include<iostream>
#include"Galaxy.h"
#include"ErrorLogger.h"
#include <FullSky.h>
#include <SparseSky.h>

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

//using namespace rf;

Galaxy::Galaxy() : n_Ring(0), fISRFloadedfile("") {

   init_pointers();

}

Galaxy::Galaxy(const Galaxy & old) {
   //Copy all the variables
   copy_variables(old);
   //Do a deep copy of all the pointers that are initialized in the old version
   copy_pointers_deep(old);
}

void Galaxy::init_pointers() {

  fISRFrf = 0;

  ISRF = ISRF_energy_density = IC_iso_emiss = IC_iso_skymap = IC_aniso_skymap = 0;

  IC_iso_rings_skymap = IC_aniso_rings_skymap = 0;
  
}

void Galaxy::copy_variables(const Galaxy& old) {

  z_min = old.z_min;
  z_max = old.z_max;
  dz = old.dz;
  r_min = old.r_min;
  r_max = old.r_max;
  dr = old.dr;
  x_min = old.x_min;
  x_max = old.x_max;
  dx = old.dx;
  y_min = old.y_min;
  y_max = old.y_max;
  dy = old.dy;
  p_min = old.p_min;
  p_max = old.p_max;
  p_factor = old.p_factor;
  
  n_spatial_dimensions = old.n_spatial_dimensions;
  
  n_pgrid = old.n_pgrid;
  n_xgrid = old.n_xgrid;
  n_rgrid = old.n_rgrid;
  n_zgrid = old.n_zgrid;
  n_ygrid = old.n_ygrid;

  p.resize(old.p.size());
  p = old.p;
  x.resize(old.x.size());
  x = old.x;
  y.resize(old.y.size());
  y = old.y;
  z.resize(old.z.size());
  z = old.z;
  r.resize(old.r.size());
  r = old.r;
  
  n_E_gammagrid = old.n_E_gammagrid;
  E_gamma_min = old.E_gamma_min;
  E_gamma_max = old.E_gamma_max;
  E_gamma_factor = old.E_gamma_factor;
  E_gamma.resize(old.E_gamma.size());
  E_gamma = old.E_gamma;
  
  n_lat = old.n_lat;
  n_long = old.n_long;
  long_min = old.long_min;
  long_max = old.long_max;
  lat_min = old.lat_min;
  lat_max = old.lat_max;
  d_long = old.d_long;
  d_lat = old.d_lat;
  
  n_nu_synchgrid = old.n_nu_synchgrid;
  nu_synch_min = old.nu_synch_min;
  nu_synch_max = old.nu_synch_max;
  nu_synch_factor = old.nu_synch_factor;
  nu_synch.resize(old.nu_synch.size());
  nu_synch = old.nu_synch;
  
  n_H2 = old.n_H2;
  n_HI = old.n_HI;
  n_HII = old.n_HII;
  HIR = old.HIR;
  COR = old.COR;
  n_Ring = old.n_Ring;
  B_field = old.B_field;
  B_field_regular = old.B_field_regular;
  B_field_random = old.B_field_random;  
  n_ISRF_components = old.n_ISRF_components;
  fISRFFactors.resize(old.fISRFFactors.size());
  fISRFFactors = old.fISRFFactors;
  nu_ISRF.resize(old.nu_ISRF.size());
  nu_ISRF = old.nu_ISRF;
  bremss_emiss = old.bremss_emiss;
  bremss_ionized_emiss = old.bremss_ionized_emiss;
  pi0_decay_emiss = old.pi0_decay_emiss;
  bremss_skymap = old.bremss_skymap;
  bremss_ionized_skymap = old.bremss_ionized_skymap;
  pi0_decay_skymap = old.pi0_decay_skymap;
  bremss_H2R_skymap = old.bremss_H2R_skymap;
  bremss_HIR_skymap = old.bremss_HIR_skymap;
  bremss_HIIR_skymap = old.bremss_HIIR_skymap;
  pi0_decay_H2R_skymap = old.pi0_decay_H2R_skymap;
  pi0_decay_HIR_skymap = old.pi0_decay_HIR_skymap;
  pi0_decay_HIIR_skymap = old.pi0_decay_HIIR_skymap;
  synchrotron_emiss = old.synchrotron_emiss;
  synchrotron_skymap = old.synchrotron_skymap;
  ionization_rate = old.ionization_rate;
  DM_emiss = old.DM_emiss;
  DM_skymap = old.DM_skymap;

  fTotalISRFNumberDensity = old.fTotalISRFNumberDensity;
  
  R_bins.resize(old.R_bins.size());
  R_bins = old.R_bins;

  fProtonLuminosity.resize(old.fProtonLuminosity.size());
  fProtonLuminosity = old.fProtonLuminosity;

  fHeliumLuminosity.resize(old.fHeliumLuminosity.size());
  fHeliumLuminosity = old.fHeliumLuminosity;

  fNucleiLuminosity.resize(old.fNucleiLuminosity.size());
  fNucleiLuminosity = old.fNucleiLuminosity;

  fPrimaryElectronLuminosity.resize(old.fPrimaryElectronLuminosity.size());
  fPrimaryElectronLuminosity = old.fPrimaryElectronLuminosity;

  fSecondaryElectronLuminosity.resize(old.fSecondaryElectronLuminosity.size());
  fSecondaryElectronLuminosity = old.fSecondaryElectronLuminosity;

  fSynchrotronLuminosity.resize(old.fSynchrotronLuminosity.size());
  fSynchrotronLuminosity = old.fSynchrotronLuminosity;

  fICLuminosity.resize(old.fICLuminosity.size());
  fICLuminosity = old.fICLuminosity;

  fPi0DecayLuminosity.resize(old.fPi0DecayLuminosity.size());
  fPi0DecayLuminosity = old.fPi0DecayLuminosity;

  fBremsstrahlungLuminosity.resize(old.fBremsstrahlungLuminosity.size());
  fBremsstrahlungLuminosity = old.fBremsstrahlungLuminosity;

  fDMLuminosity.resize(old.fDMLuminosity.size());
  fDMLuminosity = old.fDMLuminosity;

  fISRFLuminosity.resize(old.fISRFLuminosity.size());
  fISRFLuminosity = old.fISRFLuminosity;

  //Note that we specifically don't copy the Radiation field since it does not have a copy constructor
  //The same goes for the ISRF filename
  
}

void Galaxy::copy_pointers_deep(const Galaxy& old) {

  //Does not free allocated memory so this can be used in copy constructor
  init_pointers(); //Don't want un-initialized pointers.
  if (old.ISRF) {
    ISRF = new Distribution[old.n_ISRF_components];
    for (int i = 0; i < old.n_ISRF_components; ++i) {
      ISRF[i] = old.ISRF[i];
    }
  }
  if (old.ISRF_energy_density) {
    ISRF_energy_density = new Distribution[old.n_ISRF_components];
    for (int i = 0; i < old.n_ISRF_components; ++i) {
      ISRF_energy_density[i] = old.ISRF_energy_density[i];
    }
  }
  if (old.IC_iso_emiss) {
    IC_iso_emiss = new Distribution[old.n_ISRF_components];
    for (int i = 0; i < old.n_ISRF_components; ++i) {
      IC_iso_emiss[i] = old.IC_iso_emiss[i];
    }
  }

  if (old.IC_iso_skymap) {
    IC_iso_skymap = new Distribution[old.n_ISRF_components];
    for (int i = 0; i < old.n_ISRF_components; ++i) {
      IC_iso_skymap[i] = old.IC_iso_skymap[i];
    }
  }
  if (old.IC_aniso_skymap) {
    IC_aniso_skymap = new Distribution[old.n_ISRF_components];
      for (int i = 0; i < old.n_ISRF_components; ++i) {
	IC_aniso_skymap[i] = old.IC_aniso_skymap[i];
      }
  }

  //Copy the BaseSky pointers
  if (old.hpCOR.size() > 0) {
     hpCOR.resize(old.hpCOR.size());
     for (size_t i = 0; i < old.hpCOR.size(); ++i)
        hpCOR[i] = old.hpCOR[i]->clone();
  }
  if (old.hpHIR.size() > 0) {
     hpHIR.resize(old.hpHIR.size());
     for (size_t i = 0; i < old.hpHIR.size(); ++i)
        hpHIR[i] = old.hpHIR[i]->clone();
  }

  if (old.bremss_hp_skymap.get() != nullptr)
     bremss_hp_skymap = old.bremss_hp_skymap->clone();
  if (old.bremss_ionized_hp_skymap.get() != nullptr)
     bremss_ionized_hp_skymap = old.bremss_ionized_hp_skymap->clone();
  if (old.pi0_decay_hp_skymap.get() != nullptr)
     pi0_decay_hp_skymap = old.pi0_decay_hp_skymap->clone();

  if (old.synchrotron_hp_skymap.get() != nullptr)
     synchrotron_hp_skymap = old.synchrotron_hp_skymap->clone();
  if (old.synchrotron_Q_hp_skymap.get() != nullptr)
     synchrotron_Q_hp_skymap = old.synchrotron_Q_hp_skymap->clone();
  if (old.synchrotron_U_hp_skymap.get() != nullptr)
     synchrotron_U_hp_skymap = old.synchrotron_U_hp_skymap->clone();
  if (old.synchrotron_P_hp_skymap.get() != nullptr)
     synchrotron_P_hp_skymap = old.synchrotron_P_hp_skymap->clone();
  if (old.synchrotron_polang_hp_skymap.get() != nullptr)
     synchrotron_polang_hp_skymap = old.synchrotron_polang_hp_skymap->clone();
  if (old.synchrotron_polfra_hp_skymap.get() != nullptr)
     synchrotron_polfra_hp_skymap = old.synchrotron_polfra_hp_skymap->clone();
  if (old.free_free_hp_skymap.get() != nullptr)
     free_free_hp_skymap = old.free_free_hp_skymap->clone();

  if (old.DM_hp_skymap.get() != nullptr)
     DM_hp_skymap = old.DM_hp_skymap->clone();

  if (old.IC_iso_hp_skymap.size() > 0) {
    IC_iso_hp_skymap.resize(old.n_ISRF_components);
    for (int i = 0; i < old.n_ISRF_components; ++i) {
      IC_iso_hp_skymap[i] = old.IC_iso_hp_skymap[i]->clone();
    }
  }
  if (old.IC_aniso_hp_skymap.size() > 0) {
    IC_aniso_hp_skymap.resize(old.n_ISRF_components);
    for (int i = 0; i < old.n_ISRF_components; ++i) {
      IC_aniso_hp_skymap[i] = old.IC_aniso_hp_skymap[i]->clone();
    }
  }
  if (old.IC_iso_rings_hp_skymap.size() > 0) {
    IC_iso_rings_hp_skymap.resize(old.n_ISRF_components);
    for (int i = 0; i < old.n_ISRF_components; ++i) {
      IC_iso_rings_hp_skymap[i].resize(old.n_Ring);
      for (int iR(0); iR < old.n_Ring; ++iR) {
         IC_iso_rings_hp_skymap[i][iR] = old.IC_iso_rings_hp_skymap[i][iR]->clone();
      }
    }
  }
  if (old.IC_aniso_rings_hp_skymap.size() > 0) {
    IC_aniso_rings_hp_skymap.resize(old.n_ISRF_components);
    for (int i = 0; i < old.n_ISRF_components; ++i) {
      IC_aniso_rings_hp_skymap[i].resize(old.n_Ring);
      for (int iR(0); iR < old.n_Ring; ++iR) {
         IC_aniso_rings_hp_skymap[i][iR] = old.IC_iso_rings_hp_skymap[i][iR]->clone();
      }
    }
  }
  if (old.IC_iso_rings_skymap) {
     IC_iso_rings_skymap = new Distribution[old.n_ISRF_components];
     for (int i = 0; i < old.n_Ring; ++i) {
        IC_iso_rings_skymap[i] = old.IC_iso_rings_skymap[i];
     }
  }
  if (old.IC_aniso_rings_skymap) {
     IC_aniso_rings_skymap = new Distribution[old.n_ISRF_components];
     for (int i = 0; i < old.n_Ring; ++i) {
        IC_aniso_rings_skymap[i] = old.IC_aniso_rings_skymap[i];
     }
  }
  if (old.bremss_H2R_hp_skymap.size() > 0){
    bremss_H2R_hp_skymap.resize(old.n_Ring);
    for (int i = 0; i < old.n_Ring; ++i) {
      bremss_H2R_hp_skymap[i] = old.bremss_H2R_hp_skymap[i]->clone();
    }
  }
  if (old.bremss_HIR_hp_skymap.size() > 0){
    bremss_HIR_hp_skymap.resize(old.n_Ring);
    for (int i = 0; i < old.n_Ring; ++i) {
      bremss_HIR_hp_skymap[i] = old.bremss_HIR_hp_skymap[i]->clone();
    }
  }
  if (old.bremss_HIIR_hp_skymap.size() > 0){
    bremss_HIIR_hp_skymap.resize(old.n_Ring);
    for (int i = 0; i < old.n_Ring; ++i) {
      bremss_HIIR_hp_skymap[i] = old.bremss_HIIR_hp_skymap[i]->clone();
    }
  }
  if (old.pi0_decay_H2R_hp_skymap.size() > 0){
    pi0_decay_H2R_hp_skymap.resize(old.n_Ring);
    for (int i = 0; i < old.n_Ring; ++i) {
      pi0_decay_H2R_hp_skymap[i] = old.pi0_decay_H2R_hp_skymap[i]->clone();
    }
  }
  if (old.pi0_decay_HIR_hp_skymap.size() > 0){
    pi0_decay_HIR_hp_skymap.resize(old.n_Ring);
    for (int i = 0; i < old.n_Ring; ++i) {
      pi0_decay_HIR_hp_skymap[i] = old.pi0_decay_HIR_hp_skymap[i]->clone();
    }
  }
  if (old.pi0_decay_HIIR_hp_skymap.size() > 0){
    pi0_decay_HIIR_hp_skymap.resize(old.n_Ring);
    for (int i = 0; i < old.n_Ring; ++i) {
      pi0_decay_HIIR_hp_skymap[i] = old.pi0_decay_HIIR_hp_skymap[i]->clone();
    }
  }

  //Note that we specifically don't copy the Radiation field since it does not have a copy constructor
  //The same goes for the filename
  
}

Galaxy::~Galaxy() {

  free_memory();

}

void Galaxy::free_memory() {

  delete[] ISRF;
  delete[] ISRF_energy_density;
  
  delete[] IC_iso_emiss;
  //delete IC_aniso_emiss;
  delete[] IC_iso_skymap;
  delete[] IC_aniso_skymap;
  delete[] IC_iso_rings_skymap;
  delete[] IC_aniso_rings_skymap;
  
  delete fISRFrf;
}

void Galaxy::init(double r_min_, double r_max_, double dr_, 
                  double z_min_, double z_max_, double dz_) {

  INFO("Initializing 2D");
  n_spatial_dimensions=2;       //2D    
  
  r_min=r_min_;
  r_max=r_max_;
  dr   =   dr_;
  z_min=z_min_;
  z_max=z_max_;
  dz   =   dz_;
  
  n_rgrid=(int)((r_max-r_min)/dr + 1.5);
  n_zgrid=(int)((z_max-z_min)/dz + 1.5);
  
  r.resize(n_rgrid);
  z.resize(n_zgrid);

  int ir,iz;    
  for(ir=0; ir<n_rgrid; ir++) r[ir]=r_min+ir*dr;
  for(iz=0; iz<n_zgrid; iz++) z[iz]=z_min+iz*dz;
  
  n_HI.init(n_rgrid, n_zgrid, 1);
  n_H2.init(n_rgrid, n_zgrid, 1);
  n_HII.init(n_rgrid, n_zgrid, 1);
  B_field.init(n_rgrid, n_zgrid, 1);
  B_field_regular.init(n_rgrid, n_zgrid, 1);
  B_field_random.init(n_rgrid, n_zgrid, 1);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void Galaxy::init( double x_min_, double x_max_, double dx_, 
                   double y_min_, double y_max_, double dy_, 
                   double z_min_, double z_max_, double dz_) {  

  INFO("Initializing 3D");
  n_spatial_dimensions=3;       //3D    
  
  x_min=x_min_;
  x_max=x_max_;
  dx   =   dx_;
  y_min=y_min_;
  y_max=y_max_;
  dy   =   dy_;
  z_min=z_min_;
  z_max=z_max_;
  dz   =   dz_;
  
  n_xgrid=(int)((x_max-x_min)/dx + 1.5);
  n_ygrid=(int)((y_max-y_min)/dy + 1.5);
  n_zgrid=(int)((z_max-z_min)/dz + 1.5);
  
  x.resize(n_xgrid);
  y.resize(n_ygrid);
  z.resize(n_zgrid);

  int ix,iy,iz;    
  for(ix=0; ix<n_xgrid; ix++) x[ix]=x_min+ix*dx; 
  for(iy=0; iy<n_ygrid; iy++) y[iy]=y_min+iy*dy; 
  for(iz=0; iz<n_zgrid; iz++) z[iz]=z_min+iz*dz; 
  
  n_HI.init(n_xgrid, n_ygrid, n_zgrid, 1);
  n_H2.init(n_xgrid, n_ygrid, n_zgrid, 1);
  n_HII.init(n_xgrid, n_ygrid, n_zgrid, 1);	
  B_field.init(n_xgrid, n_ygrid, n_zgrid, 1);
  B_field_regular.init(n_xgrid, n_ygrid, n_zgrid, 1);
  B_field_random.init(n_xgrid, n_ygrid, n_zgrid, 1);
  
  
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void Galaxy::print() {

  cout<<"Galaxy: n_spatial_dimensions="<<n_spatial_dimensions<<endl;
  cout<<"        n_ISRF_components   ="<<n_ISRF_components   <<endl;
  
  if(n_spatial_dimensions==2) {

    cout<<"r_min="<<r_min<<endl;
    cout<<"r_max="<<r_max<<endl;
    cout<<"dr   ="<<dr   <<endl;
    cout<<"z_min="<<z_min<<endl;
    cout<<"z_max="<<z_max<<endl;
    cout<<"dz   ="<<dz   <<endl;
    
    cout<<"n_rgrid="<<n_rgrid<<endl;
    cout<<"n_zgrid="<<n_zgrid<<endl;
    
    cout<<"r grid:"<<endl;
    for(int ir=0; ir<n_rgrid; cout<<r[ir++]<<" "); cout<<endl;
    cout<<"z grid:"<<endl;
    for(int iz=0; iz<n_zgrid; cout<<z[iz++]<<" "); cout<<endl;
    
    cout<<"n_HI[0][0].s[0]="<<n_HI.d2[0][0].s[0]<<endl;
  }//(n_spatial_dimensions==2

  if(n_spatial_dimensions==3) {

    cout<<"x_min="<<x_min<<endl;
    cout<<"x_max="<<x_max<<endl;
    cout<<"dx   ="<<dx   <<endl;
    cout<<"y_min="<<y_min<<endl;
    cout<<"y_max="<<y_max<<endl;
    cout<<"dy   ="<<dy   <<endl;
    cout<<"z_min="<<z_min<<endl;
    cout<<"z_max="<<z_max<<endl;
    cout<<"dz   ="<<dz   <<endl;
    
    cout<<"n_xgrid="<<n_xgrid<<endl;
    cout<<"n_ygrid="<<n_ygrid<<endl;
    cout<<"n_zgrid="<<n_zgrid<<endl;
    
    cout<<"x grid:"<<endl;
    for(int ix=0; ix<n_xgrid; cout<<x[ix++]<<" ");  cout<<endl;
    cout<<"y grid:"<<endl;
    for(int iy=0; iy<n_ygrid; cout<<y[iy++]<<" ");  cout<<endl;
    cout<<"z grid:"<<endl;
    for(int iz=0; iz<n_zgrid; cout<<z[iz++]<<" ");  cout<<endl;
 
    cout<<"n_HI[0][0][0].s[0]="<<n_HI.d3[0][0][0].s[0]<<endl;
    cout<<"n_H2[0][0][0].s[0]="<<n_H2.d3[0][0][0].s[0]<<endl;
  
  } //(n_spatial_dimensions==3

}

Galaxy & Galaxy::operator = (const Galaxy & old) {
  //Avoid self-assignment
  if (this != &old) {
    free_memory();
    copy_variables(old);
    copy_pointers_deep(old);
  }
  return (*this);
}

std::unique_ptr< SM::BaseSky<double> > Galaxy::createFullSky(size_t order) const
{
   return std::unique_ptr< SM::BaseSky<double> > ( new SM::FullSky<double> (
            SM::SpectralBinning(std::vector<double>(std::begin(E_gamma), std::end(E_gamma))),
            order, RING, SM::CoordSys::GAL, 0.0));
}

std::unique_ptr< SM::BaseSky<double> > Galaxy::createSparseSky(size_t order) const
{
   return std::unique_ptr< SM::BaseSky<double> > ( new SM::SparseSky<double> (
            SM::SpectralBinning(std::vector<double>(std::begin(E_gamma), std::end(E_gamma))),
            order, RING, SM::CoordSys::GAL, 0.0));
}

