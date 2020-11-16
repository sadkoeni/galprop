
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * electrons_normalize.cc *                      galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"
#include <cstring>

#include "ErrorLogger.h"
#include <sstream>

int Galprop::electrons_normalize() {

  INFO("Entry");
  
  // identify the primary electrons 
  int ielectrons=-1;
  for(int i=0;i<n_species;i++) 
    if ("primary_electrons" == gcr[i].name) 
      //strcmp(gcr[i].name,"primary_electrons")==0) 
      ielectrons=i;
  
  if(ielectrons==-1){WARNING("primary electrons not found!"); return 1;}
  ostringstream buf;
  buf<<"  primary electrons found as species #"<<ielectrons;
  INFO(buf.str());
  
  double vnorm = 0;//AWS20101202
  
  if (1 == galdef.electron_norm_type) { //AWS20101202
    
    const double r0 = 8.5; // solar Galactocentric radius, kpc
    // deriving the normalization grid point  
    const int ip = int(log(galdef.electron_norm_Ekin/galdef.Ekin_min)/log(galdef.Ekin_factor) + 0.5); //IMOS20060420
    const int iz = int((0.-galdef.z_min)/galdef.dz + 0.5); // z=0, Galactic plane IMOS20060420

    // Do nothing if the grid points are out of bounds
    if (ip < 0 || ip >= gcr[ielectrons].n_pgrid) {
       ERROR("Electron normalization energy out of bounds, aborting normalization");
       return 1;
    }
    if (iz < 0 || iz >= gcr[ielectrons].n_zgrid) {
       ERROR("Galactic plane not within grid, aborting normalization");
       return 1;
    }
    
    double v1(0),v2(0),v3(0),v4(0),v5(0),v6(0);
    if (2 == galdef.n_spatial_dimensions) {

      const int ir = int((r0-galdef.r_min)/galdef.dr + 0.5);//IMOS20060420

      if (ir < 0 || ir >= gcr[ielectrons].n_rgrid) {
         ERROR("Solar position not within grid, aborting normalization");
         return 1;
      }

      buf.str("");
      buf<<"Grid point for normalization: ir r[ir] iz z[iz] ip Ekin[ip] "<<ir<<" " <<gcr[ielectrons].r[ir]<<" " <<iz <<" "<<gcr[ielectrons].z[iz]<<" "<<ip<<" "<< gcr[ielectrons].Ekin[ip];
      INFO(buf.str());
      
      v1 = gcr[ielectrons].cr_density.d2[ir  ][iz].s[ip];
      v2 = gcr[ielectrons].cr_density.d2[ir+1][iz].s[ip];
      v3 = gcr[ielectrons].cr_density.d2[ir  ][iz].s[ip+1];
      v4 = gcr[ielectrons].cr_density.d2[ir+1][iz].s[ip+1];
      v5 = v1+(r0-gcr[ielectrons].r[ir])/galdef.dr*(v2-v1); // r0 ip
      v6 = v3+(r0-gcr[ielectrons].r[ir])/galdef.dr*(v4-v3); // r0 ip+1
      
    }//n_spatial_dimensions==2
    
    if (3 == galdef.n_spatial_dimensions) {
      
      const int ix = int((r0-galdef.x_min)/galdef.dx + 0.5); //IMOS20060420
      const int iy = int((0.-galdef.y_min)/galdef.dy + 0.5); //IMOS20060420

      if (ix < 0 || ix >= gcr[ielectrons].n_xgrid) {
         ERROR("Solar position not within grid, aborting normalization");
         return 1;
      }
      if (iy < 0 || iy >= gcr[ielectrons].n_ygrid) {
         ERROR("Solar position not within grid, aborting normalization");
         return 1;
      }
      
      buf.str("");
      buf<<"Grid point for normalization: ix x[ix] iy y[iy] iz z[iz] ip Ekin[ip] "<<ix<<" " <<gcr[ielectrons].x[ix]<<" "<<iy<<" "<<gcr[ielectrons].y[iy]<<" " <<iz <<" "<<gcr[ielectrons].z[iz]<<" "<<ip<<" "<< gcr[ielectrons].Ekin[ip];            //AWS20001121
      INFO(buf.str());
      
      v1 = gcr[ielectrons].cr_density.d3[ix  ][iy][iz].s[ip];   //AWS20001121
      v2 = gcr[ielectrons].cr_density.d3[ix+1][iy][iz].s[ip];   //AWS20001121
      v3 = gcr[ielectrons].cr_density.d3[ix  ][iy][iz].s[ip+1]; //AWS20001121
      v4 = gcr[ielectrons].cr_density.d3[ix+1][iy][iz].s[ip+1]; //AWS20001121
      v5 = v1+(r0-gcr[ielectrons].x[ix])/galdef.dx*(v2-v1); // r0 ip
      v6 = v3+(r0-gcr[ielectrons].x[ix])/galdef.dx*(v4-v3); // r0 ip+1

      /*for (auto ix(0); ix < gcr[ielectrons].n_xgrid; ++ix)
	for (auto iy(0); iy < gcr[ielectrons].n_ygrid; ++iy)
	  for (auto iz(0); iz < gcr[ielectrons].n_zgrid; ++iz) {

	    std::cout << ix << " " << iy << " " << iz << " ";
	    for (auto ip(0); ip < gcr[ielectrons].n_pgrid; ++ip)
	      std::cout << gcr[ielectrons].cr_density.d3[ix][iy][iz].s[ip] << " ";
	    std::cout << std::endl;

	  }
	
      exit(0);
      */
    }//n_spatial_dimensions==3
    
    vnorm = exp(log(v5) + log(galdef.electron_norm_Ekin/gcr[ielectrons].Ekin[ip])/log(galdef.Ekin_factor)*log(v6/v5));//AWS20101202: vnorm defined at start
    
    // Abort if vnorm is not positive
    if (v5 <= 0 || v6 <= 0) {
       ERROR("Cannot normalize electrons and positrons, non-positive flux of electrons detected.");
       return 1;
    }

    buf.str("");                                                                                                                                                //AWS20101202
    buf<<"       galdef  electron flux at solar position r0 = "<<r0<<" kpc: " <<galdef.electron_norm_flux<< " at energy " <<galdef.electron_norm_Ekin<<" MeV";  //AWS20101202
    INFO(buf.str());                                                                                                                                            //AWS20101202
    buf.str("");                                                                                                                                                //AWS20101202
    buf<<"normalizing to electron flux at solar position: v1 v2 v3 v4 v5 v6 vnorm  "<<v1<<" " <<v2<<" " <<v3 <<" "<<v4 <<" "<<v5<<" "<< v6<<" "<<vnorm;         //AWS20101202
    INFO(buf.str());
    
  } //if galdef.electron_norm_type == 1      
  
  if (2 == galdef.electron_norm_type || 3 == galdef.electron_norm_type) { //AWS20101202
    
    double CR_luminosity=0;
    Particle particle;                                                                      
    particle=gcr[ielectrons];
    particle.create_transport_arrays(); 
    fill_transport_arrays(particle);
    
    for (int ip = 0; ip < gcr[ielectrons].n_pgrid; ip++) {

      if(galaxy.n_spatial_dimensions==2)
	for(int ir=0;ir<gcr[ielectrons].n_rgrid;ir++)
	  for(int iz=0;iz<gcr[ielectrons].n_zgrid;iz++)
	    {
	      if(galdef.electron_norm_type ==2) // particles s-1
		CR_luminosity+= particle.primary_source_function.d2[ir][iz].s[ip]
		  /particle.beta[ip]*particle.Ekin[ip]
		  *2.0*Pi*particle.r[ir]*galdef.dr*galdef.dz;
	      
	      if(galdef.electron_norm_type ==3) // erg s-1
		CR_luminosity+= particle.primary_source_function.d2[ir][iz].s[ip]
		  /particle.beta[ip]*particle.Ekin[ip]*particle.Ekin[ip]*1.0e6*eV_to_erg
		  *2.0*Pi*particle.r[ir]*galdef.dr*galdef.dz;
	    }
      
      if(galaxy.n_spatial_dimensions==3)
	for(int ix=0;ix<gcr[ielectrons].n_xgrid;ix++)
	  for(int iy=0;iy<gcr[ielectrons].n_ygrid;iy++)
	    for(int iz=0;iz<gcr[ielectrons].n_zgrid;iz++)
	      {
		// weights for non-symmetric case
		double sym_weight=1.0;
		// weights for fully symmetric case
		if(galdef.use_symmetry==1 && iz >0) sym_weight= 8.;
		if(galdef.use_symmetry==1 && iz==0) sym_weight= 4.;// to avoid double-counting at z=0
		
		if(galdef.electron_norm_type ==2) // particles s-1
		  CR_luminosity+= particle.primary_source_function.d3[ix][iy][iz].s[ip] 
		    /particle.beta[ip]*particle.Ekin[ip]
		    *galdef.dx*galdef.dy*galdef.dz*sym_weight;
		
		if(galdef.electron_norm_type ==3) // erg s-1
		  CR_luminosity+= particle.primary_source_function.d3[ix][iy][iz].s[ip] 
		    /particle.beta[ip]*particle.Ekin[ip]*particle.Ekin[ip]*1.0e6*eV_to_erg
		    *galdef.dx*galdef.dy*galdef.dz*sym_weight;
		
	      }		
    }//ip
    
    CR_luminosity *= 4.0*Pi/c*pow(kpc2cm,3)*log(gcr[ielectrons].Ekin_factor);
    vnorm =  CR_luminosity;
    
    buf.str("");                                                    
    if(galdef.electron_norm_type==2) buf<<"normalizing to electron luminosity in particles s-1:  vnorm =" <<vnorm;      
    if(galdef.electron_norm_type==3) buf<<"normalizing to electron luminosity in       erg s-1:  vnorm =" <<vnorm;    
    INFO(buf.str());
    
  } //if galdef.electron_norm_type == 2 or 3
  
  buf.str("");                                                    
  buf<<"galdef.electron_norm_type = "<<galdef.electron_norm_type<<"  vnorm =" <<vnorm;         
  INFO(buf.str());

  // Abort if vnorm is not positive
  if (vnorm <= 0) {
     ERROR("Cannot normalize electrons and positrons, non-positive flux of electrons detected.");
     return 1;
  }
  
  // normalize primary electrons and positrons
  if(galdef.electron_norm_type!=0) { //AWS20101202
    for (int i = 0; i < n_species; i++) {
      
      if ("primary_electrons" == gcr[i].name || "primary_positrons" == gcr[i].name) {
	//galdef.electron_source_normalization *= galdef.electron_norm_flux/vnorm; // IMOS20031016
	gcr[i].cr_density *= galdef.electron_norm_flux/vnorm;
	gcr[i].normalization_factor = galdef.electron_norm_flux/vnorm; //AWS20010121
      }
    }
  }//if
  
  if(galdef.verbose>=2){
    INFO("primary electrons after normalization:");
    gcr[ielectrons].cr_density.print();
  }  
  
  INFO("Exit");
  return 0;
}
