
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * create_galaxy.cc *                            galprop package * 5/22/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"

#include <ErrorLogger.h>
#include <sstream>

#include"B_field_3D_model.h" //AWS20080326

int Galprop::create_galaxy() {//AWS20050816
		INFO("STARTED");
		std::cout<<"Pos1"<<std::endl;

  INFO("Entry");
  INFO("Pos2");

  //cout<<" >>>> create_galaxy"<<endl;
  
  int stat=0;
  double dzz=0.02; // average with this resolution in kpc
  INFO("Pos3");
  
  if (2 == galdef.n_spatial_dimensions){
  	INFO("Pos4");
    galaxy.init(galdef.r_min, galdef.r_max, galdef.dr,  
		galdef.z_min, galdef.z_max, galdef.dz);
  	INFO("Pos5");
  }
  
  if (3 == galdef.n_spatial_dimensions) {
    
    if(galdef.use_symmetry==1) 
      galdef.x_min=galdef.y_min=galdef.z_min=0.; // IMOS20020419
    
    galaxy.init(galdef.x_min, galdef.x_max, galdef.dx,
		galdef.y_min, galdef.y_max, galdef.dy,
		galdef.z_min, galdef.z_max, galdef.dz);
  
  }                  
  
  if ("p" == galdef.p_Ekin_grid) {
    
    galaxy.p_min = galdef.p_min;
    galaxy.p_max = galdef.p_max;
    galaxy.p_factor = galdef.p_factor;

    galaxy.n_pgrid = int(log(galaxy.p_max/galaxy.p_min)/log(galaxy.p_factor) + 1.9);
 
    galaxy.p.resize(galaxy.n_pgrid);

    for (int i = 0; i < galaxy.n_pgrid; ++i)
      galaxy.p[i] = exp(log(galaxy.p_min) + i*log(galaxy.p_factor));

  }

  if ("Ekin" == galdef.p_Ekin_grid) {

    galaxy.p_min = galdef.Ekin_min;
    galaxy.p_max = galdef.Ekin_max;
    galaxy.p_factor = galdef.Ekin_factor;

    galaxy.n_pgrid = int(log(galaxy.p_max/galaxy.p_min)/log(galaxy.p_factor) + 1.9);
    
    galaxy.p.resize(galaxy.n_pgrid);
    
    for (int i = 0; i < galaxy.n_pgrid; ++i)
      galaxy.p[i] = exp(log(galaxy.p_min) + i*log(galaxy.p_factor));

  }

  galaxy.E_gamma_min = galdef.E_gamma_min;
  galaxy.E_gamma_max = galdef.E_gamma_max;
  galaxy.E_gamma_factor = galdef.E_gamma_factor;
  
  galaxy.n_E_gammagrid = int(log(galaxy.E_gamma_max/galaxy.E_gamma_min)/log(galaxy.E_gamma_factor) + 1.001);
  
  galaxy.E_gamma.resize(galaxy.n_E_gammagrid);
  
  for (int iEgamma = 0; iEgamma < galaxy.n_E_gammagrid; ++iEgamma)
    galaxy.E_gamma[iEgamma] = 
      exp(log(galaxy.E_gamma_min) + iEgamma*log(galaxy.E_gamma_factor));

  //Moved from create_galaxy  Gulli20071003
  
  galaxy.nu_synch_min = galdef.nu_synch_min;
  galaxy.nu_synch_max = galdef.nu_synch_max;
  galaxy.nu_synch_factor = galdef.nu_synch_factor;
  galaxy.n_nu_synchgrid = int(log(galaxy.nu_synch_max/galaxy.nu_synch_min)/log(galaxy.nu_synch_factor) + 1.001);
  
  galaxy.nu_synch.resize(galaxy.n_nu_synchgrid);
  
  for (int iS = 0; iS < galaxy.n_nu_synchgrid; ++iS)
    galaxy.nu_synch[iS] = exp(log(galaxy.nu_synch_min) + iS*log(galaxy.nu_synch_factor));
  
  // GAS and B-FIELD DISTRIBUTION
  
  nH_set_model(galdef); //AWS20090814


  galaxy.n_HI =1.e-6; //  non-zero values to avoid problem in energy loss logarithm
  galaxy.n_H2 =1.e-6;
  galaxy.n_HII=1.e-6;
  
  if (2 == galdef.n_spatial_dimensions) { // use average at y=0

    galaxy.fPairAbsorptionInvLength.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);
    galaxy.fPairAbsorptionInvLength = 0;
      
    for (int ir = 0; ir < galaxy.n_rgrid; ++ir) {
      
      for (int iz = 0; iz < galaxy.n_zgrid; ++iz) {

	galaxy.n_HI.d2[ir][iz].s[0] = 
	  nH_av(galaxy.r[ir], 0, galaxy.z[iz], galaxy.dz, dzz, nHI3D); //IMOS20080114
	//		nHI_av (galaxy.r[ir],0.0,galaxy.z[iz],galaxy.dz,dzz);
	      
	galaxy.n_H2.d2[ir][iz].s[0]=
	  nH_av(galaxy.r[ir], 0, galaxy.z[iz], galaxy.dz, dzz, fX_CO(galaxy.r[ir]), nH23D); //IMOS20080114
	//		nH2_av (galaxy.r[ir],0.0,galaxy.z[iz],galaxy.dz,dzz);
	      
	galaxy.n_HII.d2[ir][iz].s[0]=
	  nH_av(galaxy.r[ir], 0, galaxy.z[iz], galaxy.dz, dzz, nHII3D); //IMOS20080114
	//		nHII_av(galaxy.r[ir],0.0,galaxy.z[iz],galaxy.dz,dzz); 
	      
	galaxy.B_field.d2[ir][iz].s[0]=
	  B_field_model(galaxy.r[ir], 0, galaxy.z[iz], galdef.B_field_model);

	// 3D model total field  
	if (galdef.synchrotron>=2)  //AWS20101109
	{  
	  int debug=0; if(galdef.verbose==-1100)debug=1; //AWS20101109
	  galaxy.B_field.d2[ir][iz].s[0] = B_field_3D_model_tot(galdef.B_field_name, galdef.B_field_parameters, galaxy.r[ir], galaxy.z[iz],debug)/1.0e4;     // Gauss->Tesla   //AWS20101109
	  galaxy.B_field_regular.d2[ir][iz].s[0] = B_field_3D_model_abs(galdef.B_field_name, galdef.B_field_parameters, galaxy.r[ir], galaxy.z[iz],debug,B_REGULAR)/1.0e4;     // Gauss->Tesla   //AWS20101109
	  galaxy.B_field_random.d2[ir][iz].s[0] = B_field_3D_model_abs(galdef.B_field_name, galdef.B_field_parameters, galaxy.r[ir], galaxy.z[iz],debug,B_RANDOM)/1.0e4;     // Gauss->Tesla   //AWS20101109

	}

      }

    }
    
  }
  
  if (3 == galdef.n_spatial_dimensions) {

    galaxy.fPairAbsorptionInvLength.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);
    galaxy.fPairAbsorptionInvLength = 0;
 
#pragma omp parallel for schedule(static) default(shared)
    for (int ix = 0; ix < galaxy.n_xgrid; ++ix) {

      for (int iy = 0; iy < galaxy.n_ygrid; ++iy) {

	const double r = sqrt(galaxy.x[ix]*galaxy.x[ix] + galaxy.y[iy]*galaxy.y[iy]);

	for (int iz = 0; iz < galaxy.n_zgrid; ++iz) {

	  
	  galaxy.n_HI.d3[ix][iy][iz].s[0] =
	    nH_av(galaxy.x[ix], galaxy.y[iy], galaxy.z[iz], galaxy.dz, dzz, nHI3D); //IMOS20080114
	  //		    nHI_av (galaxy.x[ix],galaxy.y[iy],galaxy.z[iz],galaxy.dz,dzz);
		  
	  galaxy.n_H2.d3[ix][iy][iz].s[0] =
	    nH_av(galaxy.x[ix], galaxy.y[iy], galaxy.z[iz], galaxy.dz, dzz, fX_CO(r), nH23D); //IMOS20080114
	  //		    nH2_av (galaxy.x[ix],galaxy.y[iy],galaxy.z[iz],galaxy.dz,dzz);
		  
	  galaxy.n_HII.d3[ix][iy][iz].s[0] =
	    nH_av(galaxy.x[ix], galaxy.y[iy], galaxy.z[iz], galaxy.dz, dzz, nHII3D);  //IMOS20080114
	  //		    nHII_av(galaxy.x[ix],galaxy.y[iy],galaxy.z[iz],galaxy.dz,dzz);  
		  
	  galaxy.B_field.d3[ix][iy][iz].s[0] =
	    B_field_model(galaxy.x[ix], galaxy.y[iy], galaxy.z[iz], galdef.B_field_model);

	  // 3D model total field       
	  if (galdef.synchrotron>=2)                       //AWS20101109
	  {
	    int debug=0; if(galdef.verbose==-1100)debug=1; //AWS20101109
	    galaxy.B_field.d3[ix][iy][iz].s[0] = B_field_3D_model_tot(galdef.B_field_name, galdef.B_field_parameters, galaxy.x[ix], galaxy.y[iy], galaxy.z[iz],debug)/1.0e4;  //AWS20101109  
	    galaxy.B_field_regular.d3[ix][iy][iz].s[0] = B_field_3D_model_abs(galdef.B_field_name, galdef.B_field_parameters, galaxy.x[ix], galaxy.y[iy],galaxy.z[iz],debug,B_REGULAR)/1.0e4;  //AWS20101109  
	    galaxy.B_field_random.d3[ix][iy][iz].s[0] = B_field_3D_model_abs(galdef.B_field_name, galdef.B_field_parameters, galaxy.x[ix], galaxy.y[iy], galaxy.z[iz],debug,B_RANDOM)/1.0e4;  //AWS20101109  

	  }
	   
	}
	
      }
    
    }
  
  }

  if (galdef.verbose>=1) {

      INFO("galaxy.n_HI:   ");galaxy.n_HI.print();
      INFO("galaxy.n_H2:   ");galaxy.n_H2.print();
      INFO("galaxy.n_HII:  ");galaxy.n_HII.print();
      INFO("galaxy.B_field:");galaxy.B_field.print();
      INFO("galaxy.B_field_regular:");galaxy.B_field_regular.print();
      INFO("galaxy.B_field_random:");galaxy.B_field_random.print();

  }
  
  // SKYMAP PARAMETERS
  // moved to before gas surveys since required there                                    AWS20050913 
  
  galaxy.d_long = galdef.d_long;
  galaxy.long_min = galdef.long_min;
  galaxy.long_max = galdef.long_max;
  
  galaxy.d_lat = galdef.d_lat;
  galaxy.lat_min = galdef.lat_min;
  galaxy.lat_max = galdef.lat_max;
  
  galaxy.n_long = int((galaxy.long_max-galaxy.long_min)/galaxy.d_long + 1.001);
  galaxy.n_lat = int((galaxy. lat_max-galaxy. lat_min)/galaxy.d_lat  + 1.001);
  
  ostringstream buf;
  buf<<"    gamma-ray, synchrotron skymap arrays: n_long,nlat="<<galaxy.n_long<<" "<<galaxy.n_lat;
  DEBUGLOG(buf.str());
  
  
  // READING THE INTERSTELLAR RADIATION FIELD
  
  if (galdef.ISRF_factors.size() >= 3) {

    galaxy.fISRFFactors.resize(3);
   
    for (unsigned long i = 0; i < 3; ++i)
      galaxy.fISRFFactors[i] = galdef.ISRF_factors[i];

  }

  stat = read_isrf(galdef.ISRF_filetype);
  
  if (gen_isrf_energy_density() != 0) 
    return 1;
  
  if(galdef.verbose >= 1){
     INFO("      galaxy.print");
     galaxy.print();
  }

  INFO("Exit");

  //cout<<" <<<< create_galaxy"<<endl;
  return stat;

}
