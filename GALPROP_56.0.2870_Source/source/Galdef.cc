//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Galdef.cc *                                   galprop package * 10/12/2003
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <ErrorLogger.h>

#include "Galdef.h"
#include <Parameters.h>

#include <constants.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace utl;

int Bufferlength=1000;//IMOS20080114

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Galdef::Galdef() :
  ISRF_healpixOrder ( 0 )
  {
  
}

Galdef::~Galdef() {

}

int Galdef::read ( const string& theVersion,
                   const string& runNumber,
                   const string& galdefDirectory ) {

    INFO ( "Entry" );

    gdDirectory = galdefDirectory;

    //strcpy(version,version_);
    //strcpy(run_no,run_no_);

    version = theVersion;
    run_no = runNumber;

    AssignDefaultParameters();

    ostringstream bufID;
    bufID << version << "_" << run_no;
    //string tmpStringBuf = bufID.str();
    //tmpStringBuf.resize(99);
    //strcpy(galdef_ID, tmpStringBuf.c_str());
    galdef_ID = bufID.str();

    ostringstream galdef_file;
    galdef_file << galdefDirectory << "galdef_" << galdef_ID;

    //char galdef_file[100];
    //strcpy(galdef_file,galdef_directory);
    //strcat(galdef_file,"galdef_");
    //strcat(galdef_file,galdef_ID);

    ostringstream oBuf;
    oBuf << "Reading from " << galdef_file.str();
    INFO ( oBuf.str() );

    //char  filename[100];
    //strcpy(filename,galdef_file);

    const string filename = galdef_file.str();

    //Use Parameters class to parse the galdef file
    ifstream is(filename.c_str());
    if (! is.good()) {
      ostringstream buf;
      buf << "Galdef file " << filename << " not found";
      INFO ( buf.str() );
      return -1;
    }

    Parameters pars(is);
    is.close();

    return AssignParameters( pars );

}

int Galdef::AssignParameters( const Parameters &PARAMETERS_NAME) {

    //Start by reading and setting the verbosity
    //-2 means Errors and Fatal only
    //-1 adds Warnings
    //0 adds info
    //everything else adds debug (requires compilation with debug enabled).
    //Special values as before 
    SETPAROPT(verbose);
    if ( verbose == -2 ) {
      INFO("Setting logger to errors only");
      utl::ErrorLogger::GetInstance().SetSeverity(utl::ErrorLogger::eError);
    } else if ( verbose == -1 ){
      INFO("Setting logger to errors and warnings only");
      utl::ErrorLogger::GetInstance().SetSeverity(utl::ErrorLogger::eWarning);
    } else if ( verbose == 0 ) {
      INFO("Setting logger to normal (errors, warnings, and info)");
      utl::ErrorLogger::GetInstance().SetSeverity(utl::ErrorLogger::eInfo);
    } else {
      INFO("Setting logger to debug.  Requires compilation with --enable-debug.");
      utl::ErrorLogger::GetInstance().SetSeverity(utl::ErrorLogger::eDebug);
      utl::ErrorLogger::GetInstance().SetVerbosity(utl::ErrorLogger::eVerbose);
    }

    SETPAROPT(n_spatial_dimensions);

    //radial grid
    SETPAROPT(r_min);
    SETPAROPT(r_max);
    SETPAROPT(dr);

    //z-grid
    SETPAROPT(z_min);
    SETPAROPT(z_max);
    SETPAROPT(dz);

    //x-grid
    SETPAROPT(x_min);
    SETPAROPT(x_max);
    SETPAROPT(dx);

    //y-grid
    SETPAROPT(y_min);
    SETPAROPT(y_max);
    SETPAROPT(dy);

    //momentum grid
    SETPAROPT(p_min);
    SETPAROPT(p_max);
    SETPAROPT(p_factor);

    //kinetic energy grid
    SETPAROPT(Ekin_min);
    SETPAROPT(Ekin_max);
    SETPAROPT(Ekin_factor);

    //p||Ekin option
    SETPAROPT(p_Ekin_grid);

    //gamma-ray energy grid
    SETPAROPT(E_gamma_min);
    SETPAROPT(E_gamma_max);
    SETPAROPT(E_gamma_factor);

    //mode of integration over particle spectrum
    SETPAROPT(integration_mode);

    //synchrotron grid
    SETPAROPT(nu_synch_min);
    SETPAROPT(nu_synch_max);
    SETPAROPT(nu_synch_factor);

    //longitude-latitude grid
    SETPAROPT(long_min);
    SETPAROPT(long_max);
    SETPAROPT(lat_min);
    SETPAROPT(lat_max);

    SETPAROPT(d_long);
    SETPAROPT(d_lat);

    //Healpix order
    SETPAROPT(healpix_order);

    //line of sight integration IMOS20080114, only works for old integration mode
    SETPAROPT(lat_substep_number);
    SETPAROPT(LoS_step);
    SETPAROPT(LoS_substep_number);

    //Control new los integration code
    SETPAROPT(los_integration_mode);
    SETPAROPT(los_integration_accuracy);
    SETPAROPT(LoS_minStep);

    // HEALPix order for anisotropic IC calculation
    SETPAROPT(anisoHealpixOrder);

    //space diffusion
    SETPAROPT(D0_xx);

    // diffusion coefficient parametrization
    SETPAROPT(D_rigid_br);
    SETPAROPT(D_g_1);
    SETPAROPT(D_g_2);
    SETPAROPT(D_eta);

    // Reference rigidity defaults to break rigidity for backwards compatibility
    D_rigid_ref = D_rigid_br;
    SETPAROPT(D_rigid_ref);

    SETPAROPT(Dxx_plane_scale);
    SETPAROPT(Dxx_plane_scale_height);

    SETPAROPT(B_dep_diffusion);
    // Set the minimum and maximum value for the diffusion coeff at 10 GV, used when D_xx is linked with B-field
    SETPAROPT(D_xx_max);
    SETPAROPT(D_xx_min);


    //diffusive reacceleration: Alfven speed
    SETPAROPT(v_Alfven);
    SETPAROPT(diff_reacc);
    //damping
    SETPAROPT(damping_p0);
    SETPAROPT(damping_max_path_L);
    SETPAROPT(damping_const_G);
    SETPAROPT(damping_const_K);

    // convection
    SETPAROPT(convection);
    SETPAROPT(v0_conv);
    SETPAROPT(dvdz_conv);

    //injection spectra limitations
    SETPAROPT(rigid_min);
    SETPAROPT(rigid_max);
    SETPAROPT(inj_Ekin_min);
    SETPAROPT(inj_Ekin_max);

    //Source classes
    source_classes.clear();
    SETPAROPT(source_class_files);
    if (source_class_files.size() == 0) {
       utl::Parameters newPars(PARAMETERS_NAME);
       newPars.setParameter("SourceClassType", "Compatibility");
       source_classes.push_back(SourceClass::create(std::move(newPars)));
    } else {
       //Try to read full path and then in GALDEF directory
       for (size_t i(0); i < source_class_files.size(); ++i) {

          INFO("Reading source class information");

          std::string file = source_class_files[i];
          INFO("Trying "+file);
          ifstream is(file.c_str());
          if (! is.good()) {
             //May not be needed
             if (is.is_open()) 
                is.close();
             file = gdDirectory + "/" + source_class_files[i];
             INFO("Trying "+file);
             is.open(file.c_str());

             if (! is.good()) {
                std::ostringstream buf;
                buf<<"Failed to open parameter file \""<<source_class_files[i]<<"\" for source class";
                ERROR(buf.str());
                return -1;
             }
          }

          std::ostringstream buf;
          buf<<"Reading source class information from \""<<file<<"\"";
          INFO(buf.str());

          source_classes.push_back(SourceClass::create(utl::Parameters(is)));
          is.close();
       }
    }

    //other parameters
    SETPAROPT(He_H_ratio);


    //X_CO
    SETPAROPT(n_X_CO);
    SETPAROPT(X_CO);
    SETPAROPT(X_CO_parameters);
    
    //Need 4 values for X_CO_parameters if n_X_CO = 2
    if (n_X_CO == 2){
      if (X_CO_parameters.size() < 4) {
	FATAL("X_CO_parameters should contain 4 values if n_X_CO is 2");
	return -1;
      }
    }
    
    //X_CO parameters for n_X_CO 1, using tabulated values
    SETPAROPT(X_CO_values);
    SETPAROPT(X_CO_radius);
    if (X_CO_values.size() != X_CO_radius.size()) {
      FATAL("X_CO_values and X_CO_radius should have the same size");
      return -1;
    }

    if (n_X_CO == 1 || n_X_CO == 3) {
       if (X_CO_values.size() == 0) {
          FATAL("Need at least 1 value in X_CO_values and X_CO_radius for n_X_CO = 1 and n_X_CO = 3");
          return -1;
       }
    }

    SETPAROPT(propagation_X_CO);

    SETPAROPT(nHI_model);
    SETPAROPT(nH2_model);
    SETPAROPT(nHII_model);

    SETPAROPT(HI_xmlFilename);
    SETPAROPT(H2_xmlFilename);
    SETPAROPT(HII_xmlFilename);

    SETPAROPT(HII_Te);              //AWS20110701
    SETPAROPT(HII_clumping_factor); //AWS20110701

    //CO gas file
    SETPAROPT(COR_filename);

    //HI gas file
    SETPAROPT(HIR_filename);

    SETPAROPT(fragmentation);
    SETPAROPT(momentum_losses);
    SETPAROPT(radioactive_decay);
    SETPAROPT(K_capture);
    SETPAROPT(ionization_rate);

    SETPAROPT(ionization_losses);
    SETPAROPT(coulomb_losses);
    SETPAROPT(bremss_losses);
    SETPAROPT(IC_losses);
    SETPAROPT(sync_losses);

    //time grid
    SETPAROPT(start_timestep);
    SETPAROPT(end_timestep);
    SETPAROPT(timestep_factor);
    SETPAROPT(timestep_repeat);
    SETPAROPT(timestep_repeat2);

    SETPAROPT(timestep_print);
    SETPAROPT(timestep_diagnostics);
    SETPAROPT(control_diagnostics);

    SETPAROPT(solution_convergence); //AWS20110118
    SETPAROPT(solution_method);      //AWS20110118
    SETPAROPT(solution_rel_accuracy);

    //Need to assert that we have correct grid size with the vectorized solvers
    if (solution_method == 3 || solution_method == 4) {
       size_t nE = log(Ekin_max/Ekin_min)/log(Ekin_factor) + 1;
       if (nE % 8 > 3) {
          nE = (nE/8+1)*8;
       } else {
          nE = (nE/8)*8;
       }
       Ekin_factor = exp(log(Ekin_max/Ekin_min)/(nE-1))*0.999999;  //Have it the tiniest bit smaller to avoid rounding issues

       std::ostringstream oss;
       oss<<"Adjusting Ekin_factor to "<<Ekin_factor<<" for alignment purposes";
       WARNING(oss.str());
       oss.str("");

       if (n_spatial_dimensions == 2) // === 2D ===
       {
          size_t nr = (r_max - r_min)/dr + 1;
          if (nr % 8 > 2) { //Favor increased resolution rather than decreased
             nr = (nr/8 + 1)*8;
          } else {
             nr = (nr/8)*8;
          }
          dr = (r_max - r_min)/(nr-1)*0.999999;//Have it the tiniest bit smaller to avoid rounding issues

          oss<<"Adjusting dr to "<<dr<<" for alignment purposes";
          WARNING(oss.str());
          oss.str("");
       }
       else // === 3D ===
       {
          size_t nx = (x_max - x_min)/dx + 1;
          if (nx % 8 > 2) { //Favor increased resolution rather than decreased
             nx = (nx/8 + 1)*8;
          } else {
             nx = (nx/8)*8;
          }
          dx = (x_max - x_min)/(nx-1)*0.999999;//Have it the tiniest bit smaller to avoid rounding issues

          oss<<"Adjusting dx to "<<dx<<" for alignment purposes";
          WARNING(oss.str());
          oss.str("");
       }
    }


    //other parameters
    SETPAROPT(network_iterations);
    SETPAROPT(network_iter_compl);

    SETPAROPT(prop_r);
    SETPAROPT(prop_x);
    SETPAROPT(prop_y);
    SETPAROPT(prop_z);

    SETPAROPT(prop_p);

    SETPAROPT(use_symmetry);

    SETPAROPT(vectorized);


    //nuclei to include & abundances
    SETPAROPT(max_Z);

    use_Z.resize(max_Z+1, 0);
    for ( int iZ=1; iZ<=max_Z; iZ++ ) {
      ostringstream parname;
      parname << "use_Z_" << iZ;
      SETPARSTROPT(use_Z[iZ], parname.str());
    }
    
    SETPAROPT(total_cross_section);
    
    SETPAROPT(cross_section_option);

    SETPAROPT(t_half_limit);

    
    //CR species to include
    SETPAROPT(primary_electrons);
    SETPAROPT(primary_positrons);
    SETPAROPT(secondary_positrons);
    SETPAROPT(secondary_electrons);
    SETPAROPT(knock_on_electrons);
    SETPARNAMEOPT(secondary_antiprotons,secondary_antiproton);
    SETPARNAMEOPT(tertiary_antiprotons,tertiary_antiproton);
    SETPAROPT(secondary_protons);

    SETPARNAMEOPT(pair_production,pairproduction);

    //Gamma-ray calculations
    SETPAROPT(gamma_rays);
    SETPAROPT(pi0_decay);
    SETPAROPT(IC_isotropic);
    SETPAROPT(IC_anisotropic);
    SETPAROPT(bremss);
    SETPAROPT(synchrotron);
    SETPAROPT(free_free_absorption); //AWS20110701

    // Global luminosity calculation
    SETPAROPT(globalLuminosities);

    // Local bubble
    SETPAROPT(local_bubble_radius);
    SETPAROPT(local_bubble_gas_fraction);
    SETPAROPT(local_bubble_source_fraction);

    //HI, CO  surveys                                                               //AWS20051309
    //HI_survey=8;                                                     //default
    //strcpy(parstring,                              "HI_survey"           );
    //stat= read_galdef_parameter(filename,parstring,&HI_survey            );
    //CO_survey=8;                                                     //default
    //strcpy(parstring,                              "CO_survey"           );
    //stat= read_galdef_parameter(filename,parstring,&CO_survey            );

    //magnetic field
    SETPAROPT(B_field_model);

    SETPAROPT(B_field_name);

    SETPAROPT(B_field_parameters);

    //    if ( B_field_parameters.size() !=10 ) {FATAL("need exactly 10 B-field parameters in this version!"); exit ( 1 );}             //AWS20110520: no such restriction needed

    //strcpy(workstring,"-1,-2,-3,-4,-5,-6,-7,-8,-9,-10");                          //    default if absent in galdef file         //AWS20080312
    //stat= read_galdef_parameter(filename,parstring, workstring           );                                                      //AWS20080312
    //// this is not flexible, need to find a better way to read in arrays                                                         //AWS20080313
    //sscanf(workstring,"%le,%le,%le,%le,%le,%le,%le,%le,%le,%le",                                                                 //AWS20080312
//	 &B_field_parameters[0],&B_field_parameters[1],&B_field_parameters[2],&B_field_parameters[3],&B_field_parameters[4],   //AWS20080312
//	 &B_field_parameters[5],&B_field_parameters[6],&B_field_parameters[7],&B_field_parameters[8],&B_field_parameters[9])  ;//AWS20080312


    //ISRF file
    SETPAROPT(ISRF_file);

    //ISRF file type
    SETPAROPT(ISRF_filetype);

    // ISRF healpix order for skymaps
    SETPAROPT(ISRF_healpixOrder);

    //ISRF factors for IC calculation
    SETPAROPT(ISRF_factors);
    if (ISRF_factors.size() < 3) {
      ISRF_factors.resize(3,1.0);
    }

    //spectra normalization
    SETPAROPT(proton_norm_Ekin);
    SETPAROPT(proton_norm_flux);
    SETPAROPT(proton_norm_type);  //AWS20101201

    SETPAROPT(electron_norm_Ekin);
    SETPAROPT(electron_norm_flux);
    SETPAROPT(electron_norm_type);//AWS20101201

    //dark matter (DM) parameters  IMOS20050912
    SETPAROPT(DM_positrons);
    SETPAROPT(DM_electrons);
    SETPAROPT(DM_antiprotons);
    SETPAROPT(DM_gammas);

    SETPAROPT(DM_double0);
    SETPAROPT(DM_double1);
    SETPAROPT(DM_double2);
    SETPAROPT(DM_double3);
    SETPAROPT(DM_double4);
    SETPAROPT(DM_double5);
    SETPAROPT(DM_double6);
    SETPAROPT(DM_double7);
    SETPAROPT(DM_double8);
    SETPAROPT(DM_double9);

    SETPAROPT(DM_int0);
    SETPAROPT(DM_int1);
    SETPAROPT(DM_int2);
    SETPAROPT(DM_int3);
    SETPAROPT(DM_int4);
    SETPAROPT(DM_int5);
    SETPAROPT(DM_int6);
    SETPAROPT(DM_int7);
    SETPAROPT(DM_int8);
    SETPAROPT(DM_int9);

    //output controls

    vector<double> cameraLocation;

    SETPAROPT(cameraLocation);

    fCameraLocation.clear();
    fCameraLocation.insert(fCameraLocation.begin(), cameraLocation.begin(), cameraLocation.end());

    if (fCameraLocation.size() < 3) {
      
      fCameraLocation.resize(3);
      
      fCameraLocation[0] = Rsun;
      fCameraLocation[1] = 0;
      fCameraLocation[2] = 0;

    } else if (fCameraLocation.size() > 3) {

      INFO("Camera location specified with more than 3 entries.");

    }

    SETPAROPT(skymap_format);
    SETPAROPT(output_gcr_full);

    SETPAROPT(warm_start);

    SETPAROPT(test_suite);

    print();

    INFO ( "Exit" );
    return 0;

}

void Galdef::print() {

    ostringstream os;

    os<<"  ======= galdef: "<<galdef_ID  <<endl;
    os<<"  version  ="<<version <<endl;
    os<<"  run_no   ="<<run_no      <<endl;

    os<<"  n_spatial_dimensions ="<<n_spatial_dimensions   <<endl;

    os<<"  r_min    ="<<r_min   <<endl;
    os<<"  r_max    ="<<r_max   <<endl;
    os<<"  dr       ="<<dr      <<endl;
    os<<"  x_min    ="<<x_min   <<endl;
    os<<"  x_max    ="<<x_max   <<endl;
    os<<"  dx       ="<<dx      <<endl;
    os<<"  y_min    ="<<y_min   <<endl;
    os<<"  y_max    ="<<y_max   <<endl;
    os<<"  dy       ="<<dy      <<endl;
    os<<"  z_min    ="<<z_min   <<endl;
    os<<"  z_max    ="<<z_max   <<endl;
    os<<"  dz       ="<<dz      <<endl;
    os<<"  p_min    ="<<p_min   <<endl;
    os<<"  p_max    ="<<p_max   <<endl;
    os<<"  p_factor ="<<p_factor<<endl;
    os<<"  Ekin_min ="<<   Ekin_min      <<endl;
    os<<"  Ekin_max ="<<   Ekin_max      <<endl;
    os<<"  Ekin_factor ="<<Ekin_factor   <<endl;

    os<<"  p_Ekin_grid ="<<p_Ekin_grid   <<endl;

    os<<"  E_gamma_min   ="<<   E_gamma_min      <<endl;
    os<<"  E_gamma_max   ="<<   E_gamma_max      <<endl;
    os<<"  E_gamma_factor ="<<   E_gamma_factor   <<endl;
    os<<"  integration_mode ="<< integration_mode <<endl;

    os<<"  nu_synch_min    ="<<  nu_synch_min      <<endl;
    os<<"  nu_synch_max    ="<<  nu_synch_max      <<endl;
    os<<"  nu_synch_factor ="<<  nu_synch_factor   <<endl;

    os<<"  long_min      ="<<   long_min         <<endl;
    os<<"  long_max      ="<<   long_max         <<endl;
    os<<"  lat_min       ="<<    lat_min         <<endl;
    os<<"  lat_max       ="<<    lat_max         <<endl;
    os<<"  d_long        ="<<   d_long           <<endl;
    os<<"  d_lat         ="<<   d_lat            <<endl;
    os<<"  healpix_order ="<<   healpix_order    <<endl;
    os<<"  lat_substep_number ="<< lat_substep_number <<endl;
    os<<"  LoS_step      ="<<      LoS_step           <<endl;
    os<<"  LoS_minStep   ="<<      LoS_minStep        <<endl;
    os<<"  LoS_substep_number ="<< LoS_substep_number <<endl;
    os<<"  los_integration_mode ="<< los_integration_mode <<endl;
    os<<"  los_integration_accuracy ="<< los_integration_accuracy <<endl;
    os<<"  anisoHealpixOrder ="<< anisoHealpixOrder <<endl;

    os<<"  D0_xx       ="<<D0_xx         <<endl;
    os<<"  D_rigid_ref ="<<D_rigid_ref   <<endl;
    os<<"  D_rigid_br  ="<<D_rigid_br    <<endl;
    os<<"  D_g_1       ="<<D_g_1         <<endl;
    os<<"  D_g_2       ="<<D_g_2         <<endl;
    os<<"  D_eta       ="<<D_eta         <<endl;
    os<<"  Dxx_plane_scale="<<Dxx_plane_scale<<endl;
    os<<"  Dxx_plane_scale_height="<<Dxx_plane_scale_height<<endl;
    os<<"  B_dep_diffusion  ="<<B_dep_diffusion <<endl;
    os<<"  D_xx_max    ="<<D_xx_max      <<endl;
    os<<"  D_xx_min    ="<<D_xx_min      <<endl;
    os<<"  diff_reacc  ="<<diff_reacc    <<endl;
    os<<"  v_Alfven    ="<<v_Alfven      <<endl;
    os<<"  convection  ="<<convection    <<endl; //AWS20010323
    os<<"  v0_conv     ="<<v0_conv       <<endl; //AWS20010323
    os<<"  dvdz_conv   ="<<dvdz_conv     <<endl; //AWS20010323

    os<<"  damping_p0         ="<<damping_p0        <<endl;//IMOS20060330
    os<<"  damping_const_G    ="<<damping_const_G   <<endl;//IMOS20060330
    os<<"  damping_const_K    ="<<damping_const_K   <<endl;//IMOS20060330
    os<<"  damping_max_path_L ="<<damping_max_path_L<<endl;//IMOS20060330

    os<<"  rigid_min     ="<<rigid_min       <<endl;
    os<<"  rigid_max     ="<<rigid_max       <<endl;
    os<<"  inj_Ekin_min     ="<<inj_Ekin_min       <<endl;
    os<<"  inj_Ekin_max     ="<<inj_Ekin_max       <<endl;

    os<<"  He_H_ratio        ="<<He_H_ratio          <<endl;
    os<<"  n_X_CO            ="<<n_X_CO              <<endl;//IMOS20080114
    os<<"  X_CO              ="<<X_CO                <<endl;//IMOS20080114
    if ( X_CO_parameters.size() > 0) {
	    os<<"  X_CO_parameters   =";
	    for (int i = 0; i < X_CO_parameters.size(); ++i) {
		    os<<X_CO_parameters[i]<<" ";
	    }
	    os<<endl;
    }
    if ( X_CO_values.size() > 0 ) {
        os<<"   X_CO_values      =";
        for ( int i = 0; i < X_CO_values.size(); ++i )
            os<<X_CO_values[i]<<", ";
        os<<endl;
        os<<"   X_CO_radius      =";
        for ( int i = 0; i < X_CO_radius.size(); ++i )
            os<<X_CO_radius[i]<<", ";
        os<<endl;
    }
    os<<"  propagation_X_CO  ="<<propagation_X_CO    <<endl;//AWS20090623
    os<<"  nHI_model         ="<<nHI_model           <<endl;//AWS20090814
    os<<"  nH2_model         ="<<nH2_model           <<endl;//AWS20090814
    os<<"  nHII_model        ="<<nHII_model          <<endl;//AWS20090814
    os<<"  HI_xmlFilename    ="<<HI_xmlFilename      <<endl;
    os<<"  H2_xmlFilename    ="<<H2_xmlFilename      <<endl;
    os<<"  HII_xmlFilename   ="<<HII_xmlFilename     <<endl;
    os<<"  HII_Te             ="<<HII_Te             <<endl;//AWS20110701
    os<<"  HII_clumping_factor="<<HII_clumping_factor<<endl;//AWS20110701

    os<<"  COR_filename      ="<<COR_filename        <<endl;//IMOS20080114
    os<<"  HIR_filename      ="<<HIR_filename        <<endl;//IMOS20080114

    os<<"  fragmentation     ="<<fragmentation       <<endl;
    os<<"  momentum_losses   ="<<momentum_losses     <<endl;
    os<<"  radioactive_decay ="<<radioactive_decay   <<endl;
    os<<"  K_capture         ="<<K_capture           <<endl;
    os<<"  ionization_rate   ="<<ionization_rate     <<endl;  // IMOS20060420

    os<<"  ionization_losses ="<<ionization_losses   <<endl;
    os<<"  coulomb_losses    ="<<coulomb_losses      <<endl;
    os<<"  bremss_losses     ="<<bremss_losses       <<endl;
    os<<"  IC_losses         ="<<IC_losses           <<endl;
    os<<"  sync_losses       ="<<sync_losses         <<endl;

    os<<"  x_min    ="<<x_min   <<endl;
    os<<"  x_max    ="<<x_max   <<endl;
    os<<"  dx       ="<<dx      <<endl;
    os<<"  y_min    ="<<y_min   <<endl;
    os<<"  y_max    ="<<y_max   <<endl;
    os<<"  dy       ="<<dy      <<endl;

    os<<"  start_timestep   ="<< start_timestep<<endl;
    os<<"  end_timestep     ="<<   end_timestep<<endl;
    os<<"  timestep_factor  ="<<timestep_factor<<endl;
    os<<"  timestep_repeat  ="<<timestep_repeat<<endl;
    os<<"  timestep_repeat2 ="<<timestep_repeat2<<endl;
    os<<"  timestep_print   ="<<timestep_print <<endl;
    os<<"  timestep_diagnostics  ="<<timestep_diagnostics <<endl;
    os<<"  control_diagnostics  " << control_diagnostics <<endl;
    os<<"  solution_convergence " << solution_convergence<<endl;
    os<<"  solution_method      " << solution_method     <<endl;
    os<<"  solution_rel_accuracy" << solution_rel_accuracy<<endl;

    os<<"  network_iterations ="<<network_iterations<<endl;
    os<<"  network_iter_compl ="<<network_iter_compl<<endl;

    os<<"  prop_r   ="<<prop_r  <<endl;
    os<<"  prop_x   ="<<prop_x  <<endl;
    os<<"  prop_y   ="<<prop_y  <<endl;
    os<<"  prop_z   ="<<prop_z  <<endl;
    os<<"  prop_p   ="<<prop_p  <<endl;

    os<<"  use_symmetry  ="<<use_symmetry<<endl;
    os<<"  vectorized    ="<<vectorized  <<endl;

    //   os<<"  HI_survey             ="<<HI_survey            <<endl;
    //   os<<"  CO_survey             ="<<CO_survey            <<endl;

    os<<"  B_field_model         ="<<B_field_model        <<endl;                                                              //AWS20080313
    os<<"  B_field_name          ="<<B_field_name         <<endl;                                                              //AWS20080313
    os<<"  B_field_parameters    =";
    for ( int i=0;i<B_field_parameters.size();i++ ) {
       os<< B_field_parameters[i]<<"  ";
    }
    os<<endl;   //AWS20080313

    os<<"  ISRF_file             ="<<ISRF_file            <<endl;
    os<<"  ISRF_filetype         ="<<ISRF_filetype        <<endl;                                                              //AWS20091013
    os<<"  ISRF_healpixOrder     ="<<ISRF_healpixOrder    <<endl;                                                              //AWS20091013
    os<<"  ISRF_factors          ="<<ISRF_factors[0]<<","<<ISRF_factors[1]<<","<<ISRF_factors[2]   <<endl;

    if ( ISRF_filetype == 1 ) {FATAL(" ERROR: invalid ISRF_filetype, see README for details"); exit ( 0 );}                        //AWS20091013

    os<<"  proton_norm_Ekin              ="<<  proton_norm_Ekin   <<endl;
    os<<"  proton_norm_flux              ="<<  proton_norm_flux   <<endl;
    os<<"  proton_norm_type              ="<<  proton_norm_type   <<endl;
    os<<"  electron_norm_Ekin            ="<<electron_norm_Ekin   <<endl;
    os<<"  electron_norm_flux            ="<<electron_norm_flux   <<endl;
    os<<"  electron_norm_type            ="<<electron_norm_type   <<endl;

    os<<"  max_Z                 ="<<max_Z                <<endl;
    for ( int iZ=1; iZ<=max_Z; iZ++ ) os<<"use_Z_="<<iZ<<" "<<use_Z[iZ]<<endl;

    os<<"  total_cross_section  =" <<total_cross_section  <<endl;  //AWS20010620

    os<<"  cross_section_option  ="<<cross_section_option <<endl;
    os<<"  t_half_limit          ="<<t_half_limit         <<endl;  //AWS20010731


    os<<"  primary_electrons     ="<<primary_electrons    <<endl;
    os<<"  primary_positrons     ="<<primary_positrons    <<endl;
    os<<"  secondary_positrons   ="<<secondary_positrons  <<endl;
    os<<"  secondary_electrons   ="<<secondary_electrons  <<endl;
    os<<"  knock_on_electrons    ="<<knock_on_electrons   <<endl;  //IMOS20060504
    os<<"  secondary_antiprotons ="<<secondary_antiprotons<<endl;
    os<<"  tertiary_antiprotons  ="<<tertiary_antiprotons <<endl;  // IMOS20000605.5
    os<<"  secondary_protons     ="<<secondary_protons    <<endl;  // IMOS20000605.6
    os<<"  pair_production        ="<<pair_production       <<endl; // TAP20110402

    os<<"  gamma_rays            ="<<gamma_rays           <<endl;
    os<<"  pi0_decay             ="<<pi0_decay            <<endl;  // AWS20050218
    os<<"  IC_isotropic          ="<<IC_isotropic         <<endl;  // IMOS20060420
    os<<"  IC_anisotropic        ="<<IC_anisotropic       <<endl;
    os<<"  bremss                ="<<bremss               <<endl;  // IMOS20060420
    os<<"  synchrotron           ="<<synchrotron          <<endl;
    os<<"  free_free_absorption  ="<<free_free_absorption <<endl;  //AWS20110701
    os<<"  globalLuminosities    ="<<globalLuminosities   <<endl;

    os<<"  local_bubble_radius          ="<<local_bubble_radius         <<endl;
    os<<"  local_bubble_gas_fraction    ="<<local_bubble_gas_fraction   <<endl;
    os<<"  local_bubble_source_fraction ="<<local_bubble_source_fraction<<endl;

// DM: IMOS20050912
    os<<"  DM_positrons          ="<<DM_positrons         <<endl;
    os<<"  DM_electrons          ="<<DM_electrons         <<endl;
    os<<"  DM_antiprotons        ="<<DM_antiprotons       <<endl;
    os<<"  DM_gammas             ="<<DM_gammas            <<endl;

    os<<"  DM_double0            ="<<DM_double0           <<endl;
    os<<"  DM_double1            ="<<DM_double1           <<endl;
    os<<"  DM_double2            ="<<DM_double2           <<endl;
    os<<"  DM_double3            ="<<DM_double3           <<endl;
    os<<"  DM_double4            ="<<DM_double4           <<endl;
    os<<"  DM_double5            ="<<DM_double5           <<endl;
    os<<"  DM_double6            ="<<DM_double6           <<endl;
    os<<"  DM_double7            ="<<DM_double7           <<endl;
    os<<"  DM_double8            ="<<DM_double8           <<endl;
    os<<"  DM_double9            ="<<DM_double9           <<endl;

    os<<"  DM_int0               ="<<DM_int0              <<endl;
    os<<"  DM_int1               ="<<DM_int1              <<endl;
    os<<"  DM_int2               ="<<DM_int2              <<endl;
    os<<"  DM_int3               ="<<DM_int3              <<endl;
    os<<"  DM_int4               ="<<DM_int4              <<endl;
    os<<"  DM_int5               ="<<DM_int5              <<endl;
    os<<"  DM_int6               ="<<DM_int6              <<endl;
    os<<"  DM_int7               ="<<DM_int7              <<endl;
    os<<"  DM_int8               ="<<DM_int8              <<endl;
    os<<"  DM_int9               ="<<DM_int9              <<endl;

    os<<"  skymap_format        ="<<skymap_format         <<endl;
    os<<"  output_gcr_full      ="<<output_gcr_full       <<endl;
    os<<"  warm_start           ="<<warm_start            <<endl;

    os<<"  verbose    ="<<verbose    <<endl;
    os<<"  test_suite ="<<test_suite <<endl;

    INFO(os.str());
    //   exit(0);
}

void Galdef::AssignDefaultParameters () {

    n_spatial_dimensions = 2;

    r_min =   0.0;
    r_max =  20.0;
    dr    =   1.0;

    x_min = -15.0;
    x_max = +15.0;
    dx    =   1.0;

    y_min = -15.0;
    y_max = +15.0;
    dy    =   1.0;

    z_min =  -4.0;
    z_max =  +4.0;
    dz    =   0.2;

    p_min =   10.0;
    p_max =   1.e5;
    p_factor= 1.3 ;

    Ekin_min             =1.0e1;
    Ekin_max             =1.0e8;
    Ekin_factor          =1.4;

    p_Ekin_grid = "Ekin";

    E_gamma_min          = 1.e1;
    E_gamma_max          = 1.e6;
    E_gamma_factor       = 1.e1;
    integration_mode     = 1;

    nu_synch_min         = 1.0e6;
    nu_synch_max         = 1.0e10;
    nu_synch_factor      = 2.0;

    long_min             =   .0;
    long_max             =360.0;
    lat_min              =-90.0;
    lat_max              =+90.0;
    d_long               =  0.5 ;
    d_lat                =  0.5 ;
    lat_substep_number   =  1;
    LoS_step             =  0.01;
    LoS_minStep          =  1e-4;
    LoS_substep_number   =  1;
    los_integration_mode =  0;
    los_integration_accuracy =  0;
    anisoHealpixOrder    =  2;

    D0_xx                =5.75e28;
    D_rigid_ref          =4.0e3;
    D_rigid_br           =4.0e3;
    D_g_1                = 0.34;
    D_g_2                = 0.34;
    D_eta                = 1.0;
    Dxx_plane_scale      = -1.0;
    Dxx_plane_scale_height = 0.2;
    B_dep_diffusion      = 0;
    D_xx_max             = 1e30;
    D_xx_min             = 1e28;
    diff_reacc           =1 ;
    v_Alfven             =36.;

    damping_p0           = 1.e6;
    damping_const_G      = 0.02;
    damping_const_K      = 1.e5;
    damping_max_path_L   = 3.e21;

    convection           =0 ;
    v0_conv              =0.;
    dvdz_conv            =10.;

    rigid_min            =0;
    rigid_max            =numeric_limits<double>::max();
    inj_Ekin_min         =0;
    inj_Ekin_max         =numeric_limits<double>::max();

    He_H_ratio           = 0.11;
    n_X_CO               = 0;
    X_CO                 = 1.9E20;

    propagation_X_CO     = 0;
    nHI_model            = 2;
    nH2_model            = 2;
    nHII_model           = 1;
    HII_Te               = 7000.;//AWS20110701
    HII_clumping_factor  = 1.0  ;//AWS20110701

    COR_filename         = "";
    HIR_filename         = "";

    B_field_model        = 050100020;
    B_field_name         = "galprop_original";
    //    B_field_parameters.resize(10); //AWS20110520 no such restriction needed
    for ( int i = 0; i < B_field_parameters.size(); ++i ) {
       B_field_parameters[i] = -i-1;
    }

    ISRF_file            = "";

    ISRF_filetype        = 0;
    ISRF_factors.resize(3,1.0);
//  ISRF_factors[0]      = 1.0;
//  ISRF_factors[1]      = 1.0;
//  ISRF_factors[2]      = 1.0;

    fragmentation        =1;
    momentum_losses      =1;
    radioactive_decay    =1;
    K_capture            =1;
    ionization_rate      =0;

    ionization_losses    = true;
    coulomb_losses       = true;
    bremss_losses        = true;
    IC_losses            = true;
    sync_losses          = true;

    start_timestep       =  1.e9;
    end_timestep         =  1.e2;
    timestep_factor      =  0.25;
    timestep_repeat      = 20;
    timestep_repeat2     = 0;
    timestep_print       =10000;
    timestep_diagnostics =10000;
    control_diagnostics  =0;
    solution_convergence =0; 
    solution_method      =4; 
    solution_rel_accuracy=1e-3; 

    network_iterations   = 1;
    network_iter_compl   = 1;

    prop_r               = 1;
    prop_x               = 1;
    prop_y               = 1;
    prop_z               = 1;
    prop_p               = 1;

    use_symmetry         = 0;

    vectorized           = 0;
    //g.source_parameters[1] = 0.5;
    //g.source_parameters[2] = 1.0;
    //g.source_parameters[3] = 20.0;

    proton_norm_Ekin     = 1.00e+5;
    proton_norm_flux     = 4.90e-9;
    proton_norm_type     = 1;           //AWS20101202

    electron_norm_Ekin   = 34.5e3;
    electron_norm_flux   = .40e-9;
    electron_norm_type   = 1;           //AWS20101202

    max_Z                = 2;
//  use_Z[0]             = 1;
//  use_Z[1]             = 1;

//  isotopic_abundance[1][1]= 1.06e+06;
//  isotopic_abundance[1][2]=     34.8;
//  isotopic_abundance[2][3]=    9.033;
//  isotopic_abundance[2][4]= 7.199e+04;

    total_cross_section  = 1;
    cross_section_option = 011;

    t_half_limit         = 1.0e4;

    primary_electrons    = 1;
    primary_positrons    = 0;
    secondary_positrons  = 0;
    secondary_electrons  = 0;
    knock_on_electrons   = 0;
    secondary_antiprotons= 0;
    tertiary_antiprotons = 0;
    secondary_protons    = 0;
    
    pair_production       = 0;

    gamma_rays           = 0;
    pi0_decay            = 1;
    IC_isotropic         = 1;
    IC_anisotropic       = 0;
    synchrotron          = 1;
    free_free_absorption = 0;
    bremss               = 1;
    globalLuminosities   = 0;

    local_bubble_radius          = 0;
    local_bubble_gas_fraction    = 0;
    local_bubble_source_fraction = 0;

    DM_positrons         = 0;
    DM_electrons         = 0;
    DM_antiprotons       = 0;
    DM_gammas            = 0;

    fCameraLocation.resize(3);
    fCameraLocation[0] = Rsun;
    fCameraLocation[1] = 0;
    fCameraLocation[2] = 0;

    skymap_format        = 0;
    output_gcr_full      = 0;
    warm_start           = 0;

    verbose              = 0;
    test_suite           = 0;

}










