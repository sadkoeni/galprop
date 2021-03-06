#include "galprop_classes.h"
#include "galprop_internal.h"

//#include <fort_interface.h>

#include <ErrorLogger.h>
#include <Timer.h>
#include <Nuclei_Interface.h>
#include <Processes_Interface.h>

#include <config.h>

#include <cassert>
#include <iostream>
#include <sstream>

//using namespace std;

Galprop* gGalprop = nullptr;

Galprop::Galprop() {

  gGalprop = this;
  gcr = 0;
  //WarpArmModelFunction = 0;

}

Galprop::~Galprop() {

   delete[] gcr;
   //delete WarpArmModelFunction;

}

int Galprop::Run(std::string galdefPath,
		 std::string fitsPath,
		 std::string outputPath,
		 std::string outputPrefix,
		 std::string runNumber) {

  if (configure.init(galdefPath, fitsPath, outputPath, outputPrefix)) {

    FATAL("Internal error. Fix data paths!");
    return 1;

  }

  if (galdef.read(configure.fVersion, runNumber, configure.fGaldefDirectory)) {

	std::cout<<"Pos1"<<std::endl;
    FATAL("Internal error. Problem reading from galdef file!");
	std::cout<<"Pos2"<<std::endl;
    return 1;

  } 
  
  assert(2 == galdef.n_spatial_dimensions || 3 == galdef.n_spatial_dimensions);

  // electron propagation test vs analytical formula; assign galdef parameters IMOS20061030

  if (99 == abs(galdef.DM_int0)) {

    std::cout<<std::endl<<" >>>> Running electron propagation test"<<std::endl<<std::endl;
    galdef.n_spatial_dimensions =2;
    galdef.D_rigid_br           =1.e3;
    galdef.diff_reacc           =0;
    galdef.convection           =0;
    galdef.momentum_losses      =1;
    galdef.max_Z                =1;
    galdef.primary_electrons    =0;
    galdef.secondary_positrons  =0;
    galdef.secondary_electrons  =0;
    galdef.knock_on_electrons   =0;
    galdef.secondary_antiprotons=0;
    galdef.tertiary_antiprotons =0;
    galdef.gamma_rays           =0;
    galdef.DM_positrons         =1; //will contain analytical solution
    galdef.DM_electrons         =1; //will contain numerical  solution
    galdef.DM_antiprotons       =0;
    galdef.DM_gammas            =0;
    std::cout<<"** Test values are re-assigned in Galprop.cc: **"      <<std::endl;
    std::cout<<"  DM_int0               "<<galdef.DM_int0              <<std::endl;
    std::cout<<"  DM_double6            "<<galdef.DM_double6           <<std::endl;
    std::cout<<"  DM_double7            "<<galdef.DM_double7           <<std::endl;
    std::cout<<"  DM_double8            "<<galdef.DM_double8           <<std::endl;
    std::cout<<"  DM_double9            "<<galdef.DM_double9           <<std::endl;
    std::cout<<"  DM_positrons          "<<galdef.DM_positrons         <<std::endl;
    std::cout<<"  DM_electrons          "<<galdef.DM_electrons         <<std::endl;
    std::cout<<"  DM_antiprotons        "<<galdef.DM_antiprotons       <<std::endl;
    std::cout<<"  DM_gammas             "<<galdef.DM_gammas            <<std::endl;
    std::cout<<"  D_rigid_br            "<<galdef.D_rigid_br           <<std::endl;
    std::cout<<"  diff_reacc            "<<galdef.diff_reacc           <<std::endl;
    std::cout<<"  convection            "<<galdef.convection           <<std::endl;
    std::cout<<"  momentum_losses       "<<galdef.momentum_losses      <<std::endl;
    std::cout<<"  max_Z                 "<<galdef.max_Z                <<std::endl;
    std::cout<<"  primary_electrons     "<<galdef.primary_electrons    <<std::endl;
    std::cout<<"  secondary_positrons   "<<galdef.secondary_positrons  <<std::endl;
    std::cout<<"  secondary_electrons   "<<galdef.secondary_electrons  <<std::endl;
    std::cout<<"  knock_on_electrons    "<<galdef.knock_on_electrons   <<std::endl;
    std::cout<<"  secondary_antiproton  "<<galdef.secondary_antiprotons<<std::endl;
    std::cout<<"  tertiary_antiproton   "<<galdef.tertiary_antiprotons <<std::endl;
    std::cout<<"  gamma_rays            "<<galdef.gamma_rays           <<std::endl<<std::endl;
    std::cout<<"  n_spatial_dimensions  "<<galdef.n_spatial_dimensions <<std::endl;
  
  }
    
  // global test
  if (0 < galdef.test_suite) { 

    //test_suite(); 
    return 0; 

  }
  
  //reading all isotopic cross sections, nuclear reaction network, nuclear fits
	std::cout<<"Pos2"<<std::endl;
  read_nucdata(configure.fGaltoolslibDataPath);
  
  //initialization of the Barashenkov & Polanski cross section code 
	std::cout<<"Pos3"<<std::endl;
  processes::sigtap_cc(-1, configure.fGaltoolslibDataPath);             // IMOS20010511 AWS20010620
  
  //initialization of the Webber's routine
	std::cout<<"Pos4"<<std::endl;
  nuclei::set_sigma_cc(configure.fGaltoolslibDataPath);

  //initialization of the [Kachelriess, Moskalenko & Ostapchenko, arXiv:1502.04158] cross sections
	std::cout<<"Pos5"<<std::endl;
  processes::aprtab_cc(configure.fGaltoolslibDataPath);

  utl::Timer::initialize();

  int stat;
	std::cout<<"Pos6"<<std::endl;
  TIME_FUNCTION(stat,create_galaxy);
	std::cout<<"Pos7"<<std::endl;
  if (0 != stat) {
     FATAL("Error when creating galaxy, aborting");
     return 1;
  }

  TIME_FUNCTION(stat,create_gcr);
  if (0 != stat) {
     FATAL("Error when creating gcr, aborting");
     return 1;
  }

  // major routine
  TIME_FUNCTION(stat,propagate_particles);
  if (0 != stat) {
     FATAL("Error when propagating particles, aborting");
     return 1;
  }

  //cr_luminosity();
  
  //TIME_FUNCTION(stat,store_gcr);
  //if (0 != stat) 
  // return 1;                                      //GJ20100104
  //deleting all cross section arrays etc.
  cleanup_nucdata();
  
  //   if(print_BC() !=0) return 1;
  //   if(store_gcr() !=0) return 1; //IMOS20030129
  
  //if (galdef.output_gcr_full) {

  //TIME_FUNCTION(stat,store_gcr_full);
  //if (0 != stat) 
  //  return 1;
  
  //}

  //exit(0);

  // Movable viewing position for skymap generation

  valarray<double> cameraLocation(3);

  cameraLocation[0] = galdef.fCameraLocation[0]; 
  cameraLocation[1] = galdef.fCameraLocation[1];
  cameraLocation[2] = galdef.fCameraLocation[2];

  // Fill energy/frequency grids

  //Moved from create_galaxy   Gulli20070810
  /*galaxy.E_gamma_min = galdef.E_gamma_min;
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
  */
  // Generate emissivities first

  // Initialise dimensions of emissivity and ionisation rate arrays -- the 
  // emissivity arrays are always initialised independent of whether we are
  // calculating for a particular process or not because for the gamma rays
  // we want to guarantee that all components will be available (some maybe
  // zero because they are not calculated) for the pair absorption source 
  // function. 
  // Profiling indicates that the emissivity arrays have a minimal memory
  // footprint overhead, even in 3D, so the persistence of these arrays is 
  // not a resource issue.

  if (2 == galdef.n_spatial_dimensions) {

    //if (galdef.pi0_decay)
    galaxy.pi0_decay_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);  
    
    //if (galdef.bremss) {
    
    galaxy.bremss_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);
      
    galaxy.bremss_ionized_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);
      
    //}
    
    //if (galdef.IC_isotropic || galdef.IC_anisotropic) {
    
    //Moved from create_galaxy  //Gulli20070810
    galaxy.IC_iso_emiss = new Distribution[galaxy.n_ISRF_components];
    
    for (int i = 0; i < galaxy.n_ISRF_components; ++i) {
      
      galaxy.IC_iso_emiss[i].init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);
      
    }
    
    //}
    
    //if (galdef.synchrotron) {
    
    galaxy.synchrotron_emiss.init  (galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);
    galaxy.synchrotron_Q_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);                  //AWS20100708
    galaxy.synchrotron_U_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);                  //AWS20100708

    galaxy.free_free_emiss.init    (galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);                  //AWS20110905
    
    //}
    
    //if (galdef.DM_gammas) {                           // IMOS20050912
    
    galaxy.DM_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid); 
    
    //}
    
    if (galdef.ionization_rate) 
      galaxy.ionization_rate.init(galaxy.n_rgrid, galaxy.n_zgrid, 1);
    
  }

  if (3 == galdef.n_spatial_dimensions) {

    //if (galdef.pi0_decay)
    galaxy.pi0_decay_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);  
    
    //if (galdef.bremss) {
    
    galaxy.bremss_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);
    
    galaxy.bremss_ionized_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);
    
    //}
    
    //if (galdef.IC_isotropic || galdef.IC_anisotropic) {
      
    //Moved from create_galaxy  //Gulli20070810
    galaxy.IC_iso_emiss = new Distribution[galaxy.n_ISRF_components];
    
    for (int i = 0; i < galaxy.n_ISRF_components; ++i) {
      
      galaxy.IC_iso_emiss[i].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);
      
    }
      
    //}
    
    //if (galdef.synchrotron) {
      
    galaxy.synchrotron_emiss  .init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);
    galaxy.synchrotron_Q_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);                  //AWS20100708
    galaxy.synchrotron_U_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);                  //AWS20100708

    galaxy.free_free_emiss    .init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);                  //AWS20110905
    
    //}
  
    //if (galdef.DM_gammas) {                           // IMOS20050912
    
    galaxy.DM_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid); 
    
    //}
   
    if (galdef.ionization_rate) 
      galaxy.ionization_rate.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, 1);
    
  }

  // Generate emissivities

  if (galdef.pi0_decay) {
    
    TIME_FUNCTION(stat, gen_pi0_decay_emiss);
    //if (stat != 0) return 1;
    TIME_FUNCTION(stat, store_pi0_decay_emiss);
    //if (stat != 0) return 1;
    
  }
  
  if (galdef.bremss) {
    
    TIME_FUNCTION(stat, gen_bremss_emiss);
    //if (stat != 0) return 1;
    //TIME_FUNCTION(stat, store_bremss_emiss);
    //if (stat != 0) return 1;
        
  }
  
  if (galdef.IC_isotropic || galdef.IC_anisotropic) {
    
    auto isGood = false;
    
    TIME_FUNCTION(isGood, gen_IC_emiss); 
    if (!isGood) {
       FATAL("Error when generating IC emission, aborting");
       return 1;
    }
    
    //const std::string icType = "isotropic";
    
    //TIME_SUBROUTINE(store_IC_emiss);
    
  }

  if (galdef.DM_gammas) {

    TIME_FUNCTION(stat, gen_DM_emiss);
    //if (stat != 0) return 1;
    TIME_FUNCTION(stat, store_DM_emiss);
    //if (stat != 0) return 1;

  }

  valarray<double> pairProductionLuminosity(0., galaxy.n_pgrid);

  if (galdef.pair_production) {

    INFO("Start calculating pair production contribution");
        
    Particle particle;

    auto found = false;

    for (int i = 0; i < n_species; ++i) {

      if ("secondary_positrons" == gcr[i].name ||
	  "secondary_electrons" == gcr[i].name) {
	
	particle = gcr[i];
	particle.name = "secondary_positrons_pair"; // Can be either -- we'll be taking half of the propagated solution and adding to secondary electrons and positrons
	found = true;

	//std::cout << "found" << std::endl;
	
	//for (int iR = 0; iR < 1; ++iR)
	//for (int iZ = galaxy.n_zgrid/2; iZ < galaxy.n_zgrid/2+1; ++iZ)
	//for (int iE = 0; iE < galaxy.n_pgrid; ++iE)
	//  std::cout << "sec1: " << iR << " " << iZ << " " << iE << " " << galaxy.p[iE] << " " << particle.cr_density.d2[iR][iZ].s[iE] << std::endl;


	particle.cr_density = 0;
	particle.delete_transport_arrays();
	//particle.print();
	break;

      }

    }

    if (found) {

      particle.create_transport_arrays();
      fill_transport_arrays(particle);
      particle.primary_source_function = 0; // It gets filled in fill_transport_arrays ...

      TIME_SUBROUTINE(gen_secondary_pair_source, particle);

      // If we're computing the global luminosities (later), we need to 
      // obtain _now_ the global luminosity for the pair-produced secondary
      // e+/-. Otherwise, these source functions will be destroyed at exit
      // from this control block, and will not be available later on. When we
      // rework this code we should think of a better way to do this ...

      if (2 == galaxy.n_spatial_dimensions) {

	for (int iR = 0; iR < galaxy.n_rgrid; ++iR)
	  for (int iZ = 0; iZ < galaxy.n_zgrid; ++iZ)
	    for (int iP = 0; iP < galaxy.n_pgrid; ++iP) {
	      pairProductionLuminosity[iP] += particle.secondary_source_function.d2[iR][iZ].s[iP]*2.*utl::kPi*galaxy.r[iR]*galaxy.dr*galaxy.dz;
	      //std::cout << "PP: " << iR << " " << iZ << " " << iP << " " << particle.secondary_source_function.d2[iR][iZ].s[iP] << " " << pairProductionLuminosity[iP] << std::endl;
	    }

      }

      if (3 == galaxy.n_spatial_dimensions) {

	for (int iX = 0; iX < galaxy.n_xgrid; ++iX)
	  for (int iY = 0; iY < galaxy.n_ygrid; ++iY)
	    for (int iZ = 0; iZ < galaxy.n_zgrid; ++iZ)
	      for (int iP = 0; iP < galaxy.n_pgrid; ++iP)
		pairProductionLuminosity[iP] += particle.secondary_source_function.d3[iX][iY][iZ].s[iP]*galaxy.dx*galaxy.dy*galaxy.dz;

      }

      pairProductionLuminosity *= 2.*4.*utl::kPi*1./(utl::kSpeedOfLight_SI*utl::m/utl::cm)*pow(utl::kpc/utl::cm, 3.); // Conversion MeV^-1 cm^-2 s^-2 sr^-1 -> MeV^-1 s^-1 and the extra 2 in front is to account for both e+ and e-

      pairProductionLuminosity *= MeV_to_erg*pow(particle.Ekin, 2.); // Conversion MeV^-1 s^-1 -> MeV s^-1 -> erg s^-1

      propel(particle);
    
      //for (int iR = 0; iR < 1; ++iR)
      //for (int iZ = galaxy.n_zgrid/2; iZ < galaxy.n_zgrid/2+1; ++iZ)
      //  for (int iE = 0; iE < galaxy.n_pgrid; ++iE)
      //    std::cout << "sec2: " << iR << " " << iZ << " " << iE << " " << galaxy.p[iE] << " " << particle.cr_density.d2[iR][iZ].s[iE] << " " << particle.secondary_source_function.d2[iR][iZ].s[iE] << std::endl;

      //exit(0);
  
      // Find secondary positrons and electrons, add half of the propagated
      // particle flux from above to each particle type.

      for (int i = 0; i < n_species; ++i) {

	if ("secondary_positrons" == gcr[i].name) {

	  gcr[i].cr_density += particle.cr_density; // Cross section is only for a single species

	}

	if ("secondary_electrons" == gcr[i].name) {

	  gcr[i].cr_density += particle.cr_density; // Cross section is only for a single species

	}

      }

      // Rerun IC and brem calculations with the additional
      // contribution to the secondary electrons and positrons

      if (galdef.bremss)     
	TIME_FUNCTION(stat, gen_bremss_emiss);
      
      if (galdef.IC_isotropic || galdef.IC_anisotropic) {
    
	TIME_SUBROUTINE(gen_IC_emiss); 

      }

      INFO("Finished calculating pair production contribution");

    } else {

      INFO("No electrons found for propagation calculation. Specify any of primary electrons, or secondary electrons/positrons to calculate pair absorption contribution.");

    }

  }

  // Store emissivities of IC and Brem, if calculated

  if (galdef.bremss) {
  
    TIME_FUNCTION(stat, store_bremss_emiss);
    //if (stat != 0) return 1;
        
  }
  
  if (galdef.IC_isotropic || galdef.IC_anisotropic) {
        
    TIME_SUBROUTINE(store_IC_emiss);
    
  }

  // Continue with calculation of other components

  if (galdef.synchrotron) {

    TIME_FUNCTION(stat, gen_synch_emiss);
    //if (stat != 0) return 1;
    TIME_FUNCTION(stat, store_synch_emiss);
    //if (stat != 0) return 1;       //AWS20080314
  
  }

  // Generate ionisation rate

  if (galdef.ionization_rate) {

    TIME_FUNCTION(stat, gen_ionization_rate);
    //if (stat != 0) return 1;
    TIME_FUNCTION(stat, store_ionization_rate);
    //if (stat != 0) return 1;
      
    //Delete Arrays  Gulli20071003
    
    galaxy.ionization_rate.delete_array();
    
  }

  // Generate luminosities
  
  if (galdef.globalLuminosities) {

    TIME_FUNCTION(stat, gen_luminosity);
    //if (stat != 0) return 1;

    //for (int i = 0; i < galaxy.n_pgrid; ++i)
    //std::cout << i << " " << galaxy.p[i] << " " << galaxy.fProtonLuminosity[i] << " " << galaxy.fHeliumLuminosity[i] << " " << galaxy.fNucleiLuminosity[i] << " " << galaxy.fPrimaryElectronLuminosity[i] << " " << galaxy.fSecondaryElectronLuminosity[i] << " " << pairProductionLuminosity[i] << std::endl;

    galaxy.fSecondaryElectronLuminosity += pairProductionLuminosity; // This is zero if not calculated, so a trivial addition in that case

    TIME_SUBROUTINE(store_luminosity); 
    
    /*if (-999 == galdef.verbose) {

      for (int iP = 0; iP < galaxy.n_pgrid; ++iP) {
	
	std::cout << iP << " " << galaxy.p[iP] << " " << galaxy.fProtonLuminosity[iP] << " " << galaxy.fHeliumLuminosity[iP] << " " << galaxy.fNucleiLuminosity[iP] << " " << galaxy.fPrimaryElectronLuminosity[iP] << " " << galaxy.fSecondaryElectronLuminosity[iP] << " " << galaxy.fSecondaryElectronLuminosity[iP] - pairProductionLuminosity[iP] << " " << pairProductionLuminosity[iP] << std::endl;

      }

      for (int iP = 0; iP < galaxy.nu_ISRF.size(); ++iP) {

	std::cout << iP << " " << galaxy.nu_ISRF[iP]*h_planck*erg_to_eV << " " << galaxy.fISRFLuminosity[iP] << std::endl;

      }

      for (int iP = 0; iP < galaxy.n_E_gammagrid; ++iP) {

	std::cout << iP << " " << galaxy.E_gamma[iP] << " " << galaxy.fPi0DecayLuminosity[iP] << " " << galaxy.fBremsstrahlungLuminosity[iP] << " " << galaxy.fICLuminosity[iP] << " " << galaxy.fDMLuminosity[iP] << std::endl;

      }

      exit(0);

    }
    */
    
  }
  
  // Store CR data

  TIME_FUNCTION(stat,store_gcr);
  //if (0 != stat) 
  // return 1;  

  if (galdef.output_gcr_full) {
    
    TIME_FUNCTION(stat,store_gcr_full);
    //if (0 != stat) 
    //  return 1;
  
  }

  // Generate skymaps second

  // Generate skymaps for camera location(s) provided -- split according 
  // to gas and IC/sync. We could do this in one pass but the memory 
  // consumption can be fairly large if we read in the gas maps for 
  // a high healpix order *and* want to do anisotropic IC at the same time. 
  // So we create the gas maps in-memory first, do the LOS integrations 
  // for those components, delete the gas maps, etc., and then move on to
  // the IC/sync
  
  // One step integration function with new integration method if selected
  //if (1 == galdef.los_integration_mode) {

  TIME_SUBROUTINE(gen_skymaps);

  galaxy.pi0_decay_emiss.delete_array();

  galaxy.bremss_emiss.delete_array();
  galaxy.bremss_ionized_emiss.delete_array();
 
  for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
	
    galaxy.IC_iso_emiss[icomp].delete_array();
    
  }

  galaxy.synchrotron_emiss.delete_array();
  galaxy.synchrotron_Q_emiss.delete_array(); //AWS20100708
  galaxy.synchrotron_U_emiss.delete_array(); //AWS20100708
  
  galaxy.DM_emiss.delete_array();  //Gulli20070810
  
  std::ostringstream bufProc;
  bufProc << "Completed processing galdef_" << configure.fVersion << "_" << runNumber;
  INFO(bufProc.str());

  INFO("Exit: Galprop main procedure");
  
  return 0;

}

