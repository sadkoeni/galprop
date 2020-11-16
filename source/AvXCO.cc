#include "galprop_classes.h"
#include "galprop_internal.h"

#include <ErrorLogger.h>
#include <BaseSkyFitsIO.h>

#include <iostream>
#include <fstream>


int Galprop::AvXCO(std::string galdefPath,
		   std::string fitsPath,
		   std::string outputPath,
		   std::string outputPrefix,
		   std::string runNumber) {

  if (configure.init(galdefPath, fitsPath, outputPath, outputPrefix)) {

    FATAL("Internal error. Fix data paths!");
    return 1;

  }

  if (galdef.read(configure.fVersion, runNumber, configure.fGaldefDirectory)) {

    FATAL("Internal error. Problem reading from galdef file!");
    return 1;

  } 

	/*
  if ( 0 != create_galaxy() ) {

     FATAL("Internal error. Problem allocating memory.");
     return 1;

  }
	*/

  //Read in the CO gas maps
  read_gas_maps("COR");

  //Set up the los integrators
  std::vector< std::unique_ptr< SM::LOSfunction<double> > > funcs(3);
  funcs[0].reset(new GasFunction("CO", 0, *this)); //Plain CO
  funcs[1].reset(new GasFunction("H2", 0, *this)); //H2, X corrected CO
  funcs[2].reset(new GasFunction("H2", 1.0, *this)); //H2, scaled with radius to get effective radius

  std::vector<double> Rbins(&galaxy.R_bins[0], &galaxy.R_bins[0]+galaxy.n_Ring);
  Rbins.push_back(galaxy.R_bins[2*galaxy.n_Ring-1]);
  SM::LOSintegrator<double> losInt(galdef.r_max, galdef.z_min, galdef.z_max, Rbins, galdef.fCameraLocation, galdef.LoS_step, galdef.los_integration_accuracy, galdef.LoS_minStep);

  if (3 == galdef.skymap_format) {
     std::vector< std::unique_ptr< SM::BaseSky<double> > > avXCOhp(galaxy.n_Ring), avRhp(galaxy.n_Ring);

     for (int i = 0; i < galaxy.n_Ring; ++i ) {
        avXCOhp[i] = galaxy.hpCOR[i]->clone();
        avRhp[i] = galaxy.hpCOR[i]->clone();
     }

#pragma omp parallel for schedule(dynamic) default(shared)
     for ( int ii = 0; ii < galaxy.hpCOR[0]->Npix(); ++ii ) {
        auto co = galaxy.hpCOR[0]->GetCoordinate(ii);
        double l, b;
        co.getCoordinates(l, b, SM::CoordSys::GAL);
        l *= utl::kConvertRadiansToDegrees;
        b *= utl::kConvertRadiansToDegrees;

	double dl = 90./galaxy.hpCOR[0]->Nside();
	vector< vector<double> > gas = losInt.integrate(l, b, funcs);

	for ( size_t j = 0; j < gas[0].size(); ++j ) {
	   if ( gas[0][j] != 0 ) {
#pragma omp critical (avXCO)
              {
                 avXCOhp[j]->SetValue(ii, 0, gas[1][j]/gas[0][j]);
                 avRhp  [j]->SetValue(ii, 0, gas[2][j]/gas[1][j]);
              }
	   }
	}
     }

     const std::string filestart = configure.fOutputDirectory + configure.fOutputPrefix;
     const std::string fileend = "_healpix_" + galdef.galdef_ID + ".gz";

     for ( int i = 0; i < galaxy.n_Ring; ++i ) {
        std::ostringstream ost;
        ost << filestart << "averageXCO_ring_"<<i+1<<fileend;
        SM::writeToFits(*avXCOhp[i], ost.str(), true, true, "Radius", "kpc");
        ost.str("");

        ost << filestart << "averageXCORadius_ring_"<<i+1<<fileend;
        SM::writeToFits(*avRhp[i], ost.str(), true, true, "Radius", "kpc");
     }

     // Calculate the average XCO and the corresponding radius, weighted by the
     // real CO map.
     std::vector<double> avXCO(galaxy.n_Ring, 0.0), avR(galaxy.n_Ring, 0.0);
     for (size_t i_ring = 0; i_ring < galaxy.n_Ring; ++i_ring) {
	double weight(0), weightR(0);
        for ( auto it = galaxy.hpCOR[i_ring]->begin(); it != galaxy.hpCOR[i_ring]->end(); ++it) {
           const size_t ihp = it.healpixIndex();

           avXCO[i_ring] += avXCOhp[i_ring]->GetValue(ihp,0) * (*it);
           avR  [i_ring] += avRhp  [i_ring]->GetValue(ihp,0)*avXCOhp[i_ring]->GetValue(ihp,0) * (*it);
           weight += *it;
           weightR += avXCOhp[i_ring]->GetValue(ihp,0) * (*it);

	}
	avXCO[i_ring] /= weight;
	avR  [i_ring] /= weightR;
     }

     // Possibly store it as a fits file in the future, but a simple text file
     // will due for now
     const std::string filename = filestart + "averageXCO_" + galdef.galdef_ID + ".txt";
     std::fstream ofs(filename.c_str(), std::fstream::trunc | std::fstream::out);
     for (size_t i_ring = 0; i_ring < galaxy.n_Ring; ++i_ring) {
	ofs << avR[i_ring] << "\t" << avXCO[i_ring] << std::endl;
     }
     ofs.close();
  }

  return 0;
} 
