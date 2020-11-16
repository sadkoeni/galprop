#include "Galprop.h"
#include "galprop_classes.h"
#include "galprop_internal.h"
#include "constants.h"
#include "fitsio.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include <valarray>
#include <vector>

#include <PhysicalConstants.h>
#include <ErrorLogger.h>

#include <RadiationField.h>

// interfaces to nH.cc
double nHI(double Rkpc, double Zkpc);
double nH2(double Rkpc, double Zkpc);
double nH2(double Rkpc, double Zkpc, double H2toCO);
double nHII(double Rkpc, double Zkpc);

// Compute luminosity of Galaxy and other integral properties like gas mass
// CR density  gcr.cr_density is in c/4pi * n(p) [ cm s^-1  * cm^-3 MeV^-1]
// Primary source function is    in c/4pi        [ cm s^-1  * cm^-3 MeV^-1 s-1]
// Luminosity (cm^-3)= Integral n(Ekin). A.Ekin. Ekin dlog(Ekin)
// because each nucleus has KE=A.Ekin
//
// The particle spectra are assumed to be on equal kinetic energy per 
// nucleon grids which is the standard for galprop.
// UNITS OF c/4pi*n(p) = (1/A)flux(Ekin) = (1/A)c/4pi*beta*n(Ekin)
// so factor required = A.4pi/c/beta 

int Galprop::gen_luminosity() {

  INFO("Entry");

  galaxy.fProtonLuminosity.resize(galaxy.n_pgrid);
  galaxy.fHeliumLuminosity.resize(galaxy.n_pgrid);
  galaxy.fNucleiLuminosity.resize(galaxy.n_pgrid);
  galaxy.fPrimaryElectronLuminosity.resize(galaxy.n_pgrid); 
  galaxy.fSecondaryElectronLuminosity.resize(galaxy.n_pgrid); // e+/- total
  
  galaxy.fProtonLuminosity = 0;
  galaxy.fHeliumLuminosity = 0;
  galaxy.fNucleiLuminosity = 0;
  galaxy.fPrimaryElectronLuminosity = 0;
  galaxy.fSecondaryElectronLuminosity = 0;

  galaxy.fSynchrotronLuminosity.resize(galaxy.n_nu_synchgrid);
  galaxy.fICLuminosity.resize(galaxy.n_E_gammagrid);
  galaxy.fPi0DecayLuminosity.resize(galaxy.n_E_gammagrid);
  galaxy.fBremsstrahlungLuminosity.resize(galaxy.n_E_gammagrid);
  galaxy.fDMLuminosity.resize(galaxy.n_E_gammagrid);
  galaxy.fISRFLuminosity.resize(galaxy.nu_ISRF.size());

  galaxy.fSynchrotronLuminosity = 0;
  galaxy.fICLuminosity = 0;
  galaxy.fPi0DecayLuminosity = 0;
  galaxy.fBremsstrahlungLuminosity = 0;
  galaxy.fDMLuminosity = 0;
  galaxy.fISRFLuminosity = 0;

  // ISRF luminosity is calculated if the ISRF filetype is version 3 or later.
  // Otherwise it will be zero. This is not a bug, but by design (we no longer
  // support the earlier ISRF versions).

  if (3 == galdef.ISRF_filetype) {
    
    // Construct the radiation field and integrate over the surface of the
    // bounding volume for the total luminosity emitted by the galaxy
    
    const std::string fitsDirectory = configure.fFITSDataDirectory;
    const std::string isrfFilename = galdef.ISRF_file;
    const std::string filename = fitsDirectory + isrfFilename;
    
    ostringstream buf1;
    buf1 << "Reading ISRF from " << filename;
    INFO(buf1.str());
    
    if ( filename != galaxy.fISRFloadedfile ) {
       delete galaxy.fISRFrf;
       galaxy.fISRFrf = new rf::RadiationField(filename, galaxy.nu_ISRF, galdef.ISRF_healpixOrder);
       galaxy.fISRFloadedfile = filename;
    }

    rf::RadiationField &rf = *galaxy.fISRFrf;
    
    valarray<double> targetE(0., galaxy.nu_ISRF.size());
    targetE = h_planck*erg_to_eV*galaxy.nu_ISRF;

    if (2 == galaxy.n_spatial_dimensions) {

      const valarray<double>& rBound = rf.GetBoundaryR();
      const valarray<double>& zBound = rf.GetBoundaryZ();
      
      const double dR = galaxy.dr;
      const double dZ = galaxy.dz;
      
      const double rMax = std::min(galaxy.r_max, rBound[rBound.size()-1]);// : std::min(sqrt(galaxy.x_max*galaxy.x_max + galaxy.y_max*galaxy.y_max), rBound[rBound.size()-1]));
      
      const double zMax = std::min(galaxy.z_max, zBound[zBound.size()-1]);
      const double zMin = std::max(galaxy.z_min, -zBound[zBound.size()-1]);
      
      const int nRBins = int((rMax)/dR);
      
      const int nZBins = int((zMax - zMin)/dZ);
      
      assert (nRBins >= 1);
      assert (nZBins >= 1);
      
      //cout << "Bounds: " << rMax << " " << zMax << " " << zMin << endl;
      
      const vec3 normalUpper(0, 0, 1), normalSide(1, 0, 0);
      
      for (int iR = 0; iR < nRBins; ++iR) {
	
	const double dA = utl::kPi*(pow((iR+1)*dR, 2.) - pow(iR*dR, 2.))*pow(kpc2cm, 2.); // kpc2 -> cm2
	
	const double rr = (iR + 0.5)*dR, zz = zBound[zBound.size()-1] - 0.5*dZ;
	
	const rf::RadiationField::ThreeVector pos(rr, 0, zz);
	
	Skymap<double> skymap = rf.GetSkymap(pos, rf::RadiationField::TOTAL, galdef.ISRF_healpixOrder);
	
	valarray<double> flux(0., galaxy.nu_ISRF.size());
	
	vector<int> pixelList;
	
	const pointing lb(utl::kPi, 0.), ub(utl::kPiOnTwo, utl::kTwoPi);
	
	skymap.query_rectangle(lb, ub, pixelList);
	
	for (int iPix = 0; iPix < pixelList.size(); ++iPix) {
	  
	  if (pixelList[iPix] >= 0) {
	    
	    const ArraySlice<double>& arr = skymap[pixelList[iPix]];
	    
	    const vec3 direction = skymap.pix2vec(pixelList[iPix]);
	    
	    const double cosTheta = dotprod(normalUpper, direction);
	    
	    //cout << "Upper: " << iPix << " " << pixelList[iPix] << " " << cosTheta << " " << direction << " " << normalUpper << endl;
	    
	    for (int iE = 0; iE < flux.size(); ++iE)
	      flux[iE] += arr[iE];//*cosTheta;
	    
	  }
	  
	}
	
	galaxy.fISRFLuminosity += flux*targetE*targetE*dA*skymap.solidAngle()*eV_to_erg*utl::kSpeedOfLight_SI*utl::m/utl::cm;
	
	/*cout << "BoundR: " << iR << " " << rr << " " << zz << endl;
	  
	  for (int iE = 0; iE < flux.size(); ++iE)
	  cout << iE << " " << flux[iE]*targetE[iE]*targetE[iE]*skymap.solidAngle() << endl;
	*/
      }
      
      for (int iZ = 0; iZ < nZBins; ++iZ) {
	
	const double dA = 2.*utl::kPi*rBound[rBound.size()-1]*dZ*pow(kpc2cm, 2.); // kpc2 -> cm2
	
	const double rr = rBound[rBound.size()-1] - 0.5*dR, zz = -zBound[zBound.size()-1] + (iZ + 0.5)*dZ;
	
	const rf::RadiationField::ThreeVector pos(rr, 0, zz);
	
	Skymap<double> skymap = rf.GetSkymap(pos, rf::RadiationField::TOTAL, galdef.ISRF_healpixOrder);
	
	valarray<double> flux(0., galaxy.nu_ISRF.size());

	vector<int> pixelList;
	
	const pointing lb(utl::kPi, -utl::kPiOnTwo), ub(0., utl::kPiOnTwo);
	
	skymap.query_rectangle(lb, ub, pixelList);
	
	for (int iPix = 0; iPix < pixelList.size(); ++iPix) {
	  
	  if (pixelList[iPix] >= 0) {
	    
	    const ArraySlice<double>& arr = skymap[pixelList[iPix]];
	    
	    const vec3 direction = skymap.pix2vec(pixelList[iPix]);
	    
	    const double cosTheta = dotprod(normalSide, direction);
	    
	    //cout << "Lower: " << iPix << " " << pixelList[iPix] << " " << cosTheta << " " << direction << " " << normalSide << endl;
	    
	    for (int iE = 0; iE < flux.size(); ++iE)
	      flux[iE] += arr[iE];//*cosTheta;
	    
	  }
	  
	}
	
	galaxy.fISRFLuminosity += flux*targetE*targetE*dA*skymap.solidAngle()*eV_to_erg*utl::kSpeedOfLight_SI*utl::m/utl::cm;
	
	/*cout << "BoundZ: " << iZ << " " << rr << " " << zz << endl;
	  
	  for (int iE = 0; iE < flux.size(); ++iE)
	  cout << iE << " " << flux[iE]*targetE[iE]*targetE[iE]*skymap.solidAngle() << endl;
	*/
      }

    } else if (3 == galaxy.n_spatial_dimensions) {

    } else {

      INFO("Not 2 or 3 dimensional. Exiting.");
      exit(-1);

    }
    
    //exit(0);

    //for (int iE = 0; iE < galaxy.fISRFLuminosity.size(); ++iE)
    //cout << iE << " " << targetE[iE] << " " << galaxy.fISRFLuminosity[iE] << endl;
    
    //exit(0);

  }
  
  double nHI_total, nH2_total, nHII_total, nH_total, volume_total;
  
  nHI_total = nH2_total = nHII_total = nH_total = volume_total = 0;

  valarray<double> nHI_column_density, nH2_column_density, nHII_column_density;

  valarray< pair<double, double> > coordinates;

  if (2 == galaxy.n_spatial_dimensions) {

    // Cosmic rays first

    for (int i = 0; i < n_species; ++i) {

      Particle particle;

      particle = gcr[i];
      particle.create_transport_arrays();
      fill_transport_arrays(particle);

      valarray<double> crLuminosity(0., galaxy.n_pgrid);

      if ("secondary_positrons" != gcr[i].name &&
	  "secondary_electrons" != gcr[i].name &&
	  "secondary_antiprotons" != gcr[i].name) {

	for (int iR = 0; iR < galaxy.n_rgrid; ++iR)
	  for (int iZ = 0; iZ < galaxy.n_zgrid; ++iZ)
	    for (int iP = 0; iP < galaxy.n_pgrid; ++iP)
	      crLuminosity[iP] += 
		particle.primary_source_function.d2[iR][iZ].s[iP]/
		particle.beta[iP]*pow(particle.Ekin[iP], 2.)*
		2.*utl::kPi*galaxy.r[iR]*galaxy.dr*galaxy.dz;

      } else {

	gen_secondary_source(particle);

	for (int iR = 0; iR < galaxy.n_rgrid; ++iR)
	  for (int iZ = 0; iZ < galaxy.n_zgrid; ++iZ)
	    for (int iP = 0; iP < galaxy.n_pgrid; ++iP)
	      crLuminosity[iP] += 
		particle.secondary_source_function.d2[iR][iZ].s[iP]/
		particle.beta[iP]*pow(particle.Ekin[iP], 2.)*
		2.*utl::kPi*galaxy.r[iR]*galaxy.dr*galaxy.dz;
	
      }

      crLuminosity *= 4.*utl::kPi/c*1.e6*eV_to_erg*pow(kpc2cm, 3.);

      if ("Hydrogen_1" == gcr[i].name)
	galaxy.fProtonLuminosity += crLuminosity;
      else if ("Helium_4" == gcr[i].name)
	galaxy.fHeliumLuminosity += crLuminosity;
      else if ("primary_electrons" == gcr[i].name)
	galaxy.fPrimaryElectronLuminosity += crLuminosity;
      else if ("secondary_positrons" == gcr[i].name ||
	       "secondary_electrons" == gcr[i].name) 
	galaxy.fSecondaryElectronLuminosity += crLuminosity;
      else
	galaxy.fNucleiLuminosity += crLuminosity;

    }

    nHI_column_density.resize(galaxy.n_rgrid);
    nH2_column_density.resize(galaxy.n_rgrid);
    nHII_column_density.resize(galaxy.n_rgrid);
    coordinates.resize(galaxy.n_rgrid);

    nHI_column_density = nH2_column_density = nHII_column_density = 0;

    const int irmax = galaxy.n_rgrid - 1; // NB last point
    //  irmax=6; // for testing

    for (int iR = 0; iR < irmax; ++iR) 
      coordinates[iR] = pair<double, double>(galaxy.r[iR], 0.);

    ostringstream maxRbuf;

    maxRbuf << "Maximum radius for luminosity = " << galaxy.r[irmax];

    INFO(maxRbuf.str());

    for (int ir = 0; ir < irmax; ++ir) {

      for (int iz = 0; iz < galaxy.n_zgrid; ++iz) {
	
	const double volume = utl::kPi*(pow(galaxy.r[ir+1], 2.) - pow(galaxy.r[ir], 2.))*galaxy.dz*pow(kpc2cm, 3.) ; // not exact of course because emissivity is at r[ir]. volume in kpc^3 -> cm^3

	const double nHI_ = nHI3D(galaxy.r[ir], 0, galaxy.z[iz]);
	const double nH2_ = nH23D(galaxy.r[ir], 0, galaxy.z[iz], fX_CO(galaxy.r[ir]))*2.0; // Gulli20100218 Using fX_CO for X ratio. NB two atoms per molecule!
	const double nHII_ = nHII3D(galaxy.r[ir], 0, galaxy.z[iz]);
	const double nH = nHI_ + nH2_;// + nHII_;
     
	nHI_total += nHI_*volume; 
	nH2_total += nH2_*volume;
	nHII_total += nHII_*volume;
	nH_total += nH*volume;

	nHI_column_density[ir] += nHI_*galaxy.dz*kpc2cm;
	nH2_column_density[ir] += nH2_*galaxy.dz*kpc2cm;
	nHII_column_density[ir] += nHII_*galaxy.dz*kpc2cm;

	volume_total += volume;	

	const double volumeFactor = volume*4.*utl::kPi, gasFactor = volumeFactor*(nH + nHII_);

	for (int ip = 0; ip < galaxy.n_E_gammagrid; ++ip) {
	  
	  // Emissivity is per sr so need 4pi factor
	  
	  const double neutralBrem = std::isnan(galaxy.bremss_emiss.d2[ir][iz].s[ip]) ? 0. : galaxy.bremss_emiss.d2[ir][iz].s[ip];

	  const double ionisedBrem = std::isnan(galaxy.bremss_ionized_emiss.d2[ir][iz].s[ip]) ? 0. : galaxy.bremss_ionized_emiss.d2[ir][iz].s[ip];
	  
	  galaxy.fBremsstrahlungLuminosity[ip] += (neutralBrem*nH + ionisedBrem*nHII_)*volumeFactor*galaxy.E_gamma[ip]*galaxy.E_gamma[ip];

	  galaxy.fPi0DecayLuminosity[ip] += galaxy.pi0_decay_emiss.d2[ir][iz].s[ip]*gasFactor*galaxy.E_gamma[ip]*galaxy.E_gamma[ip];
	  //galaxy.fDMLuminosity[ip] += galaxy.DM_emiss.d2[ir][iz].s[ip]*volumeFactor; // Fix this shortly ... TAP20111101

	  for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp)
	    galaxy.fICLuminosity[ip] += galaxy.IC_iso_emiss[icomp].d2[ir][iz].s[ip]*volumeFactor*galaxy.E_gamma[ip]*galaxy.E_gamma[ip];

	}

	for (int inu = 0; inu < galaxy.n_nu_synchgrid; ++inu) {
	  
	  galaxy.fSynchrotronLuminosity[inu] += galaxy.synchrotron_emiss.d2[ir][iz].s[inu]*volumeFactor;
	
	}
	
      }

    }

  }

  if (3 == galaxy.n_spatial_dimensions) {

    // Cosmic rays first
    
    for (int i = 0; i < n_species; ++i) {

      Particle particle;

      particle = gcr[i];
      particle.create_transport_arrays();
      fill_transport_arrays(particle);

      valarray<double> crLuminosity(0., galaxy.n_pgrid);

      if ("secondary_positrons" != gcr[i].name &&
	  "secondary_electrons" != gcr[i].name &&
	  "secondary_antiprotons" != gcr[i].name) {

	for (int iX = 0; iX < galaxy.n_xgrid; ++iX)
	  for (int iY = 0; iY < galaxy.n_ygrid; ++iY)
	    for (int iZ = 0; iZ < galaxy.n_zgrid; ++iZ)
	      for (int iP = 0; iP < galaxy.n_pgrid; ++iP)
		crLuminosity[iP] += 
		  particle.primary_source_function.d3[iX][iY][iZ].s[iP]/
		  particle.beta[iP]*pow(particle.Ekin[iP], 2.)*
		  galaxy.dx*galaxy.dy*galaxy.dz;
	
      } else {

	gen_secondary_source(particle);

	for (int iX = 0; iX < galaxy.n_xgrid; ++iX)
	  for (int iY = 0; iY < galaxy.n_ygrid; ++iY)
	    for (int iZ = 0; iZ < galaxy.n_zgrid; ++iZ)
	      for (int iP = 0; iP < galaxy.n_pgrid; ++iP)
		crLuminosity[iP] += 
		  particle.secondary_source_function.d3[iX][iY][iZ].s[iP]/
		  particle.beta[iP]*pow(particle.Ekin[iP], 2.)*
		  galaxy.dx*galaxy.dy*galaxy.dz;
	
      }

      crLuminosity *= 4.*utl::kPi/c*1.e6*eV_to_erg*pow(kpc2cm, 3.);

      if ("Hydrogen_1" == gcr[i].name)
	galaxy.fProtonLuminosity += crLuminosity;
      else if ("Helium_4" == gcr[i].name)
	galaxy.fHeliumLuminosity += crLuminosity;
      else if ("primary_electrons" == gcr[i].name)
	galaxy.fPrimaryElectronLuminosity += crLuminosity;
      else if ("secondary_positrons" == gcr[i].name ||
	       "secondary_electrons" == gcr[i].name) 
	galaxy.fSecondaryElectronLuminosity += crLuminosity;
      else
	galaxy.fNucleiLuminosity += crLuminosity;

    }

    //nHI_column_density.resize(galaxy.n_rgrid);
    //nH2_column_density.resize(galaxy.n_rgrid);
    //nHII_column_density.resize(galaxy.n_rgrid);
    //coordinates.resize(galaxy.n_rgrid);

    //nHI_column_density = nH2_column_density = nHII_column_density = 0;

    //const int irmax = galaxy.n_rgrid - 1; // NB last point
    //  irmax=6; // for testing

    //for (int iR = 0; iR < irmax; ++iR) 
    //coordinates[iR] = pair<double, double>(galaxy.r[iR], 0.);

    //ostringstream maxRbuf;

    //maxRbuf << "Maximum radius for luminosity = " << galaxy.r[irmax];

    //INFO(maxRbuf.str());
    
    for (int iX = 0; iX < galaxy.n_xgrid; ++iX) 
      for (int iY = 0; iY < galaxy.n_ygrid; ++iY) 
	for (int iZ = 0; iZ < galaxy.n_zgrid; ++iZ) {
	  
	  const double r = sqrt(galaxy.x[iX]*galaxy.x[iX] + galaxy.y[iY]*galaxy.y[iY]);
	
	  const double volume = galaxy.dx*galaxy.dy*galaxy.dz*pow(kpc2cm, 3.); // volume in kpc^3 -> cm^3
	  
	  const double nHI_ = nHI3D(galaxy.x[iX], galaxy.y[iY], galaxy.z[iZ]);
	  const double nH2_ = nH23D(galaxy.x[iX], galaxy.y[iY], galaxy.z[iZ], fX_CO(r))*2.0; // Gulli20100218 Using fX_CO for X ratio. NB two atoms per molecule!
	  const double nHII_ = nHII3D(galaxy.x[iX], galaxy.y[iY], galaxy.z[iZ]);
	  const double nH = nHI_ + nH2_;// + nHII_;
     
	  nHI_total += nHI_*volume; 
	  nH2_total += nH2_*volume;
	  nHII_total += nHII_*volume;
	  nH_total += nH*volume;
	  
	  //nHI_column_density[ir] += nHI_*galaxy.dz*kpc2cm;
	  //nH2_column_density[ir] += nH2_*galaxy.dz*kpc2cm;
	  //nHII_column_density[ir] += nHII_*galaxy.dz*kpc2cm;
	  
	  volume_total += volume;	
	  
	  const double volumeFactor = volume*4.*utl::kPi, gasFactor = volumeFactor*(nH + nHII_);

	  for (int ip = 0; ip < galaxy.n_E_gammagrid; ++ip) {
	  
	    // Emissivity is per sr so need 4pi factor above
	    
	    const double neutralBrem = std::isnan(galaxy.bremss_emiss.d3[iX][iY][iZ].s[ip]) ? 0. : galaxy.bremss_emiss.d3[iX][iY][iZ].s[ip];

	    const double ionisedBrem = std::isnan(galaxy.bremss_ionized_emiss.d3[iX][iY][iZ].s[ip]) ? 0. : galaxy.bremss_ionized_emiss.d3[iX][iY][iZ].s[ip];

	    galaxy.fBremsstrahlungLuminosity[ip] += (neutralBrem*nH + ionisedBrem*nHII_)*volumeFactor*galaxy.E_gamma[ip]*galaxy.E_gamma[ip];
	  
	    galaxy.fPi0DecayLuminosity[ip] += galaxy.pi0_decay_emiss.d3[iX][iY][iZ].s[ip]*gasFactor*galaxy.E_gamma[ip]*galaxy.E_gamma[ip];
	    //galaxy.fDMLuminosity[ip] += galaxy.DM_emiss.d2[ir][iz].s[ip]*volumeFactor; // Fix this shortly ... TAP20111101
	    
	    for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp)
	      galaxy.fICLuminosity[ip] += galaxy.IC_iso_emiss[icomp].d3[iX][iY][iZ].s[ip]*volumeFactor*galaxy.E_gamma[ip]*galaxy.E_gamma[ip];
	    
	  }
	  
	  for (int inu = 0; inu < galaxy.n_nu_synchgrid; ++inu) {
	    
	    galaxy.fSynchrotronLuminosity[inu] += galaxy.synchrotron_emiss.d3[iX][iY][iZ].s[inu]*volumeFactor;
	
	  }
	
	}

  }

  galaxy.fSynchrotronLuminosity *= galaxy.nu_synch; // erg Hz^-1 s^-1 -> erg s^-1
  galaxy.fICLuminosity *= MEV2ERG; // MeV s^-1 -> erg s^-1
  galaxy.fPi0DecayLuminosity *= MEV2ERG;
  galaxy.fBremsstrahlungLuminosity *= MEV2ERG;
  galaxy.fDMLuminosity *= MEV2ERG;
 
  valarray<double> gamma_total_luminosity_spectrum(0., galaxy.n_E_gammagrid);
  
  gamma_total_luminosity_spectrum = galaxy.fBremsstrahlungLuminosity + galaxy.fPi0DecayLuminosity + galaxy.fICLuminosity + galaxy.fDMLuminosity;
  
  // Because emissivity is MeV^2 cm-3 sr-1 MeV-1 s-1 
  // the units are  MeV s-1  as in Strong etal 2000 ApJ 537, 763 Fig 20.
  
  /*valarray<double> gamma_total_integral_luminosity(0., galaxy.n_E_gammagrid);
  valarray<double> gamma_bremss_integral_luminosity(0., galaxy.n_E_gammagrid);
  valarray<double> gamma_pi0_decay_integral_luminosity(0., galaxy.n_E_gammagrid);
  valarray<double> gamma_IC_integral_luminosity(0., galaxy.n_E_gammagrid);
  valarray<double> gamma_DM_integral_luminosity(0., galaxy.n_E_gammagrid);
  
  valarray<double> gamma_total_integral_photons(0., galaxy.n_E_gammagrid);
  valarray<double> gamma_bremss_integral_photons(0., galaxy.n_E_gammagrid);
  valarray<double> gamma_pi0_decay_integral_photons(0., galaxy.n_E_gammagrid);
  valarray<double> gamma_IC_integral_photons(0., galaxy.n_E_gammagrid);
  valarray<double> gamma_DM_integral_photons(0., galaxy.n_E_gammagrid);
  
  for (int ip = 0; ip < galaxy.n_E_gammagrid; ++ip)
    for (int ipp = ip; ipp < galaxy.n_E_gammagrid; ++ipp) {
      
      gamma_bremss_integral_luminosity[ip] += galaxy.fBremsstrahlungLuminosity[ipp];
      
      gamma_pi0_decay_integral_luminosity[ip] += galaxy.fPi0DecayLuminosity[ipp];
      
      gamma_IC_integral_luminosity[ip] += galaxy.fICLuminosity[ipp];
      
      gamma_DM_integral_luminosity[ip] += galaxy.fDMLuminosity[ipp];
      
      gamma_total_integral_luminosity[ip] = 
	gamma_bremss_integral_luminosity[ip] + 
	gamma_pi0_decay_integral_luminosity[ip] + 
	gamma_IC_integral_luminosity[ip] + 
	gamma_DM_integral_luminosity[ip];
      
      gamma_bremss_integral_photons[ip] += galaxy.fBremsstrahlungLuminosity[ipp]/galaxy.E_gamma[ipp];
      
      gamma_pi0_decay_integral_photons[ip] += galaxy.fPi0DecayLuminosity[ipp]/galaxy.E_gamma[ipp];
      
      gamma_IC_integral_photons[ip] += galaxy.fICLuminosity[ipp]/galaxy.E_gamma[ipp];
      
      gamma_DM_integral_photons[ip] += galaxy.fDMLuminosity[ipp]/galaxy.E_gamma[ipp];
      
      gamma_total_integral_photons[ip] += gamma_bremss_integral_photons[ip] + gamma_pi0_decay_integral_photons[ip] + gamma_IC_integral_photons[ip] + gamma_DM_integral_photons[ip]; // *E instead of *E^2
      
    }
  
  // integral Ef(E)dE = integral E^2f(E) dlog(E) = sum E^2 f(E) delta logE.  delta logE = log(E factor)
  gamma_total_integral_luminosity *= log(galaxy.E_gamma_factor);
  gamma_total_integral_luminosity *= MEV2ERG;// from constants.h
  
  gamma_bremss_integral_luminosity *= log(galaxy.E_gamma_factor);
  gamma_bremss_integral_luminosity *= MEV2ERG;// from constants.h
  
  gamma_pi0_decay_integral_luminosity *= log(galaxy.E_gamma_factor);
  gamma_pi0_decay_integral_luminosity *= MEV2ERG;// from constants.h
  
  gamma_IC_integral_luminosity *= log(galaxy.E_gamma_factor);
  gamma_IC_integral_luminosity *= MEV2ERG;// from constants.h
  
  // integral  f(E)dE = integral Ef(E) dlog(E) = sum E f(E) delta logE.  delta logE = log(E factor)
  gamma_total_integral_photons *= log(galaxy.E_gamma_factor);
  gamma_bremss_integral_photons *= log(galaxy.E_gamma_factor);
  gamma_pi0_decay_integral_photons *= log(galaxy.E_gamma_factor);
  gamma_IC_integral_photons *= log(galaxy.E_gamma_factor);
  
  valarray<double> synchrotron_integral_luminosity(0., galaxy.n_nu_synchgrid);
  
  for (int inu = 0; inu < galaxy.n_nu_synchgrid; ++inu)
    for (int iinu = inu; iinu < galaxy.n_nu_synchgrid; ++iinu)
      synchrotron_integral_luminosity[inu] += galaxy.fSynchrotronLuminosity[iinu]*galaxy.nu_synch[iinu];
  
  synchrotron_integral_luminosity *= log(galaxy.nu_synch_factor);
  */  
  // ======================== print out the results =====================================
  
  // integral properties of Galaxy
  
  ostringstream gasBuf1, gasBuf2, volBuf;
  
  volBuf << "Total galactic volume: " << volume_total << " cm^3";
  
  INFO(volBuf.str());
  
  gasBuf1 << "Total gas atoms: " 
	  << nHI_total << " (HI) " 
	  << nH2_total << " (H2) " 
	  << nHII_total << " (HII) " 
	  << nH_total << " (total)";
  gasBuf2 << "Total gas mass (Msun): " 
	  << nHI_total*utl::kHydrogenMass_SI/utl::kSolarMass_SI << " (HI) " 
	  << nH2_total*utl::kHydrogenMass_SI/utl::kSolarMass_SI << " (H2) " 
	  << nHII_total*utl::kHydrogenMass_SI/utl::kSolarMass_SI << " (HII) " 
	  << nH_total*utl::kHydrogenMass_SI/utl::kSolarMass_SI << " (total) ";
  
  INFO(gasBuf1.str());
  INFO(gasBuf2.str());
  
  const double sigmaFactor = utl::kHydrogenMass_SI/utl::kSolarMass_SI*pow(kpc2cm, 2.)*1.0e-6; // atoms cm-2 to Msun/pc^2
  
  for (int iGas = 0; iGas < nHI_column_density.size(); ++iGas) {
    
    ostringstream numBuf, massBuf;
    /*
      numBuf << "z-column density for " << (2 == galaxy.n_spatial_dimensions ? " R = " << galaxy.r[ir] << " kpc, " 
      << nHI_column_density[ir] << " (HI atoms cm^-2) "
      << nH2_column_density[ir] << " (H2 atmos cm^-2) "
      << nHII_column_density[ir] << " (HII atoms cm^-2)";
      
      massBuf << "z-column density: R = " << galaxy.r[ir] << " kpc, " 
      << nHI_column_density[ir]*sigmaFactor << " (HI Msun pc^-2) "
      << nH2_column_density[ir]*sigmaFactor << " (H2 Msun pc^-2) "
      << nHII_column_density[ir]*sigmaFactor << " (HII Msun pc^-2)";
      
      INFO(numBuf.str());
      INFO(massBuf.str());
      
      }
    */
    
  }   

  /*ostringstream crBuf;
  crBuf << "Cosmic-ray luminosity (erg s^-1): " << endl
	<< "Energy "
	<< "Protons "
	<< "Helium "
	<< "Nuclei "
	<< "Primary electrons " 
	<< "Secondary electrons" << endl;
 
  for (int iP = 0; iP < galaxy.n_pgrid; ++iP)
    crBuf  << galaxy.p[iP] << " "
	   << galaxy.fProtonLuminosity[iP] << " "
	   << galaxy.fHeliumLuminosity[iP] << " "
	   << galaxy.fNucleiLuminosity[iP] << " "
	   << galaxy.fPrimaryElectronLuminosity[iP] << " "
	   << galaxy.fSecondaryElectronLuminosity[iP] << endl;
   
  INFO(crBuf.str());

  ostringstream gammaBuf;
  gammaBuf << "Gamma-ray luminosity (erg s^-1) and photons (ph s^-1): " << endl 
	   << "Energy "
	   << "Bremss "
	   << "pi0 "
	   << "IC "
	   << "DM "
	   << "Total "
	   << "Integral bremss "
	   << "Integral pi0 "
	   << "Integral IC "
	   << "Integral DM "
	   << "Integral Total "
	   << "Photons bremss "
	   << "Photons pi0 "
	   << "Photons IC "
	   << "Photons DM "
	   << "Photons Total" << endl;
  
  for (int ip = 0; ip < galaxy.n_E_gammagrid; ++ip)
    gammaBuf << galaxy.E_gamma[ip] << " "
	     << galaxy.fBremsstrahlungLuminosity[ip] << " "
	     << galaxy.fPi0DecayLuminosity[ip] << " "
	     << galaxy.fICLuminosity[ip] << " "
	     << galaxy.fDMLuminosity[ip] << " "
	     << gamma_total_luminosity_spectrum[ip] << " "
	     << gamma_bremss_integral_luminosity[ip] << " "
	     << gamma_pi0_decay_integral_luminosity[ip] << " "
	     << gamma_IC_integral_luminosity[ip] << " "
	     << gamma_DM_integral_luminosity[ip] << " "
	     << gamma_total_integral_luminosity[ip] << " "
	     << gamma_bremss_integral_photons[ip] << " "
	     << gamma_pi0_decay_integral_photons[ip] << " "
	     << gamma_IC_integral_photons[ip] << " "
	     << gamma_DM_integral_photons[ip] << " "
	     << gamma_total_integral_photons[ip] << endl;
    
  //gammaBuf << "Energy = " << galaxy.E_gamma[ip] << " "  
  //	 << "bremss: " << gamma_bremss_luminosity_spectrum[ip] << " "
  //	 << "pi0: " << gamma_pi0_decay_luminosity_spectrum[ip] << " "
  //	 << "IC: " << gamma_IC_luminosity_spectrum[ip] << " "
  //	 << "total: " << gamma_total_luminosity_spectrum[ip] << " "
  //	 << "integral luminosity erg s-1 :"
  //	 << " bremss: "<< gamma_bremss_integral_luminosity   [ip]<<" "
  //	 <<  " pi0: "  <<gamma_pi0_decay_integral_luminosity [ip]<<" "
  //	 <<  " IC: "   <<gamma_IC_integral_luminosity        [ip]<<" "
  //	 <<  " total: "<<gamma_total_integral_luminosity     [ip]<<" "
  //	 << " integral photons s-1:"            
  //	 << " bremss: "<< gamma_bremss_integral_photons      [ip]<<" "
  //	 <<  " pi0: "  <<gamma_pi0_decay_integral_photons    [ip]<<" "
  //	 <<  " IC: "   <<gamma_IC_integral_photons           [ip]<<" "
  //	 <<  " total: "<<gamma_total_integral_photons        [ip]<<" "
  //<<endl;
  
  INFO(gammaBuf.str());
  
  ostringstream syncBuf;
  
  syncBuf << "Synchrotron luminosity (erg Hz^-1 s^-1): " << endl;
  syncBuf << "Frequency " 
	  << "Differential " 
	  << "Integral" << endl;
  
  for (int iNu = 0; iNu < galaxy.n_nu_synchgrid; ++iNu)
    syncBuf << galaxy.nu_synch[iNu] << " "
	    << galaxy.fSynchrotronLuminosity[iNu] << " " 
	    << synchrotron_integral_luminosity[iNu] << endl;
  
  INFO(syncBuf.str());
  */  
  //for (int inu = 0; inu < galaxy.n_nu_synchgrid; ++inu) {
  //cout<<"nu ="<<galaxy.nu_synch[inu]<<" synchrotron luminosity erg Hz-1 s-1: "<< synchrotron_luminosity_spectrum[inu]
  // <<" integral luminosity erg s-1: "<< synchrotron_integral_luminosity[inu]
  //<<endl;
  //}
  
  INFO("Exit");
  
  return 0;
 
}
