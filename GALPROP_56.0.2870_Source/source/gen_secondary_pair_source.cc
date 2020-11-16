#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <valarray>

using namespace std;

#include "galprop_classes.h"
#include "galprop_internal.h"

#include <ErrorLogger.h>
#include <PhysicalConstants.h>

// See, e.g., Astropart. Phys. 12 (2000) 217-238, or Jauch and Rohrlich

double CrossSectionPairProduction(const double x) {

  double result = 0;

  if (x < 1.) {

    const double beta = sqrt(1. - x);

    const double ln = (x < 1.e-4 ? log(2.*(1. + beta)/x) : log((1. + beta)/(1. - beta)));

    const double f = 0.5*x*((2. + x*(2. - x))*ln - 2.*beta*(1. + x));

    result = 3./8.*(utl::kThomsonCrossSection_SI*utl::m/utl::cm*utl::m/utl::cm)*f; // cm^2

  }

  return result;

}

double InteractionLengthPairProduction(const double gammaE, 
				       const valarray<double>& targetE,
				       const valarray<double>& numberDensity) {

  double result = 0;

  const int cosThetaBins = 101, bins = targetE.size();

  const double dCosTheta = 2./(cosThetaBins - 1);

  const double fac = 0.5*log(targetE[1]/targetE[0]);

  for (int i = 0; i < bins-1; ++i) {

    const double fac1 = targetE[i]*numberDensity[i];
    const double fac2 = targetE[i+1]*numberDensity[i+1];

    const double eps1 = targetE[i]/utl::kElectronMass, eps2 = targetE[i+1]/utl::kElectronMass, eps = gammaE/utl::kElectronMass;

    for (int j = 0; j < cosThetaBins-1; ++j) {

      const double cosTheta = -1. + (j+0.5)*dCosTheta;

      const double Ec1 = 0.5*eps1*eps*(1. - cosTheta), x1 = 1./Ec1;
      const double Ec2 = 0.5*eps2*eps*(1. - cosTheta), x2 = 1./Ec2;
			      
      result += fac*dCosTheta*(1. - cosTheta)*(fac1*CrossSectionPairProduction(x1) + fac2*CrossSectionPairProduction(x2));
     
    }

  }

  return result; // cm^-1

}

// Returns particles cm^2 -- this means for one species (e- or e+ only)
// Uses Eq. 32 from A&A 325, 866 (1997) -- eps1 and eps2 have been 
// swapped: here eps2 is the dimensionless gamma ray energy and eps1 is 
// the dimensionless target photon energy. In the Boettecher & Schlickeiser 
// paper they have the opposite meaning

double CrossSectionIsotropicBoettecher(const double eps2, 
				       const double eps1, 
				       const double gamma) {
  
  double result = 0;
  
  const double epsMin = eps2/4./gamma/(eps2 - gamma);

  if (eps1 >= epsMin) {

    result = 1./8.*
      0.75*(utl::kThomsonCrossSection_SI*utl::m/utl::cm*utl::m/utl::cm)*
      (4.*eps2*eps2/gamma/(eps2 - gamma)*log(4.*eps1*gamma*(eps2 - gamma)/eps2) -
       8.*eps2*eps1 +
       2.*eps2*eps2/gamma/(eps2 - gamma)*(2.*eps2*eps1 - 1.) -
       (1. - 1./(eps2*eps1))*pow(eps2*eps2/gamma/(eps2 - gamma), 2.));
    
  }

  return result;

}

void Galprop::gen_secondary_pair_source(Particle& particle) {

  INFO("Entry");

  // Use isotropic approximation for gamma rays and target photons. 
  // To get this working, calculate the gamma-ray number density using simply 
  // emissivity/r^2 integrated over the Galaxy for each point. Next step 
  // will be to generate a coarse gamma intensity skymap at each point, 
  // including the absorption on the ISRF (simply by exp(-tau) weighting).

  Distribution pairProductionSource;

  const int targetBins = galaxy.nu_ISRF.size();

  valarray<double> gamma(0., galaxy.n_pgrid), eps2(0., galaxy.n_E_gammagrid), eps1(0., targetBins), gammaE(0., galaxy.n_E_gammagrid), targetE(0., targetBins);

  gamma = galaxy.p*(utl::MeV/utl::eV)*1./utl::kElectronMass;
  eps2 = galaxy.E_gamma*(utl::MeV/utl::eV)*1./utl::kElectronMass;
  eps1 = galaxy.nu_ISRF*(utl::kPlanck_SI/utl::e_SI)*1./utl::kElectronMass;

  gammaE = eps2*utl::kElectronMass;
  targetE = eps1*utl::kElectronMass;
  
  const double eps2Delta = log(eps2[1]/eps2[0]);
  const double eps1Delta = log(eps1[1]/eps1[0]);

  //valarray<double> isotropicIntensity(0., galaxy.n_E_gammagrid);

  //const double isotropicNormalisation = 1.03e-5; // Integral intensity > 100 MeV cm^-2 s^-1 sr^-1
  //const double isotropicIndex = 2.41; // Differential index
  //const double isotropicCutoffEnergy = 10.*utl::TeV;

  //for (int i = 0; i < galaxy.n_E_gammagrid; ++i)
  //isotropicIntensity[i] = isotropicNormalisation*pow(gammaE[i], -isotropicIndex)*exp(-gammaE[i]/isotropicCutoffEnergy);

  //isotropicIntensity *= -(1. - isotropicIndex)/pow(100.*utl::MeV, 1. - isotropicIndex)*utl::MeV/utl::eV;

  //for (int i = 0; i < galaxy.n_E_gammagrid; ++i)
  //cout << i << " " << gammaE[i] << " " << isotropicIntensity[i] << endl;

  //exit(0);
  
  assert(eps2Delta > 0);
  assert(eps1Delta > 0);
  
  const double ds = 0.5; // kpc -- hard wired for the moment!

  if (2 == galdef.n_spatial_dimensions) {

    const double dTheta = 10.*utl::kConvertDegreesToRadians; // Hard wired for now ...

    const int rBins = galaxy.n_rgrid, zBins = galaxy.n_zgrid, thetaBins = utl::kPi/dTheta;

    assert (thetaBins > 0);

    INFO("Calculating inverse absorption path-length");

#pragma omp parallel for schedule(dynamic) default(shared)
    for (int iR = 0; iR < rBins; ++iR)
      for (int iZ = 0; iZ < zBins; ++iZ) {

	const valarray<double>& numberDensity = galaxy.fTotalISRFNumberDensity.d2[iR][iZ].GetSpec();

	for (int iG = 0; iG < galaxy.n_E_gammagrid; ++iG)
	  galaxy.fPairAbsorptionInvLength.d2[iR][iZ].s[iG] = 
	    InteractionLengthPairProduction(gammaE[iG], targetE, numberDensity);

      }

    // Accumulate total emissivity throughout calculation volume

    Distribution gammaRayEmissivity(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);

    INFO("Calculating gamma-ray emissivity");

#pragma omp parallel for schedule(dynamic) default(shared)
    for (int iR = 0; iR < rBins; ++iR)
      for (int iZ = 0; iZ < zBins; ++iZ) {

	valarray<double> ic(0., galaxy.n_E_gammagrid);
	
	for (int i = 0; i < galaxy.n_ISRF_components; ++i) 
	  ic += galaxy.IC_iso_emiss[i].d2[iR][iZ].GetSpec();

	for (int iE = 0; iE < galaxy.n_E_gammagrid; ++iE) {
 
	  gammaRayEmissivity.d2[iR][iZ].s[iE] = ic[iE] + 
	    galaxy.pi0_decay_emiss.d2[iR][iZ].s[iE]*(galaxy.n_HI.d2[iR][iZ].s[0] + 2.*galaxy.n_H2.d2[iR][iZ].s[0] + galaxy.n_HII.d2[iR][iZ].s[0]) + 
	    galaxy.bremss_emiss.d2[iR][iZ].s[iE]*(galaxy.n_HI.d2[iR][iZ].s[0] + 2.*galaxy.n_H2.d2[iR][iZ].s[0]) +
	    galaxy.bremss_ionized_emiss.d2[iR][iZ].s[iE]*galaxy.n_HII.d2[iR][iZ].s[0]; // MeV^-1 cm^-3 s^-1 sr^-1

	}

      }

    Distribution gammaRayIntensity(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);

    INFO("Calculating gamma-ray intensity. This may take a while ...");

    //#pragma omp parallel for schedule(dynamic) default(shared)
    for (int iR = 0; iR < rBins; ++iR) {

      //#pragma omp critical 
      //{
      //cout << iR << endl;
	//}

#pragma omp parallel for schedule(dynamic) default(shared)
      for (int iZ = 0; iZ < zBins; ++iZ) {

	const double x = galaxy.r[iR], y = 0, z = galaxy.z[iZ];

	for (int iRR = 0; iRR < rBins; ++iRR)
	  for (int iZZ = 0; iZZ < zBins; ++iZZ)
	    for (int iTheta = 0; iTheta < thetaBins; ++iTheta) {
	      
	      const double theta = (iTheta + 0.5)*dTheta, xp = galaxy.r[iRR]*cos(theta), yp = galaxy.r[iRR]*sin(theta), zp = galaxy.z[iZZ];
	      
	      const double rp2 = (x - xp)*(x - xp) + (y - yp)*(y - yp) + (z - zp)*(z - zp), invRp2 = (rp2 > 0. ? 1./rp2 : 0.);

	      const double vol = (pow(galaxy.r[iRR] + galaxy.dr, 2.) - pow(galaxy.r[iRR], 2.))*dTheta*galaxy.dz;

	      // Calculate optical depth from (x, y, z) to (xp, yp, zp) as 
	      // function of gamma ray energy
		      
	      valarray<double> opticalDepth(0., galaxy.n_E_gammagrid);

	      // Get direction cosines, integrate along line of sight

	      const double rp = sqrt(rp2);

	      const int steps = (rp/ds > 1. ? int(rp/ds) : 1);

	      const double pX = (rp > 0. ? (xp - x)/rp : 0.);
	      const double pY = (rp > 0. ? (yp - y)/rp : 0.);
	      const double pZ = (rp > 0. ? (zp - z)/rp : 0.);

	      double xx = x, yy = y, zz = z;

	      for (int i = 0; i < steps; ++i) {

		xx += pX*ds;
		yy += pY*ds;
		zz += pZ*ds;

		const double rr = sqrt(xx*xx + yy*yy);

		if (zz < galaxy.z_min || zz > galaxy.z_max)
		  continue;

		if (rr > galaxy.r_max)
		  continue;

		const double dR = (rr - galaxy.r_min)/(galaxy.r_max - galaxy.r_min);
		const double dZ = (zz - galaxy.z_min)/(galaxy.z_max - galaxy.z_min);

		const unsigned int iNR = (dR < 1. ? dR*(galaxy.n_rgrid-1) : galaxy.n_rgrid-1);
		const unsigned int iNZ = (dZ < 1. ? dZ*(galaxy.n_zgrid-1) : galaxy.n_zgrid-1);
  
		//cout << iNR << " " << iNZ << " " << rBins << " " << zBins << endl;

		for (int iG = 0; iG < galaxy.n_E_gammagrid; ++iG)
		  opticalDepth[iG] += galaxy.fPairAbsorptionInvLength.d2[iNR][iNZ].s[iG]*utl::kpc/utl::cm*ds;

	      }
	      
	      for (int iG = 0; iG < galaxy.n_E_gammagrid; ++iG) {
		gammaRayIntensity.d2[iR][iZ].s[iG] +=  
		  ((std::isnan(gammaRayEmissivity.d2[iRR][iZZ].s[iG]) ? 0. : gammaRayEmissivity.d2[iRR][iZZ].s[iG])*invRp2*vol*
		   2.* // Comes from other half of galaxy
		   (utl::kKiloParsec*utl::m/utl::cm))* // MeV^-1 cm^-2 s^-1 sr^-1
		  exp(-opticalDepth[iG]); // Absorption

	      }
	      
	    }

      }
     	    
    }

    //for (int iE = 0; iE < galaxy.n_E_gammagrid; ++iE) 
    //cout << iE << " " << gammaE[iE] << " " << gammaRayIntensity.d2[0][0].s[iE] << " " << gammaRayIntensity.d2[0][zBins/2].s[iE] << " " << gammaRayIntensity.d2[17][0].s[iE] << " " << gammaRayIntensity.d2[17][zBins/2].s[iE] << endl; //" " << isotropicIntensity[iE] << endl;

    //exit(0);

    // Construct source function throughout Galactic volume

    pairProductionSource.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_pgrid);

    INFO("Calculating pair production source distribution");

#pragma omp parallel for schedule(dynamic) default(shared)
    for (int iR = 0; iR < rBins; ++iR) 
      for (int iZ = 0; iZ < zBins; ++iZ) {

	for (int iE = 0; iE < galaxy.n_pgrid; ++iE) {

	  // Using approximation from Boettcher and Schlickeiser

	  double pairBoettecher = 0;

	  const int iGMin = std::max(0, int(log(gamma[iE]/eps2[0])/eps2Delta));

	  for (int iG = iGMin; iG < galaxy.n_E_gammagrid; ++iG) {

	    double sum = 0;

	    if (eps2[iG] > gamma[iE]) {

	      const double epsMin = eps2[iG]/4./gamma[iE]/(eps2[iG] - gamma[iE]);

	      if (epsMin <= eps1[targetBins-1]) {
		
		const int iTMin = std::max(0, int(log(epsMin/eps1[0])/eps1Delta));
     
		for (int iT = iTMin; iT < targetBins; ++iT) {

		  const double fac = targetE[iT]*eps1Delta/eps1[iT]/eps1[iT]*CrossSectionIsotropicBoettecher(eps2[iG], eps1[iT], gamma[iE]); // eV cm^2
		
		  sum += (galaxy.ISRF[0].d2[iR][iZ].s[iT] + galaxy.ISRF[1].d2[iR][iZ].s[iT] + galaxy.ISRF[2].d2[iR][iZ].s[iT])*fac; // eV^-1 cm^-3 * eV cm^2 -> cm^-1

		}

	      }

	    }

	    const double fac = gammaRayIntensity.d2[iR][iZ].s[iG]*
	      galaxy.E_gamma[iG]*eps2Delta/eps2[iG]/eps2[iG]/eps2[iG]; // MeV^-1 cm^-2 s^-1 sr^-1 * MeV -> cm^-2 s^-1 sr^-1
	    
	    pairBoettecher += fac*sum; // cm^-1 * cm^-2 s^-1 sr^-1 -> cm^-3 s^-1 sr^-1

	  }

	  pairProductionSource.d2[iR][iZ].s[iE] = 
	    pairBoettecher* // cm^-3 s^-1 sr^-1
	    utl::kSpeedOfLight_SI*(utl::m/utl::cm)* // -> cm^-2 s^-2 sr^-1
	    (1./(utl::kElectronMass*(utl::eV/utl::MeV))); // -> MeV^-1 cm^-2 s^-2 sr^-1 -- same units as used in gen_secondary_positron_source.cc

	  assert(pairProductionSource.d2[iR][iZ].s[iE] >= 0.);
 
	}

      }

  }

  if (3 == galdef.n_spatial_dimensions) {

    const int xBins = galaxy.n_xgrid, yBins = galaxy.n_ygrid, zBins = galaxy.n_zgrid;

    INFO("Calculating inverse absorption path-length");

#pragma omp parallel for schedule(dynamic) default(shared)
    for (int iX = 0; iX < xBins; ++iX)
      for (int iY = 0; iY < yBins; ++iY)
	for (int iZ = 0; iZ < zBins; ++iZ) {

	  auto& numberDensity = galaxy.fTotalISRFNumberDensity.d3[iX][iY][iZ].GetSpec();

	  for (int iG = 0; iG < galaxy.n_E_gammagrid; ++iG)
	    galaxy.fPairAbsorptionInvLength.d3[iX][iY][iZ].s[iG] = 
	      InteractionLengthPairProduction(gammaE[iG], targetE, numberDensity);

	}


    // Accumulate total emissivity throughout calculation volume

    INFO("Calculating gamma-ray emissivity");
    
    Distribution gammaRayEmissivity(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);

#pragma omp parallel for schedule(dynamic) default(shared)
    for (int iX = 0; iX < xBins; ++iX)
      for (int iY = 0; iY < yBins; ++iY) 
	for (int iZ = 0; iZ < zBins; ++iZ) {

	  valarray<double> ic(0., galaxy.n_E_gammagrid);
	
	  for (int i = 0; i < galaxy.n_ISRF_components; ++i) 
	    ic += galaxy.IC_iso_emiss[i].d3[iX][iY][iZ].GetSpec();

	  for (int iE = 0; iE < galaxy.n_E_gammagrid; ++iE) {
 
	    gammaRayEmissivity.d3[iX][iY][iZ].s[iE] = ic[iE] + 
	      galaxy.pi0_decay_emiss.d3[iX][iY][iZ].s[iE]*(galaxy.n_HI.d3[iX][iY][iZ].s[0] + 2.*galaxy.n_H2.d3[iX][iY][iZ].s[0] + galaxy.n_HII.d3[iX][iY][iZ].s[0]) + 
	      galaxy.bremss_emiss.d3[iX][iY][iZ].s[iE]*(galaxy.n_HI.d3[iX][iY][iZ].s[0] + 2.*galaxy.n_H2.d3[iX][iY][iZ].s[0]) +
	      galaxy.bremss_ionized_emiss.d3[iX][iY][iZ].s[iE]*galaxy.n_HII.d3[iX][iY][iZ].s[0]; // MeV^-1 cm^-3 s^-1 sr^-1

	  }

	}

    Distribution gammaRayIntensity(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);

    INFO("Calculating gamma-ray intensity. This may take a while ...");

    //#pragma omp parallel for schedule(dynamic) default(shared)
    for (int iX = 0; iX < xBins; ++iX)
      for (int iY = 0; iY < yBins; ++iY)
#pragma omp parallel for schedule(dynamic) default(shared)
	for (int iZ = 0; iZ < zBins; ++iZ) {

	  const double x = galaxy.x[iX], y = galaxy.y[iY], z = galaxy.z[iZ];

	  for (int iXX = 0; iXX < xBins; ++iXX)
	    for (int iYY = 0; iYY < yBins; ++iYY)
	      for (int iZZ = 0; iZZ < zBins; ++iZZ) {
		
		const double xp = galaxy.x[iXX], yp = galaxy.y[iYY], zp = galaxy.z[iZZ];
		
		const double rp2 = (x - xp)*(x - xp) + (y - yp)*(y - yp) + (z - zp)*(z - zp), invRp2 = (rp2 > 0. ? 1./rp2 : 0.);

		const double vol = galaxy.dx*galaxy.dy*galaxy.dz;
		
		// Calculate optical depth from (x, y, z) to (xp, yp, zp) as 
		// function of gamma ray energy
		
		valarray<double> opticalDepth(0., galaxy.n_E_gammagrid);

		// Get direction cosines, integrate along line of sight

		const double rp = sqrt(rp2);
	      
		const int steps = (rp/ds > 1. ? int(rp/ds) : 1);

		const double pX = (rp > 0. ? (xp - x)/rp : 0.);
		const double pY = (rp > 0. ? (yp - y)/rp : 0.);
		const double pZ = (rp > 0. ? (zp - z)/rp : 0.);

		double xx = x, yy = y, zz = z;
		
		for (int i = 0; i < steps; ++i) {
		  
		  xx += pX*ds;
		  yy += pY*ds;
		  zz += pZ*ds;

		  if (zz < galaxy.z_min || zz > galaxy.z_max)
		    continue;

		  if (xx < galaxy.x_min || xx > galaxy.x_max)
		    continue;

		  if (yy < galaxy.y_min || yy > galaxy.y_max)
		    continue;
		  	
		  const double dX = (xx - galaxy.x_min)/(galaxy.x_max - galaxy.x_min);
		  const double dY = (yy - galaxy.y_min)/(galaxy.y_max - galaxy.y_min);
		  const double dZ = (zz - galaxy.z_min)/(galaxy.z_max - galaxy.z_min);

		  const unsigned int iNX = (dX < 1. ? dX*(galaxy.n_xgrid-1) : galaxy.n_xgrid-1);
		  const unsigned int iNY = (dY < 1. ? dY*(galaxy.n_ygrid-1) : galaxy.n_ygrid-1);
		  const unsigned int iNZ = (dZ < 1. ? dZ*(galaxy.n_zgrid-1) : galaxy.n_zgrid-1);
  
		  for (int iG = 0; iG < galaxy.n_E_gammagrid; ++iG)
		    opticalDepth[iG] += galaxy.fPairAbsorptionInvLength.d3[iNX][iNY][iNZ].s[iG]*utl::kpc/utl::cm*ds;

		}
	 
		for (int iG = 0; iG < galaxy.n_E_gammagrid; ++iG) {
		  gammaRayIntensity.d3[iX][iY][iZ].s[iG] += 
		    (std::isnan(gammaRayEmissivity.d3[iXX][iYY][iZZ].s[iG]) ? 0. : gammaRayEmissivity.d3[iXX][iYY][iZZ].s[iG])*invRp2*vol*
		    (utl::kKiloParsec*utl::m/utl::cm)* // MeV^-1 cm^-2 s^-1 sr^-1
		    exp(-opticalDepth[iG]); // Takes care of absorption

		}

	      }
     	    
	}

    // Construct source function throughout Galactic volume

    pairProductionSource.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_pgrid);

    INFO("Calculating pair production source distribution");

#pragma omp parallel for schedule(dynamic) default(shared)
    for (int iX = 0; iX < xBins; ++iX) 
      for (int iY = 0; iY < yBins; ++iY)
	for (int iZ = 0; iZ < zBins; ++iZ) {

	  for (int iE = 0; iE < galaxy.n_pgrid; ++iE) {

	    // Using approximation from Boettcher and Schlickeiser
	    
	    double pairBoettecher = 0;

	    const int iGMin = std::max(0, int(log(gamma[iE]/eps2[0])/eps2Delta));
	    	    
	    for (int iG = iGMin; iG < galaxy.n_E_gammagrid; ++iG) {
	      
	      double sum = 0;
	      
	      if (eps2[iG] > gamma[iE]) {
		
		const double epsMin = eps2[iG]/4./gamma[iE]/(eps2[iG] - gamma[iE]);
		
		if (epsMin <= eps1[targetBins-1]) {
		  
		  const int iTMin = std::max(0, int(log(epsMin/eps1[0])/eps1Delta));
		  
		  for (int iT = iTMin; iT < targetBins; ++iT) {
		    
		    const double fac = targetE[iT]*eps1Delta/eps1[iT]/eps1[iT]*CrossSectionIsotropicBoettecher(eps2[iG], eps1[iT], gamma[iE]); // eV cm^2
		
		    sum += (galaxy.ISRF[0].d3[iX][iY][iZ].s[iT] + galaxy.ISRF[1].d3[iX][iY][iZ].s[iT] + galaxy.ISRF[2].d3[iX][iY][iZ].s[iT])*fac; // eV^-1 cm^-3 * eV cm2 -> cm^-1
		    
		  }
		  
		}
		
	      }
	      
	      const double fac = gammaRayIntensity.d3[iX][iY][iZ].s[iG]*
		galaxy.E_gamma[iG]*eps2Delta/eps2[iG]/eps2[iG]/eps2[iG]; // MeV^-1 cm^-2 s^-1 sr^-1 * MeV
	    
	      pairBoettecher += fac*sum; // cm^-3 s^-1 sr^-1

	    }
	    
	    pairProductionSource.d3[iX][iY][iZ].s[iE] = 
	      pairBoettecher* // cm^-3 s^-1 sr^-1
	      utl::kSpeedOfLight_SI*(utl::m/utl::cm)* // -> cm^-2 s^-2 sr^-1
	      (1./(utl::kElectronMass*(utl::eV/utl::MeV))); // -> MeV^-1 cm^-2 s^-2 sr^-1 -- same units as used in gen_secondary_positron_source.cc
	 
	    assert(pairProductionSource.d3[iX][iY][iZ].s[iE] >= 0.);
   
	  }
	  
	}
	
  }
   
  particle.secondary_source_function = pairProductionSource;

  INFO("Exit");

}
