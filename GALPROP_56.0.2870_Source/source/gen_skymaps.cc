#include "galprop_classes.h"
#include "galprop_internal.h"
#include "los_integration.h"
#include "kappa_free_free.h"

#include <ErrorLogger.h>
#include <Timer.h>
#include <iterator>
#include <algorithm>
#include <cassert>

#include <FullSky.h>

using namespace std;

// Inverse Compton source function assuming isotropic electrons and targets
// Units : photons MeV^-1 cm^2 electron^-1 photon^-1

static double IsotropicICSSourceFunction(const double epsilon2,
					 const double epsilon1,
					 const double gamma) {

  using namespace utl;

  double result = 0.;

  const double p = 4.*epsilon1*gamma;

  if (gamma > 0. &&
      gamma > epsilon2 &&
      epsilon2 <= p*gamma/(1. + p)) {

    const double invP = 1/p;

    const double invGammaEps2 = 1/(gamma - epsilon2);

    const double q = epsilon2*invP*invGammaEps2;///p/(gamma - epsilon2);

    const double pq = p*q;

    const double F =
     2.*q*log(q) + (1. + 2.*q)*(1. - q) + 0.5*pq*pq/(1. + pq)*(1. - q);

    result =
      3./4.*(kThomsonCrossSection_SI*utl::m/utl::cm*utl::m/utl::cm)*1./kElectronMass*1./(utl::eV/utl::MeV)*F/gamma/gamma/epsilon1;

    //cout << "Iso: " << result*kElectronMass << " " << kElectronMass << " " << kThomsonCrossSection_SI*utl::m/utl::cm*utl::m/utl::cm << endl;

  }

  return result;

}

// General inverse Compton source function assuming isotropic electrons
// Units : photons MeV^-1 cm^2 electron^-1 photon^-1 

static double AnisotropicICSSourceFunction(const double epsilon2,
					   const double epsilon1,
					   const double gamma,
					   const double cosZeta) {
  
  using namespace utl;

  double result = 0.;
  
  if (gamma > 1. &&
      gamma > epsilon2) {
    
    const double oneOnGamma = 1./gamma;
    
    const double beta = sqrt(1. - oneOnGamma*oneOnGamma);
    
    const double epsilon1Prime = epsilon1*gamma*(1. + beta*cosZeta);
    
    if (epsilon2 <= 2.*gamma*epsilon1Prime/(1. + 2.*epsilon1Prime)) {
      
      const double epsilon2OnGamma = epsilon2*oneOnGamma;
      
      const double oneOnEpsilon1Prime = 1./epsilon1Prime;
      
      const double gammaMinusEpsilon2 = gamma - epsilon2;
      
      const double F = 2. - 2.*epsilon2OnGamma*(oneOnEpsilon1Prime + 2.) +
        epsilon2OnGamma*epsilon2OnGamma*(oneOnEpsilon1Prime*oneOnEpsilon1Prime +
                                         2*oneOnEpsilon1Prime + 3.) -
        epsilon2OnGamma*epsilon2OnGamma*epsilon2OnGamma;
      
      result = 3./8.*(kThomsonCrossSection_SI/utl::cm/utl::cm)*F/epsilon1/
        gammaMinusEpsilon2/gammaMinusEpsilon2/(kElectronMass*utl::eV/utl::MeV);
      
    }
    
  }

  return result;
    
}

static double BlackBodyNumberDensity(const double energy, const double kT) {   

  using namespace utl;

  const double energyOnPi = energy/utl::kPi,
    constant = 1./(utl::kPlanckReduced/utl::s/utl::eV*utl::kSpeedOfLight_SI/utl::cm);

  return constant*constant*constant*
    energyOnPi*energyOnPi/(exp(energy/kT) - 1.); // eV^-1 cm^-3

}

//Class to integrate the anisotropic emissivity.  Should be used only for low healpix orders otherwise takes long time.
Galprop::AnisoICFunction::AnisoICFunction(const Galprop &gp) : 
   cameraX(gp.galdef.fCameraLocation[0]),
   cameraY(gp.galdef.fCameraLocation[1]),
   cameraZ(gp.galdef.fCameraLocation[2]),
   fgp(gp),
   cosThetaBins(500) // Set this to the number for precomputing the anisotropic cross section. Accuracy is best for 500, 5000 (need to follow up once working).
{

      const double dCosTheta = 2./cosThetaBins;

      //Distribution electrons;

      // identify the electrons/positrons IMOS20030217
      if (2 == fgp.galdef.n_spatial_dimensions) 
	electrons.init(fgp.gcr[0].n_rgrid, fgp.gcr[0].n_zgrid, fgp.gcr[0].n_pgrid);
      
      if (3 == fgp.galdef.n_spatial_dimensions) 
	electrons.init(fgp.gcr[0].n_xgrid, fgp.gcr[0].n_ygrid, fgp.gcr[0].n_zgrid, fgp.gcr[0].n_pgrid);
      
      electrons = 0.;
      
      int iE, ielectrons = -1;
      
      for (iE = 0, ielectrons = -1; iE < fgp.n_species; ++iE)  
	if (100 == 100*abs(fgp.gcr[iE].Z) + fgp.gcr[iE].A) {
	  
	  ielectrons = iE;
	  electrons += fgp.gcr[ielectrons].cr_density;
	  
	  std::ostringstream buf;
	  buf << "CR " << fgp.gcr[ielectrons].name << " found as species #" << ielectrons;
	  INFO(buf.str());
	  
	}
      
      if (-1 == ielectrons) {
	
	std::ostringstream buf;
	buf << "CR electrons/positrons not found.";
	INFO(buf.str());
	electrons.delete_array();
        throw(std::invalid_argument("CR electrons/positrons not found"));
	
      }

      nEGammaBins = fgp.galaxy.n_E_gammagrid;

      gammaE.resize(nEGammaBins, 0.);
      valarray<double> eps2(0., nEGammaBins);

      for (size_t i = 0; i < gammaE.size(); ++i)
	gammaE[i] = fgp.galaxy.E_gamma[i];

      eps2 = gammaE/Mele;

      electronBins = fgp.gcr[ielectrons].n_pgrid;

      electronE.resize(electronBins, 0.);
      valarray<double> gamma(0., electronBins);

      for (size_t i = 0; i < gamma.size(); ++i)
	electronE[i] = fgp.gcr[ielectrons].Ekin[i];

      gamma = electronE/Mele;

      targetBins = fgp.galaxy.ISRF[0].n_pgrid;

      targetE.resize(targetBins, 0.);
      valarray<double> eps1(0., targetBins), targetFreq(0., targetBins);

      for (size_t i = 0; i < targetE.size(); ++i) {

	targetFreq[i] = fgp.galaxy.nu_ISRF[i];
	targetE[i] = h_planck*fgp.galaxy.nu_ISRF[i]*erg_to_eV;  // target photon energy in eV

      }

      eps1 = targetE*1e-6/Mele;

      INFO("Precalculating cross section data");
	
      isoCrossSection.resize(nEGammaBins*electronBins*targetBins, 0.);
      anisoCrossSection.resize(nEGammaBins*electronBins*targetBins*cosThetaBins, 0.);

      for (size_t i = 0; i < gammaE.size(); ++i)
	for (size_t j = 0; j < gamma.size(); ++j)
	  for (size_t k = 0; k < targetE.size(); ++k) {
	    
	    const size_t isoIndex = i*electronBins*targetBins + j*targetBins + k;
	    
	    isoCrossSection[isoIndex] = IsotropicICSSourceFunction(eps2[i], eps1[k], gamma[j]);

	    for (size_t l = 0; l < cosThetaBins; ++l) {
	      
	      const size_t anisoIndex = i*electronBins*targetBins*cosThetaBins + j*targetBins*cosThetaBins + k*cosThetaBins + l;
	      
	      const double cosTheta = -1. + (l + 0.5)*dCosTheta;
	      
	      anisoCrossSection[anisoIndex] = AnisotropicICSSourceFunction(eps2[i], eps1[k], gamma[j], cosTheta);
	      
	    }

	  }   

      // Construct the radiation field

      /*
      const std::string fitsDirectory = fgp.configure.fFITSDataDirectory;
      const std::string isrfFilename = fgp.galdef.ISRF_file;
      const std::string filename = fitsDirectory + isrfFilename;
      
      ostringstream buf1;
      buf1 << "Reading ISRF from " << filename;
      INFO(buf1.str());
      
      if ( filename != fgp.galaxy.fISRFloadedfile ) {
         delete fgp.galaxy.fISRFrf;
         fgp.galaxy.fISRFrf = new rf::RadiationField(filename, targetFreq, fgp.galdef.ISRF_healpixOrder);
         fgp.galaxy.fISRFloadedfile = filename;
      }
      */

      rf::RadiationField& rf = *fgp.galaxy.fISRFrf;

      const unsigned int rBins3D = std::max(fgp.galaxy.n_ygrid, fgp.galaxy.n_xgrid)/2; 
      
      rBins = (2 == fgp.gcr[0].n_spatial_dimensions ? fgp.galaxy.n_rgrid : utl::kSqrtTwo*rBins3D + 1);
      
      assert (fgp.galaxy.n_ISRF_components <= 3);

      rGrid.resize(rBins, 0.);

      const double rMax3D = std::max(fgp.galaxy.y_max, fgp.galaxy.x_max);

      rMax = (2 == fgp.gcr[0].n_spatial_dimensions ? fgp.galaxy.r_max : utl::kSqrtTwo*rMax3D); 
      //	   sqrt(galaxy.x_max*galaxy.x_max + galaxy.y_max*galaxy.y_max));

      for (size_t iR = 0; iR < rBins; ++iR)
	rGrid[iR] = (2 == fgp.gcr[0].n_spatial_dimensions ? fgp.galaxy.r[iR] : rMax/rBins*iR);
      
      std::ostringstream buf2;
      buf2 << "Constructing ISRF angular distribution for healpix order " << fgp.galdef.ISRF_healpixOrder;
      INFO(buf2.str());
      
      cmbSkymap.Resize(fgp.galdef.ISRF_healpixOrder, targetFreq);
      
      cmbNumberDensity.resize(targetFreq.size(), 0.);
      
      for (size_t i = 0; i < targetE.size(); ++i) {
	
	cmbNumberDensity[i] = BlackBodyNumberDensity(targetE[i], 2.735*utl::kBoltzmann_SI/utl::e_SI); // eV^-1 cm^-3 
		
      }
      
      for (int i = 0; i < cmbSkymap.Npix(); ++i)
	for (size_t iT = 0; iT < targetE.size(); ++iT)
	  cmbSkymap[i][iT] = (cmbNumberDensity[iT]*1./utl::kFourPi*cmbSkymap.solidAngle()); // Convert to eV^-1 cm^-3 * dSA/4pi per pixel

      dirTarget.resize(cmbSkymap.Npix(), vec3(0, 0, 0));
      
      for (size_t i = 0; i < dirTarget.size(); ++i)
	dirTarget[i] = cmbSkymap.pix2coord(i).healpixAng().to_vec3();

      // Just initialise the cache -- don't care about the result

      rf::RadiationField::ThreeVector pos(0., 0., 0.);
      
      rf.GetSkymap(pos, rf::RadiationField::TOTAL, fgp.galdef.ISRF_healpixOrder);

}

std::valarray<double> Galprop::AnisoICFunction::operator () ( const double x, const double y, const double z, const vec3 &cdir ) const {
   std::valarray<double> output(0.,nEGammaBins*6);

   const double theta = std::atan2(y, x);
   //Need to rotate the camera direction by -theta and flip the sign on z as we should be pointing in the other direction
   const vec3 dir( cdir.x*std::cos(theta)+cdir.y*std::sin(theta), -cdir.x*std::sin(theta)+cdir.y*std::cos(theta), -cdir.z ), dir3D(cdir.x, cdir.y, -cdir.z);

   //Pre calculate cosZeta (This could be put in the constructor for several values of dir like in the old routine)
   std::valarray<double> cosZeta(0., cmbSkymap.Npix());
   for (size_t i = 0; i < cosZeta.size(); ++i) {
      cosZeta[i] = dotprod(dirTarget[i], dir);
      if (std::fabs(cosZeta[i]) < 1e-10)
         cosZeta[i] = 0;
   }
/*
   //Find the index into the gamma skymap
   const unsigned int gammaPix = hpTest.vec2pix(dir);
*/
      
   //Do the actual integral, linear interpolation assuming the bins are evenly distributed in a linear fashion
   const double dz = (z - fgp.galaxy.z_min) / (fgp.galaxy.z_max - fgp.galaxy.z_min);
   const size_t izl = dz <= 0. ? 0 : dz < 1. ? dz*(fgp.galaxy.n_zgrid-1) : fgp.galaxy.n_zgrid-2;
   const size_t izu = izl + 1;

   //We don't want to extrapolate, so use upper/lower boundary in case the galaxy is smaller than the integration region (highly unlikely)
   const double rz = (fgp.galaxy.z[izu] - z) / (fgp.galaxy.z[izu] - fgp.galaxy.z[izl]);
   const double lzf = rz < 0. ? 0. : rz > 1. ? 1. : rz;
   const double uzf = 1-lzf;

   //Since ISRF is on a R,z grid, we need this for both cases
   const double r = std::sqrt(x*x+y*y);

   //galaxy.r_min is always 0, so implicitly assumed here
   const double dr = r / rMax;
   const size_t irl = dr <= 0. ? 0 : dr < 1. ? dr*(rBins-1) : rBins-2;
   const size_t iru = irl + 1;

   const double rr = (rGrid[iru] - r) / (rGrid[iru] - rGrid[irl]);
   const double lrf = rr < 0. ? 0. : rr > 1. ? 1. : rr;
   const double urf = 1-lrf;

   /*const size_t index1 = irl*fgp.galaxy.n_zgrid + izl;
   const size_t index2 = iru*fgp.galaxy.n_zgrid + izl;
   const size_t index3 = irl*fgp.galaxy.n_zgrid + izu;
   const size_t index4 = iru*fgp.galaxy.n_zgrid + izu;

   const Skymap<double>& opt1 = optAngDist[index1];
   const Skymap<double>& ir1 = irAngDist[index1];
   const Skymap<double>& opt2 = optAngDist[index2];
   const Skymap<double>& ir2 = irAngDist[index2];
   const Skymap<double>& opt3 = optAngDist[index3];
   const Skymap<double>& ir3 = irAngDist[index3];
   const Skymap<double>& opt4 = optAngDist[index4];
   const Skymap<double>& ir4 = irAngDist[index4];
   */
   rf::RadiationField& rf = *fgp.galaxy.fISRFrf;

   const rf::RadiationField::ThreeVector pos(x, y, z);

   const Skymap<double> optIntensity =  (rf.Is3D() ? rf.GetSkymap(pos, rf::RadiationField::OPTICAL, fgp.galdef.ISRF_healpixOrder) : (rf.GetSkymap(pos, rf::RadiationField::DIRECT, fgp.galdef.ISRF_healpixOrder) +  rf.GetSkymap(pos, rf::RadiationField::SCATTERED, fgp.galdef.ISRF_healpixOrder)))*cmbSkymap.solidAngle();

   const Skymap<double> irIntensity = (rf.Is3D() ? rf.GetSkymap(pos, rf::RadiationField::INFRARED, fgp.galdef.ISRF_healpixOrder) : (rf.GetSkymap(pos, rf::RadiationField::TRANSIENT, fgp.galdef.ISRF_healpixOrder) +  rf.GetSkymap(pos, rf::RadiationField::THERMAL, fgp.galdef.ISRF_healpixOrder)))*cmbSkymap.solidAngle();

   std::valarray<double> optNumberDensity(0., targetE.size()), irNumberDensity(0., targetE.size());

   for (size_t iTarget = 0; iTarget < targetE.size(); ++iTarget) {

     optNumberDensity[iTarget] = optIntensity.sum(iTarget);
     irNumberDensity[iTarget] = irIntensity.sum(iTarget);
  
   }

   //Need different interpolation as well for electrons in 3D case
   double lxf(0.), uxf(0.), lyf(0.), uyf(0.);
   size_t ixl(0), ixu(0), iyl(0), iyu(0);
   if (3 == fgp.gcr[0].n_spatial_dimensions) {
      const double dx = (x - fgp.galaxy.x_min) / (fgp.galaxy.x_max - fgp.galaxy.x_min);
      ixl = dx <= 0. ? 0 : dx < 1. ? dx*(fgp.galaxy.n_xgrid-1) : fgp.galaxy.n_xgrid-2;
      ixu = ixl + 1;

      const double rx = (fgp.galaxy.x[ixu] - x) / (fgp.galaxy.x[ixu] - fgp.galaxy.x[ixl]);
      lxf = rx < 0. ? 0. : rx > 1. ? 1. : rx;
      uxf = 1-lxf;

      const double dy = (y - fgp.galaxy.y_min) / (fgp.galaxy.y_max - fgp.galaxy.y_min);
      iyl = dy <= 0. ? 0 : dy < 1. ? dy*(fgp.galaxy.n_ygrid-1) : fgp.galaxy.n_ygrid-2;
      iyu = iyl + 1;

      const double ry = (fgp.galaxy.y[iyu] - y) / (fgp.galaxy.y[iyu] - fgp.galaxy.y[iyl]);
      lyf = ry < 0. ? 0. : ry > 1. ? 1. : ry;
      uyf = 1-lyf;
   }

   const double logTargetE1over0 = log(targetE[1]/targetE[0]);

   for (size_t iElectron = 0; iElectron < electronBins; ++iElectron) {

      // All grid points for the bi-/trilinear interpolation
      
      double elecSum(0.); 
      if (2 == fgp.gcr[0].n_spatial_dimensions) {
         elecSum =
            (electrons.d2[irl][izl].s[iElectron]*lrf*lzf +
             electrons.d2[iru][izl].s[iElectron]*urf*lzf +
             electrons.d2[irl][izu].s[iElectron]*lrf*uzf +
             electrons.d2[iru][izu].s[iElectron]*urf*uzf)*
            electronE[iElectron]*log(electronE[1]/electronE[0]);
      }
      if (3 == fgp.gcr[0].n_spatial_dimensions) {
         elecSum =
            (electrons.d3[ixl][iyl][izl].s[iElectron]*lxf*lyf*lzf +
             electrons.d3[ixu][iyl][izl].s[iElectron]*uxf*lyf*lzf +
             electrons.d3[ixl][iyu][izl].s[iElectron]*lxf*uyf*lzf +
             electrons.d3[ixu][iyu][izl].s[iElectron]*uxf*uyf*lzf +
             electrons.d3[ixl][iyl][izu].s[iElectron]*lxf*lyf*uzf +
             electrons.d3[ixu][iyl][izu].s[iElectron]*uxf*lyf*uzf +
             electrons.d3[ixl][iyu][izu].s[iElectron]*lxf*uyf*uzf +
             electrons.d3[ixu][iyu][izu].s[iElectron]*uxf*uyf*uzf)*
            electronE[iElectron]*log(electronE[1]/electronE[0]);
      }
         
      for (size_t iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {
        
         double anisoSumOpt = 0, isoSumOpt = 0;
         double anisoSumIR = 0, isoSumIR = 0;
         double anisoSumCMB = 0, isoSumCMB = 0; // For checking
      
         for (size_t iTarget = 0; iTarget < targetBins; ++iTarget) {
      
           const size_t isoIndex = iEGamma*electronBins*targetBins + iElectron*targetBins + iTarget;
      
           // This is the isotropic cross section summation
           
           const double xsIso = isoCrossSection[isoIndex]*targetE[iTarget]*logTargetE1over0;
      
           isoSumOpt += xsIso*optNumberDensity[iTarget];
	   //             (optISRF[index1][iTarget]*lrf*lzf +
	   //optISRF[index2][iTarget]*urf*lzf +
	   //optISRF[index3][iTarget]*lrf*uzf +
	   //optISRF[index4][iTarget]*urf*uzf);
          
           isoSumIR += xsIso*irNumberDensity[iTarget];
	   //             (irISRF[index1][iTarget]*lrf*lzf +
	   //irISRF[index2][iTarget]*urf*lzf +
	   //irISRF[index3][iTarget]*lrf*uzf +
	   //irISRF[index4][iTarget]*urf*uzf);
       				
           isoSumCMB += xsIso*cmbNumberDensity[iTarget];
      
           // Inner loop over the target skymaps
      
           double sumSkymapOpt = 0, sumSkymapIR = 0, sumSkymapCMB = 0;
      
           for (int iSkymap = 0; iSkymap < cmbSkymap.Npix(); ++iSkymap) {
             //const unsigned int cosZetaIndex = iSkymap*hpTest.Npix() + gammaPix;

	     const auto cosTheta = (rf.Is3D() ? dotprod(dirTarget[iSkymap], dir3D) : cosZeta[iSkymap]);
	     
	     const auto cosThetaIndex = ((1. + cosTheta < 2.) ? int(0.5*(1. + cosTheta)*cosThetaBins) : cosThetaBins-1);

	     //(rf.Is3D() ? 0 : ((1. + cosZeta[iSkymap]) < 2. ? int(0.5*(1. + cosZeta[iSkymap])*cosThetaBins) : cosThetaBins-1)); 
	     
             //const unsigned int cosThetaIndex = 
	     //((1. + cosZeta[iSkymap]) < 2. ? int(0.5*(1. + cosZeta[iSkymap])*cosThetaBins) : cosThetaBins-1); 
       	       //((1. + cosZetaTest[cosZetaIndex]) < 2. ? int(0.5*(1. + cosZetaTest[cosZetaIndex])*cosThetaBins) : cosThetaBins-1); 
      
             const size_t anisoIndex = iEGamma*electronBins*targetBins*cosThetaBins + iElectron*targetBins*cosThetaBins + iTarget*cosThetaBins + cosThetaIndex;
           
             sumSkymapOpt += anisoCrossSection[anisoIndex]*optIntensity[iSkymap][iTarget];
	     //	       (opt1[iSkymap][iTarget]*lrf*lzf +
	     //opt2[iSkymap][iTarget]*urf*lzf +
	     //opt3[iSkymap][iTarget]*lrf*uzf +
	     //opt4[iSkymap][iTarget]*urf*uzf);
	     
             sumSkymapIR += anisoCrossSection[anisoIndex]*irIntensity[iSkymap][iTarget];
	     //	       (ir1[iSkymap][iTarget]*lrf*lzf +
	     //ir2[iSkymap][iTarget]*urf*lzf +
	     //ir3[iSkymap][iTarget]*lrf*uzf +
	     //ir4[iSkymap][iTarget]*urf*uzf);
             
             sumSkymapCMB += anisoCrossSection[anisoIndex]*
	       cmbSkymap[iSkymap][iTarget];
	     
           }
            
           anisoSumOpt += sumSkymapOpt*targetE[iTarget]*logTargetE1over0;
           anisoSumIR += sumSkymapIR*targetE[iTarget]*logTargetE1over0;
           anisoSumCMB += sumSkymapCMB*targetE[iTarget]*logTargetE1over0;
      
         }
      
         // All grid points for the bi-/trilinear interpolation
         
         /*double elecSum = 0; 
         if (2 == fgp.gcr[0].n_spatial_dimensions) {
            elecSum =
               (electrons.d2[irl][izl].s[iElectron]*lrf*lzf +
                electrons.d2[iru][izl].s[iElectron]*urf*lzf +
                electrons.d2[irl][izu].s[iElectron]*lrf*uzf +
                electrons.d2[iru][izl].s[iElectron]*urf*uzf)*
               electronE[iElectron]*log(electronE[1]/electronE[0]);
         }
         if (3 == fgp.gcr[0].n_spatial_dimensions) {
            elecSum =
               (electrons.d3[ixl][iyl][izl].s[iElectron]*lxf*lyf*lzf +
                electrons.d3[ixu][iyl][izl].s[iElectron]*uxf*lyf*lzf +
                electrons.d3[ixl][iyu][izl].s[iElectron]*lxf*uyf*lzf +
                electrons.d3[ixu][iyu][izl].s[iElectron]*uxf*uyf*lzf +
                electrons.d3[ixl][iyl][izu].s[iElectron]*lxf*lyf*uzf +
                electrons.d3[ixu][iyl][izu].s[iElectron]*uxf*lyf*uzf +
                electrons.d3[ixl][iyu][izu].s[iElectron]*lxf*uyf*uzf +
                electrons.d3[ixu][iyu][izu].s[iElectron]*uxf*uyf*uzf)*
               electronE[iElectron]*log(electronE[1]/electronE[0]);
         }
         */
         output[0*nEGammaBins + iEGamma] += isoSumOpt*elecSum;
         output[1*nEGammaBins + iEGamma] += isoSumIR*elecSum;
         output[2*nEGammaBins + iEGamma] += isoSumCMB*elecSum;
      
         output[3*nEGammaBins + iEGamma] += anisoSumOpt*elecSum;
         output[4*nEGammaBins + iEGamma] += anisoSumIR*elecSum;
         output[5*nEGammaBins + iEGamma] += anisoSumCMB*elecSum;
      
         //cout << iEGamma << " " << iElectron << " " << isoSumOpt << " " << isoSumIR << " " << isoSumCMB << " " << anisoSumOpt << " " << anisoSumIR << " " << anisoSumCMB << " " << elecSum << endl;
      
      }

   }

   return output;
}

//Since the pixel code will be identical for healpix and CAR, just have it a function.
static void calculatePixel(const double l, const double b, const double dl, 
      const double cameraTheta, const double cameraR,
      const vector<double> &Rbins,
      const std::unique_ptr< SM::LOSintegrator<double> > &gasInt, 
      const std::unique_ptr< SM::LOSintegrator< valarray<double> > > &emissInt, 
      const vector< std::unique_ptr< SM::LOSfunction<double> > > &gasFuncs, 
      const vector< std::unique_ptr< SM::LOSfunction< valarray<double> > > > &emissFuncs,
      const std::unique_ptr< SM::LOSfunction<valarray<double> > > &absorption,
      vector< vector< double > >&gas, vector< vector< valarray<double> > > &emiss);

//Returns the index for the Annulus at radius R
static size_t findAnnulus(const vector<double> &Rbins, const double R) {
   size_t il(0);
   size_t iu(Rbins.size());
   while (iu-il>1){
      size_t i = (il+iu)/2;
      if ( R < Rbins[i] ){
         iu = i;
      } else {
         il = i;
      }
   }
   return il;
}

void Galprop::gen_skymaps() 
{

   INFO( "Entry" );

   //Use the LOS integrators to integrate the maps.  Need two integrators, one for the gasMaps, other for emissivity
   std::unique_ptr< SM::LOSintegrator<double> > gasInt;
   std::unique_ptr< SM::LOSintegrator< std::valarray<double> > > emissInt;

   //A vector with functions to integrate
   vector< std::unique_ptr< SM::LOSfunction<double> > > gasFuncs;
   vector< std::unique_ptr< SM::LOSfunction< valarray<double> > > > emissFuncs;

   //Create a vector with the annuli boundaries, only initialized for pi0 and bremsstrahlung
   vector<double> Rbins;
   if (galdef.pi0_decay || galdef.bremss) 
   {

      const std::string cor("COR"); // IMOS20080114
      const std::string hir("HIR"); // IMOS20080114

      TIME_SUBROUTINE(read_gas_maps, hir); // IMOS20080114
      TIME_SUBROUTINE(read_gas_maps, cor); // IMOS20080114

      copy (&galaxy.R_bins[0], &galaxy.R_bins[0]+galaxy.n_Ring, back_inserter(Rbins));
      Rbins.push_back(galaxy.R_bins[2*galaxy.n_Ring-1]);

      for (size_t i = 0; i < Rbins.size(); ++i)
      {
         std::ostringstream buf;
         buf<<"Rbins["<<i<<"] = "<<Rbins[i];
         INFO(buf.str());
      }

      //Gas functions only needed for bremss and pi0 decay
      gasFuncs.emplace_back(new GasFunction("HI", 0, *this));
      gasFuncs.emplace_back(new GasFunction("H2", 0, *this));
      gasFuncs.emplace_back(new GasFunction("CO", 0, *this));
   }

   //Initialize the integrators according to geometry
   if (2 == galdef.n_spatial_dimensions) 
   {
      gasInt.reset(new SM::LOSintegrator<double>(galaxy.r_max, galaxy.z_min, galaxy.z_max, Rbins, galdef.fCameraLocation, galdef.LoS_step, galdef.los_integration_accuracy, galdef.LoS_minStep));
      emissInt.reset(new SM::LOSintegrator<std::valarray<double> >(galaxy.r_max, galaxy.z_min, galaxy.z_max, Rbins, galdef.fCameraLocation, galdef.LoS_step, galdef.los_integration_accuracy, galdef.LoS_minStep));
   } 
   else if (3 == galdef.n_spatial_dimensions)  
   {
      gasInt.reset(new SM::LOSintegrator<double>(galaxy.x_min, galaxy.x_max, galaxy.y_min, galaxy.y_max, galaxy.z_min, galaxy.z_max, Rbins, galdef.fCameraLocation, galdef.LoS_step, galdef.los_integration_accuracy, galdef.LoS_minStep));
      emissInt.reset(new SM::LOSintegrator<std::valarray<double> >(galaxy.x_min, galaxy.x_max, galaxy.y_min, galaxy.y_max, galaxy.z_min, galaxy.z_max, Rbins, galdef.fCameraLocation, galdef.LoS_step, galdef.los_integration_accuracy, galdef.LoS_minStep));
   }

   //Calculate the angle the camera makes to the global grid (from positive x axis to x,y, goes from 0 to 360)
   const double cX = galdef.fCameraLocation[0];
   const double cY = galdef.fCameraLocation[1];
   const double cameraTheta = std::atan2(cY,cX)*180./M_PI;
   const double cameraR = sqrt(cX*cX + cY+cY);

   // Add absorption if pair production is calculated
   std::unique_ptr< SM::LOSfunction<std::valarray<double> > > absorption;
   if ( galdef.pair_production ) 
   {
      std::vector<double> z(&galaxy.z[0], &galaxy.z[0]+galaxy.n_zgrid);

      if (2 == galaxy.n_spatial_dimensions) 
      {
         std::vector<double> r(&galaxy.r[0], &galaxy.r[0]+galaxy.n_rgrid);
         absorption.reset( new DistributionFunction(galaxy.fPairAbsorptionInvLength, z, r));

      } 
      else if (3 == galaxy.n_spatial_dimensions) 
      {
         std::vector<double> x(&galaxy.x[0], &galaxy.x[0]+galaxy.n_xgrid);
         std::vector<double> y(&galaxy.y[0], &galaxy.y[0]+galaxy.n_ygrid);
         absorption.reset( new DistributionFunction(galaxy.fPairAbsorptionInvLength, z, x, y));

      } 
      else 
      {

	ERROR("Absorption calculation called with dimensions != 2 or 3");
	exit(-1);

      }
   }

   //Add the correct functions, depending on what was asked in the galdef file and initialize the skymaps
   if (galdef.pi0_decay) 
   {

      emissFuncs.emplace_back(new GasEmissFunction("PION", "HI", *this));
      emissFuncs.emplace_back(new GasEmissFunction("PION", "H2", *this));
      emissFuncs.emplace_back(new GasEmissFunction("PION", "HII", *this));

      if (3 == galdef.skymap_format) 
      { // HEALPix format

         galaxy.pi0_decay_hp_skymap = galaxy.createFullSky(galdef.healpix_order);

         if (2 == galdef.gamma_rays) 
         { // split according to rings

            galaxy.pi0_decay_H2R_hp_skymap.reserve(galaxy.n_Ring);
            galaxy.pi0_decay_HIR_hp_skymap.reserve(galaxy.n_Ring);
            galaxy.pi0_decay_HIIR_hp_skymap.reserve(galaxy.n_Ring);

            //Do sparse sky for inner Galaxy, outer Galaxy, and all H2R maps
            for (int iR = 0; iR < galaxy.n_Ring; ++iR) 
            { 

               galaxy.pi0_decay_H2R_hp_skymap.emplace_back(galaxy.createSparseSky(galdef.healpix_order));

               if (galaxy.R_bins[iR] < Rsun && galaxy.R_bins[iR+1] > Rsun) {
                  galaxy.pi0_decay_HIR_hp_skymap.emplace_back(galaxy.createFullSky(galdef.healpix_order));
                  galaxy.pi0_decay_HIIR_hp_skymap.emplace_back(galaxy.createFullSky(galdef.healpix_order));
               } else {
                  galaxy.pi0_decay_HIR_hp_skymap.emplace_back(galaxy.createSparseSky(galdef.healpix_order));
                  galaxy.pi0_decay_HIIR_hp_skymap.emplace_back(galaxy.createSparseSky(galdef.healpix_order));
               }

            }

         }

         INFO("Starting pion-decay integration");

#pragma omp parallel for schedule(dynamic) default(shared)
         for ( int ii = 0; ii < galaxy.pi0_decay_hp_skymap->Npix(); ++ii ) {

           if ( ii % (galaxy.pi0_decay_hp_skymap->Npix()/100) == 0) {
              std::ostringstream buf;
              buf << "Beginning pixel " << ii << " (" << (100*ii)/galaxy.pi0_decay_hp_skymap->Npix() << "%)";
              INFO(buf.str());
           }
	   
	   SM::Coordinate co = galaxy.pi0_decay_hp_skymap->GetCoordinate ( ii );
           double l, b;
           co.getCoordinates(l, b, SM::CoordSys::GAL);
           l *= utl::kConvertRadiansToDegrees;
           b *= utl::kConvertRadiansToDegrees;
	   
	   vector< vector<double> > gas;
	   vector< vector< std::valarray<double> > > emiss;
	   
	   calculatePixel(l, b, 90./galaxy.pi0_decay_hp_skymap->Nside(), cameraTheta, cameraR, Rbins, gasInt, emissInt, gasFuncs, emissFuncs, absorption, gas, emiss);
	   
	   //Add the pixels to the maps
	   
	   //Correct the emissivities with the real gas maps
	   for (size_t iR = 0; iR < emiss[0].size(); ++iR) 
	     {

               //No need to go further if there is nothing to add
               if ( emiss[0][iR].max() <= 0 && emiss[1][iR].max() <= 0 && emiss[2][iR].max() <= 0 ) continue;

               if ( galaxy.renormGas ) 
               {
                  //Calculate the average X_CO factor to multiply the CO map with.
                  double av_X_CO;
                  if ( gas[2][iR] > 0 ) 
                  {
                     av_X_CO = gas[1][iR]/gas[2][iR];
                  } 
                  else 
                  {
                     const double R_ring = ( galaxy.R_bins[iR] +galaxy.R_bins[galaxy.n_Ring+iR] ) /2.;
                     av_X_CO = fX_CO ( R_ring );
                  }

                  //The weights to correct the emissivities
                  const double w_HI = gas[0][iR] > 0 ? galaxy.hpHIR[iR]->GetValue(ii,0) / gas[0][iR] : 1;
                  const double w_H2 = gas[1][iR] > 0 ? 2*av_X_CO*galaxy.hpCOR[iR]->GetValue(ii,0) / gas[1][iR] : 1;

                  emiss[0][iR] *= w_HI;
                  emiss[1][iR] *= w_H2;
               }

               for (size_t iE=0; iE < emiss[0][iR].size(); ++iE) 
                  galaxy.pi0_decay_hp_skymap->GetReference(ii, iE) += emiss[0][iR][iE] + emiss[1][iR][iE] + emiss[2][iR][iE];

               //Split into rings
               if (galdef.gamma_rays==2 && iR < size_t(galaxy.n_Ring) ) 
               {
                  if (emiss[0][iR][0] > 0) 
#pragma omp critical (pi0_decay_HIR_hp_skymap)
                     for (size_t iE=0; iE < emiss[0][iR].size(); ++iE) 
                        galaxy.pi0_decay_HIR_hp_skymap[iR]->GetReference(ii,iE) += emiss[0][iR][iE];
                  if (emiss[1][iR][0] > 0) 
#pragma omp critical (pi0_decay_H2R_hp_skymap)
                     for (size_t iE=0; iE < emiss[1][iR].size(); ++iE) 
                        galaxy.pi0_decay_H2R_hp_skymap[iR]->GetReference(ii,iE) += emiss[1][iR][iE];
                  if (emiss[2][iR][0] > 0) 
#pragma omp critical (pi0_decay_HIIR_hp_skymap)
                     for (size_t iE=0; iE < emiss[2][iR].size(); ++iE) 
                        galaxy.pi0_decay_HIIR_hp_skymap[iR]->GetReference(ii,iE) += emiss[2][iR][iE];
               }

            }
         }

         //Store the skymaps and clear the memory
         TIME_SUBROUTINE(store_pi0_decay_skymap);

         if (2 == galdef.gamma_rays) 
         {                           //AWS20050302

            TIME_SUBROUTINE(store_pi0_decay_H2R_skymap);
            TIME_SUBROUTINE(store_pi0_decay_HIR_skymap);
            TIME_SUBROUTINE(store_pi0_decay_HII_skymap);

         }

         galaxy.pi0_decay_hp_skymap.reset(nullptr);

         galaxy.pi0_decay_H2R_hp_skymap.resize(0);
         galaxy.pi0_decay_HIR_hp_skymap.resize(0);
         galaxy.pi0_decay_HIIR_hp_skymap.resize(0);

      } 
      else 
      {

         galaxy.pi0_decay_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);

         if (2 == galdef.gamma_rays) 
         {

            galaxy.pi0_decay_H2R_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_Ring, galaxy.n_E_gammagrid);

            galaxy.pi0_decay_HIR_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_Ring, galaxy.n_E_gammagrid);

            galaxy.pi0_decay_HIIR_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_Ring, galaxy.n_E_gammagrid);

         }

         INFO("Starting pion-decay integration");
#pragma omp parallel for schedule(dynamic) default(shared)
         for ( int i_long=0; i_long<galaxy.n_long; i_long++ ) 
         {
            if ( i_long % (galaxy.n_long/100) == 0) 
            {
               ostringstream buf;
               buf << "Beginning longitude row " << i_long << " (" << (100*i_long)/galaxy.n_long << "%)";
               INFO(buf.str());
            }
            for ( int i_lat =0; i_lat <galaxy.n_lat; i_lat++ ) 
            {
               const double l=galaxy.long_min + i_long*galaxy.d_long;
               const double b=galaxy. lat_min + i_lat *galaxy.d_lat ;

               vector< vector<double> > gas;
               vector< vector< std::valarray<double> > > emiss;

               calculatePixel(l, b, 360./galaxy.n_long, cameraTheta, cameraR, Rbins, gasInt, emissInt, gasFuncs, emissFuncs, absorption, gas, emiss);

               //Add the pixels to the maps

               //Correct the emissivities with the real gas maps
               for (size_t iR = 0; iR < emiss[0].size(); ++iR) 
               {

                  //No need to go further if there is nothing to add
                  if ( galaxy.renormGas && (gas[0][iR] <= 0 || gas[1][iR] <= 0) ) continue;

                  if ( galaxy.renormGas ) 
                  {
                     //Calculate the average X_CO factor to multiply the CO map with.
                     double av_X_CO;
                     if ( gas[2][iR] > 0 ) 
                     {
                        av_X_CO = gas[1][iR]/gas[2][iR];
                     } 
                     else 
                     {
                        const double R_ring = ( galaxy.R_bins[iR] +galaxy.R_bins[galaxy.n_Ring+iR] ) /2.;
                        av_X_CO = fX_CO ( R_ring );
                     }

                     //The weights to correct the emissivities
                     const double w_HI = galaxy.HIR.d3[i_long][i_lat][iR].s[0] / gas[0][iR];
                     const double w_H2 = 2*av_X_CO*galaxy.COR.d3[i_long][i_lat][iR].s[0] / gas[1][iR];

                     emiss[0][iR] *= w_HI;
                     emiss[1][iR] *= w_H2;
                  }

                  for ( int iE=0; iE < galaxy.n_E_gammagrid; ++iE )
                     galaxy.pi0_decay_skymap.d2[i_long][i_lat].s[iE] += emiss[0][iR][iE] + emiss[1][iR][iE] + emiss[2][iR][iE];

                  //Split into rings
                  if (galdef.gamma_rays==2 && iR < size_t(galaxy.n_Ring) ) 
                  {

                     for ( int iE=0; iE < galaxy.n_E_gammagrid; ++iE ) 
                     {

                        galaxy.pi0_decay_HIR_skymap.d3[i_long][i_lat][iR].s[iE] += emiss[0][iR][iE];
                        galaxy.pi0_decay_H2R_skymap.d3[i_long][i_lat][iR].s[iE] += emiss[1][iR][iE];
                        galaxy.pi0_decay_HIIR_skymap.d3[i_long][i_lat][iR].s[iE] += emiss[2][iR][iE];

                     }
                  }

               }
            }
         }

         //Store the skymaps and clear the memory
         TIME_SUBROUTINE(store_pi0_decay_skymap);

         if (2 == galdef.gamma_rays) 
         {                           //AWS20050302

            TIME_SUBROUTINE(store_pi0_decay_H2R_skymap);
            TIME_SUBROUTINE(store_pi0_decay_HIR_skymap);
            TIME_SUBROUTINE(store_pi0_decay_HII_skymap);

         }

         galaxy.pi0_decay_skymap.delete_array();  //Gulli20070810

         galaxy.pi0_decay_H2R_skymap.delete_array();  //Gulli20070810
         galaxy.pi0_decay_HIR_skymap.delete_array();  //Gulli20070810
         galaxy.pi0_decay_HIIR_skymap.delete_array();  //IMOS20080114

      }

      emissFuncs.resize(0);

   }


   //Do bremsstrahlung
   if (galdef.bremss) 
   {
      emissFuncs.emplace_back(new GasEmissFunction("BREMSS", "HI", *this));
      emissFuncs.emplace_back(new GasEmissFunction("BREMSS", "H2", *this));
      emissFuncs.emplace_back(new GasEmissFunction("BREMSS", "HII", *this));

      if (3 == galdef.skymap_format) 
      { // HEALPix format

         galaxy.bremss_hp_skymap = galaxy.createFullSky(galdef.healpix_order);
         galaxy.bremss_ionized_hp_skymap = galaxy.createFullSky(galdef.healpix_order);

         if (2 == galdef.gamma_rays) 
         { // split according to rings

            galaxy.bremss_H2R_hp_skymap.reserve(galaxy.n_Ring);
            galaxy.bremss_HIR_hp_skymap.reserve(galaxy.n_Ring);
            galaxy.bremss_HIIR_hp_skymap.reserve(galaxy.n_Ring);

            //Do sparse sky for inner Galaxy, outer Galaxy, and all H2R maps
            for (int iR = 0; iR < galaxy.n_Ring; ++iR) 
            { 

               galaxy.bremss_H2R_hp_skymap.emplace_back(galaxy.createSparseSky(galdef.healpix_order));

               if (galaxy.R_bins[iR] < Rsun && galaxy.R_bins[iR+1] > Rsun) {
                  galaxy.bremss_HIR_hp_skymap.emplace_back(galaxy.createFullSky(galdef.healpix_order));
                  galaxy.bremss_HIIR_hp_skymap.emplace_back(galaxy.createFullSky(galdef.healpix_order));
               } else {
                  galaxy.bremss_HIR_hp_skymap.emplace_back(galaxy.createSparseSky(galdef.healpix_order));
                  galaxy.bremss_HIIR_hp_skymap.emplace_back(galaxy.createSparseSky(galdef.healpix_order));
               }

            }

         }


         INFO("Starting bremss integration");

#pragma omp parallel for schedule(dynamic) default(shared)
         for ( int ii = 0; ii < galaxy.bremss_hp_skymap->Npix(); ++ii ) {

           if ( ii % (galaxy.bremss_hp_skymap->Npix()/100) == 0) {
              std::ostringstream buf;
              buf << "Beginning pixel " << ii << " (" << (100*ii)/galaxy.bremss_hp_skymap->Npix() << "%)";
              INFO(buf.str());
           }
	   
	   SM::Coordinate co = galaxy.bremss_hp_skymap->GetCoordinate ( ii );
           double l, b;
           co.getCoordinates(l, b, SM::CoordSys::GAL);
           l *= utl::kConvertRadiansToDegrees;
           b *= utl::kConvertRadiansToDegrees;

            vector< vector<double> > gas;
            vector< vector< std::valarray<double> > > emiss;

            calculatePixel(l, b, 90./galaxy.bremss_hp_skymap->Nside(), cameraTheta, cameraR, Rbins, gasInt, emissInt, gasFuncs, emissFuncs, absorption, gas, emiss);

            //Add the pixels to the maps

            //Correct the emissivities with the real gas maps
            for (size_t iR = 0; iR < emiss[0].size(); ++iR) 
            {

               //No need to go further if there is nothing to add
               if ( emiss[0][iR].max() <= 0 && emiss[1][iR].max() <= 0 && emiss[2][iR].max() <= 0 ) continue;

               if ( galaxy.renormGas ) 
               {
                  //Calculate the average X_CO factor to multiply the CO map with.
                  double av_X_CO;
                  if ( gas[2][iR] > 0 ) 
                  {
                     av_X_CO = gas[1][iR]/gas[2][iR];
                  } 
                  else 
                  {
                     const double R_ring = ( galaxy.R_bins[iR] +galaxy.R_bins[galaxy.n_Ring+iR] ) /2.;
                     av_X_CO = fX_CO ( R_ring );
                  }

                  //The weights to correct the emissivities
                  const double w_HI = gas[0][iR] > 0 ? galaxy.hpHIR[iR]->GetValue(ii,0) / gas[0][iR] : 1;
                  const double w_H2 = gas[1][iR] > 0 ? 2*av_X_CO*galaxy.hpCOR[iR]->GetValue(ii,0) / gas[1][iR] : 1;

                  emiss[0][iR] *= w_HI;
                  emiss[1][iR] *= w_H2;
               }

               for (size_t iE=0; iE < emiss[0][iR].size(); ++iE) {
                  galaxy.bremss_hp_skymap->GetReference(ii, iE) += emiss[0][iR][iE] + emiss[1][iR][iE] + emiss[2][iR][iE];
                  galaxy.bremss_ionized_hp_skymap->GetReference(ii, iE) += emiss[2][iR][iE];
               }

               //Split into rings
               if (galdef.gamma_rays==2 && iR < size_t(galaxy.n_Ring)) 
               {
                  if (emiss[0][iR][0] > 0) 
#pragma omp critical (bremss_HIR_hp_skymap)
                     for (size_t iE=0; iE < emiss[0][iR].size(); ++iE) 
                        galaxy.bremss_HIR_hp_skymap[iR]->GetReference(ii,iE) += emiss[0][iR][iE];
                  if (emiss[1][iR][0] > 0) 
#pragma omp critical (bremss_H2R_hp_skymap)
                     for (size_t iE=0; iE < emiss[1][iR].size(); ++iE) 
                        galaxy.bremss_H2R_hp_skymap[iR]->GetReference(ii,iE) += emiss[1][iR][iE];
                  if (emiss[2][iR][0] > 0) 
#pragma omp critical (bremss_HIIR_hp_skymap)
                     for (size_t iE=0; iE < emiss[2][iR].size(); ++iE) 
                        galaxy.bremss_HIIR_hp_skymap[iR]->GetReference(ii,iE) += emiss[2][iR][iE];
               }

            }

         }

         //Store the maps and delete
         TIME_SUBROUTINE(store_bremss_skymap);

         if (2 == galdef.gamma_rays) 
         {                    //AWS20050302

            TIME_SUBROUTINE(store_bremss_H2R_skymap);
            TIME_SUBROUTINE(store_bremss_HIR_skymap);
            TIME_SUBROUTINE(store_bremss_HII_skymap);

         }

         TIME_SUBROUTINE(store_bremss_ionized_skymap);

         galaxy.bremss_hp_skymap.reset(nullptr);
         
         galaxy.bremss_H2R_hp_skymap.clear();
         galaxy.bremss_HIR_hp_skymap.clear();
         galaxy.bremss_HIIR_hp_skymap.clear();

      } 
      else 
      {

         galaxy.bremss_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);
         galaxy.bremss_ionized_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);

         if (2 == galdef.gamma_rays) 
         {

            galaxy.bremss_H2R_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_Ring, galaxy.n_E_gammagrid);

            galaxy.bremss_HIR_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_Ring, galaxy.n_E_gammagrid);

            galaxy.bremss_HIIR_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_Ring, galaxy.n_E_gammagrid);// IMOS20080114

         }

         INFO("Starting bremss integration");
#pragma omp parallel for schedule(dynamic) default(shared)
         for ( int i_long=0; i_long<galaxy.n_long; i_long++ ) 
         {
            if ( i_long % (galaxy.n_long/100) == 0) 
            {
               ostringstream buf;
               buf << "Beginning longitude row " << i_long << " (" << (100*i_long)/galaxy.n_long << "%)";
               INFO(buf.str());
            }
            for ( int i_lat =0; i_lat <galaxy.n_lat; i_lat++ ) 
            {
               const double l=galaxy.long_min + i_long*galaxy.d_long;
               const double b=galaxy. lat_min + i_lat *galaxy.d_lat ;

               vector< vector<double> > gas;
               vector< vector< std::valarray<double> > > emiss;

               calculatePixel(l, b, 360./galaxy.n_long, cameraTheta, cameraR, Rbins, gasInt, emissInt, gasFuncs, emissFuncs, absorption, gas, emiss);

               //Add the pixels to the maps

               //Correct the emissivities with the real gas maps
               for (size_t iR = 0; iR < emiss[0].size(); ++iR) 
               {

                  //No need to go further if there is nothing to add
                  if ( galaxy.renormGas && (gas[0][iR] <= 0 || gas[1][iR] <= 0) ) continue;

                  if ( galaxy.renormGas ) 
                  {
                     //Calculate the average X_CO factor to multiply the CO map with.
                     double av_X_CO;
                     if ( gas[2][iR] > 0 ) 
                     {
                        av_X_CO = gas[1][iR]/gas[2][iR];
                     } 
                     else 
                     {
                        const double R_ring = ( galaxy.R_bins[iR] +galaxy.R_bins[galaxy.n_Ring+iR] ) /2.;
                        av_X_CO = fX_CO ( R_ring );
                     }

                     //The weights to correct the emissivities
                     const double w_HI = galaxy.HIR.d3[i_long][i_lat][iR].s[0] / gas[0][iR];
                     const double w_H2 = 2*av_X_CO*galaxy.COR.d3[i_long][i_lat][iR].s[0] / gas[1][iR];

                     emiss[0][iR] *= w_HI;
                     emiss[1][iR] *= w_H2;
                  }

                  for ( int iE=0; iE < galaxy.n_E_gammagrid; ++iE ) 
                  {

                     galaxy.bremss_skymap.d2[i_long][i_lat].s[iE] += emiss[0][iR][iE] + emiss[1][iR][iE] + emiss[2][iR][iE];
                     galaxy.bremss_ionized_skymap.d2[i_long][i_lat].s[iE] += emiss[2][iR][iE];

                  }

                  //Split into rings
                  if (2 == galdef.gamma_rays && iR < size_t(galaxy.n_Ring)) 
                  {

                     for ( int iE=0; iE < galaxy.n_E_gammagrid; ++iE ) 
                     {

                        galaxy.bremss_HIR_skymap.d3[i_long][i_lat][iR].s[iE] += emiss[0][iR][iE];
                        galaxy.bremss_H2R_skymap.d3[i_long][i_lat][iR].s[iE] += emiss[1][iR][iE];
                        galaxy.bremss_HIIR_skymap.d3[i_long][i_lat][iR].s[iE] += emiss[2][iR][iE];

                     }
                  }

               }

            }
         }

         //Store the maps and delete
         TIME_SUBROUTINE(store_bremss_skymap);

         if (2 == galdef.gamma_rays) 
         {                    //AWS20050302

            TIME_SUBROUTINE(store_bremss_H2R_skymap);
            TIME_SUBROUTINE(store_bremss_HIR_skymap);
            TIME_SUBROUTINE(store_bremss_HII_skymap);

         }

         TIME_SUBROUTINE(store_bremss_ionized_skymap);

         galaxy.bremss_skymap.delete_array();  //Gulli20070810

         galaxy.bremss_H2R_skymap.delete_array();  //Gulli20070810
         galaxy.bremss_HIR_skymap.delete_array();  //Gulli20070810
         galaxy.bremss_HIIR_skymap.delete_array();  //IMOS20080114

      }

      emissFuncs.resize(0);

   }

   //Delete the gas functions, they are not needed anymore
   gasFuncs.resize(0);
   galaxy.COR.delete_array();
   galaxy.HIR.delete_array();


   //Do IC integration if needed
   auto anisoValid = false;

   if (galdef.IC_isotropic || galdef.IC_anisotropic) 
   {

      std::vector<double> z(&galaxy.z[0], &galaxy.z[0]+galaxy.n_zgrid);

      if (2 == galaxy.n_spatial_dimensions) 
      {
         std::vector<double> r(&galaxy.r[0], &galaxy.r[0]+galaxy.n_rgrid);

         for (int i_comp = 0; i_comp < galaxy.n_ISRF_components; ++i_comp) 
            emissFuncs.emplace_back(new DistributionFunction(galaxy.IC_iso_emiss[i_comp], z, r));

      } 
      else if (3 == galaxy.n_spatial_dimensions) 
      {
         std::vector<double> x(&galaxy.x[0], &galaxy.x[0]+galaxy.n_xgrid);
         std::vector<double> y(&galaxy.y[0], &galaxy.y[0]+galaxy.n_ygrid);

         for (int i_comp = 0; i_comp < galaxy.n_ISRF_components; ++i_comp) 
            emissFuncs.emplace_back(new DistributionFunction(galaxy.IC_iso_emiss[i_comp], z, x, y));

      }

      //This has to be done in a separate step, creating an aniso to iso ratio maps.
      if (galdef.IC_anisotropic) {
	
	const auto nComps = galaxy.n_ISRF_components;
	
	//         bool haveOMP = false;
	
	//#ifdef _OPENMP
	//haveOMP = true;
	//#endif
	
	if (galdef.ISRF_filetype < 3 || nComps < 3) {
	  
	  std::ostringstream buf;
	  buf << "Anisotropic IC calculation only available for version 3 ISRF files and later with all 3 ISRF components (optical, infrared, CMB). Only isotropic IC will be calculated.";
	  INFO(buf.str());
	  galdef.IC_anisotropic = 0; //Need to set this to false because we initialize the maps in this construct
	  
	} else {
            
	  std::vector< std::unique_ptr< SM::LOSfunction<std::valarray<double> > > > anisoFuncs(1);
	  try {
	    anisoFuncs[0].reset(new AnisoICFunction(*this));
	    anisoValid = true;
	  } catch (std::invalid_argument) {
	    anisoValid = false;
	    galdef.IC_anisotropic = 0; //Need to set this to false because we initialize the maps in this routine
	  }
	  
	  if (anisoValid) {
	    
	    const auto internalHealpixOrder = galdef.anisoHealpixOrder; //To save time we calculate for a coarse grid and interpolate
	    
            auto anisoSkymapOpt = galaxy.createFullSky(internalHealpixOrder);
            auto anisoSkymapIR = galaxy.createFullSky(internalHealpixOrder);
            auto anisoSkymapCMB = galaxy.createFullSky(internalHealpixOrder);
            auto isoSkymapOpt = galaxy.createFullSky(internalHealpixOrder);
            auto isoSkymapIR = galaxy.createFullSky(internalHealpixOrder);
            auto isoSkymapCMB = galaxy.createFullSky(internalHealpixOrder);
	    
	    std::vector< std::unique_ptr< SM::BaseSky<double> >  > anisoSkymapOptRings, 
	      anisoSkymapIRRings, 
	      anisoSkymapCMBRings, 
	      isoSkymapOptRings, 
	      isoSkymapIRRings, 
	      isoSkymapCMBRings;
	    
	    if (2 == galdef.gamma_rays) {
	      
	      anisoSkymapOptRings.reserve(galaxy.n_Ring);
	      anisoSkymapIRRings.reserve(galaxy.n_Ring);
	      anisoSkymapCMBRings.reserve(galaxy.n_Ring);
	      isoSkymapOptRings.reserve(galaxy.n_Ring);
	      isoSkymapIRRings.reserve(galaxy.n_Ring);
	      isoSkymapCMBRings.reserve(galaxy.n_Ring);
	      
	      for (int iR(0); iR < galaxy.n_Ring; ++iR) {

                 //Do Sparse for inner galaxy only
                 if (galaxy.R_bins[iR+1] < Rsun) {

                    anisoSkymapOptRings.emplace_back(galaxy.createSparseSky(internalHealpixOrder));
                    anisoSkymapIRRings.emplace_back(galaxy.createSparseSky(internalHealpixOrder));
                    anisoSkymapCMBRings.emplace_back(galaxy.createSparseSky(internalHealpixOrder));
                    isoSkymapOptRings.emplace_back(galaxy.createSparseSky(internalHealpixOrder));
                    isoSkymapIRRings.emplace_back(galaxy.createSparseSky(internalHealpixOrder));
                    isoSkymapCMBRings.emplace_back(galaxy.createSparseSky(internalHealpixOrder));

                 } else {
		
                    anisoSkymapOptRings.emplace_back(galaxy.createFullSky(internalHealpixOrder));
                    anisoSkymapIRRings.emplace_back(galaxy.createFullSky(internalHealpixOrder));
                    anisoSkymapCMBRings.emplace_back(galaxy.createFullSky(internalHealpixOrder));
                    isoSkymapOptRings.emplace_back(galaxy.createFullSky(internalHealpixOrder));
                    isoSkymapIRRings.emplace_back(galaxy.createFullSky(internalHealpixOrder));
                    isoSkymapCMBRings.emplace_back(galaxy.createFullSky(internalHealpixOrder));

                 }

                 //Default value of 1 for the anisoSkymaps that are used to store the ratios later
                 anisoSkymapOptRings.back()->ApplyFunction([](double &a) { a = 1.0; });
                 anisoSkymapIRRings.back()->ApplyFunction([](double &a) { a = 1.0; });
                 anisoSkymapCMBRings.back()->ApplyFunction([](double &a) { a = 1.0; });
		
	      }
	      
	    }
	    
	    //Have a different integration step size
	    const double stepSize = 1e-3;
	    
            std::unique_ptr< SM::LOSintegrator< std::valarray<double> > > anisoInt;
	    if (2 == galdef.n_spatial_dimensions) {
	      anisoInt.reset( new SM::LOSintegrator<std::valarray<double> >(galaxy.r_max, galaxy.z_min, galaxy.z_max, Rbins, galdef.fCameraLocation, stepSize, 1e-2, 1e-4));
	    } else if (3 == galdef.n_spatial_dimensions) {
	      anisoInt.reset( new SM::LOSintegrator<std::valarray<double> >(galaxy.x_min, galaxy.x_max, galaxy.y_min, galaxy.y_max, galaxy.z_min, galaxy.z_max, Rbins, galdef.fCameraLocation, stepSize, 1e-2, 1e-4));
	    }

	    //Loop over each pixel and integrate
#pragma omp parallel for schedule(dynamic) default(shared)
	    for (auto iPix = 0; iPix < anisoSkymapOpt->Npix(); ++iPix) {

              std::ostringstream buf;
              buf << "Beginning pixel " << iPix << " (" << (100*iPix)/anisoSkymapOpt->Npix() << "%)";
              INFO(buf.str());

	      auto coord = anisoSkymapOpt->GetCoordinate(iPix);
	      
              double l, b;
              coord.getCoordinates(l, b, SM::CoordSys::GAL);
              l *= utl::kConvertRadiansToDegrees;
              b *= utl::kConvertRadiansToDegrees;
	      
	      const std::vector<std::vector<std::valarray<double> > > emiss = anisoInt->integrate(l, b, anisoFuncs);
	      
	      for ( size_t iR=0; iR < emiss[0].size(); ++iR) {
		for ( int iE=0; iE < galaxy.n_E_gammagrid; ++iE ) {
		  
		  isoSkymapOpt->GetReference(iPix,iE) += emiss[0][iR][0*galaxy.n_E_gammagrid + iE];
		  isoSkymapIR->GetReference(iPix,iE) += emiss[0][iR][1*galaxy.n_E_gammagrid + iE];
		  isoSkymapCMB->GetReference(iPix,iE) += emiss[0][iR][2*galaxy.n_E_gammagrid + iE];
		  anisoSkymapOpt->GetReference(iPix,iE) += emiss[0][iR][3*galaxy.n_E_gammagrid + iE];
		  anisoSkymapIR->GetReference(iPix,iE) += emiss[0][iR][4*galaxy.n_E_gammagrid + iE];
		  anisoSkymapCMB->GetReference(iPix,iE) += emiss[0][iR][5*galaxy.n_E_gammagrid + iE];

                }
		  
                if (2 == galdef.gamma_rays && iR < size_t(galaxy.n_Ring)) {
#pragma omp critical (anisoIC)
                   for ( int iE=0; iE < galaxy.n_E_gammagrid; ++iE ) {
                      isoSkymapOptRings[iR]->GetReference(iPix,iE) = emiss[0][iR][0*galaxy.n_E_gammagrid + iE];
                      isoSkymapIRRings[iR]->GetReference(iPix,iE) = emiss[0][iR][1*galaxy.n_E_gammagrid + iE];
                      isoSkymapCMBRings[iR]->GetReference(iPix,iE) = emiss[0][iR][2*galaxy.n_E_gammagrid + iE];
                      anisoSkymapOptRings[iR]->GetReference(iPix,iE) = emiss[0][iR][3*galaxy.n_E_gammagrid + iE];
                      anisoSkymapIRRings[iR]->GetReference(iPix,iE) = emiss[0][iR][4*galaxy.n_E_gammagrid + iE];
                      anisoSkymapCMBRings[iR]->GetReference(iPix,iE) = emiss[0][iR][5*galaxy.n_E_gammagrid + iE];
                   }
		  
                }
	      }
	    }

	    // Reuse the skymaps for memory efficiency, they are not used again anyway
            SM::ApplyFunction<double,double>(*anisoSkymapOpt, *isoSkymapOpt, [](double &a, const double &b) { a = b > 0 ? a/b : 1; });
            SM::ApplyFunction<double,double>(*anisoSkymapIR, *isoSkymapIR, [](double &a, const double &b) { a = b > 0 ? a/b : 1; });
            SM::ApplyFunction<double,double>(*anisoSkymapCMB, *isoSkymapCMB, [](double &a, const double &b) { a = b > 0 ? a/b : 1; });

            if (2 == galdef.gamma_rays ) {

               for ( int iR=0; iR < galaxy.n_Ring; ++iR) {
		    
                  SM::ApplyFunction<double,double>(*anisoSkymapOptRings[iR], *isoSkymapOptRings[iR], [](double &a, const double &b) { a = b > 0 ? a/b : 1; });
                  SM::ApplyFunction<double,double>(*anisoSkymapIRRings[iR], *isoSkymapIRRings[iR], [](double &a, const double &b) { a = b > 0 ? a/b : 1; });
                  SM::ApplyFunction<double,double>(*anisoSkymapCMBRings[iR], *isoSkymapCMBRings[iR], [](double &a, const double &b) { a = b > 0 ? a/b : 1; });

               }
            }
		
	    
            /*
             * This is only valid in very special cases and should be set up as a unit test at some point
	    if (galdef.verbose == -1843) 
	      {
		for (size_t iSpec = 0; iSpec < gammaE.size(); ++iSpec) 
                  {
		    for (int iPix = 0; iPix < anisoSkymapOpt.Npix(); ++iPix) 
		      {
			
                        //Check the pixels where the anisoSkymapCMB is not symmetric around the GC and GP
                        //Only need to check for one quadrant
                        SM::Coordinate co = anisoSkymapOpt.pix2coord(iPix);
                        if (co.l() < 180 && co.l() >= 0 && co.b() > 0 ) 
			  {
			    SM::Coordinate co1(co.l(), -co.b()), co2(-co.l(), co.b()), co3(-co.l(), -co.b());
			    if (fabs(anisoSkymapCMB[co][iSpec] - anisoSkymapCMB[co1][iSpec])/anisoSkymapCMB[co][iSpec] > 1e-10) 
			      {
				std::cout<<"AnisoCMB not symmetric, "<<co<<", "<<co1<<", "<<iSpec<<", "<<1-anisoSkymapCMB[co1][iSpec]/anisoSkymapCMB[co][iSpec]<<std::endl;
			      }
			    if (fabs(anisoSkymapCMB[co][iSpec] - anisoSkymapCMB[co2][iSpec])/anisoSkymapCMB[co][iSpec] > 1e-10) 
			      {
				std::cout<<"AnisoCMB not symmetric, "<<co<<", "<<co2<<", "<<iSpec<<", "<<1-anisoSkymapCMB[co2][iSpec]/anisoSkymapCMB[co][iSpec]<<std::endl;
			      }
			    if (fabs(anisoSkymapCMB[co][iSpec] - anisoSkymapCMB[co3][iSpec])/anisoSkymapCMB[co][iSpec] > 1e-10) 
			      {
				std::cout<<"AnisoCMB not symmetric, "<<co<<", "<<co3<<", "<<iSpec<<", "<<1-anisoSkymapCMB[co3][iSpec]/anisoSkymapCMB[co][iSpec]<<std::endl;
			      }
			  }
		      }
                  }
	      }
              */
	    
	    //Use the aniso maps to store the ratio.
	    if (3 == galdef.skymap_format) {
	      
	      galaxy.IC_aniso_hp_skymap.resize(galaxy.n_ISRF_components);
	      
              // Create the maps
              for (int i = 0; i < galaxy.n_ISRF_components; ++i)
                 galaxy.IC_aniso_hp_skymap[i] = galaxy.createFullSky(galdef.healpix_order);

	      // interpolate to the desired healpix order
              anisoSkymapOpt->Interpolate(*galaxy.IC_aniso_hp_skymap[0], false);
              anisoSkymapIR->Interpolate(*galaxy.IC_aniso_hp_skymap[1], false);
              anisoSkymapCMB->Interpolate(*galaxy.IC_aniso_hp_skymap[2], false);
	      
	      if (2 == galdef.gamma_rays) { // split according to rings
		
		galaxy.IC_aniso_rings_hp_skymap.resize(galaxy.n_ISRF_components);
		
		for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
		  
		  galaxy.IC_aniso_rings_hp_skymap[icomp].resize(galaxy.n_Ring);
		  
		}
		
		for ( int iR(0); iR < galaxy.n_Ring; ++iR) {
		  
                   //Create the maps
                   if (galaxy.R_bins[iR+1] < Rsun) {
                      for (int i = 0; i < galaxy.n_ISRF_components; ++i)
                         galaxy.IC_aniso_rings_hp_skymap[i][iR] = galaxy.createSparseSky(galdef.healpix_order);
                   } else {
                      for (int i = 0; i < galaxy.n_ISRF_components; ++i)
                         galaxy.IC_aniso_rings_hp_skymap[i][iR] = galaxy.createFullSky(galdef.healpix_order);
                   }

                   // interpolate to the desired healpix order
                   anisoSkymapOptRings[iR]->Interpolate(*galaxy.IC_aniso_rings_hp_skymap[0][iR], false);
                   anisoSkymapIRRings[iR]->Interpolate(*galaxy.IC_aniso_rings_hp_skymap[1][iR], false);
                   anisoSkymapCMBRings[iR]->Interpolate(*galaxy.IC_aniso_rings_hp_skymap[2][iR], false);

                }
		
	      }
	      	      
	    } else {
	      
	      galaxy.IC_aniso_skymap = new Distribution[galaxy.n_ISRF_components];
	      
	      for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {

		galaxy.IC_aniso_skymap[icomp].init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);
		
	      }
	      
	      if (2 == galdef.gamma_rays) { // split according to rings
		
		galaxy.IC_aniso_rings_skymap = new Distribution[galaxy.n_ISRF_components];
		
		for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
		  
		  galaxy.IC_aniso_rings_skymap[icomp].init(galaxy.n_long, galaxy.n_lat, galaxy.n_Ring, galaxy.n_E_gammagrid);
		  
		}
		
	      }
	      
	      for ( int i_long=0; i_long<galaxy.n_long; i_long++ ) {
		for ( int i_lat =0; i_lat <galaxy.n_lat; i_lat++ ) {
                  const double l=(galaxy.long_min + i_long*galaxy.d_long)*utl::kConvertDegreesToRadians;
                  const double b=(galaxy. lat_min + i_lat *galaxy.d_lat )*utl::kConvertDegreesToRadians;
		  const SM::Coordinate co(l,b, SM::CoordSys::GAL);

                  //Get the interpolated spectra
                  auto anisoRatioOpt = anisoSkymapOpt->GetSkyInterpolatedSpectra(co);
                  auto anisoRatioIR = anisoSkymapIR->GetSkyInterpolatedSpectra(co);
                  auto anisoRatioCMB = anisoSkymapCMB->GetSkyInterpolatedSpectra(co);
		  
		  for ( int iE=0; iE < galaxy.n_E_gammagrid; ++iE ) {
		    galaxy.IC_aniso_skymap[0].d2[i_long][i_lat].s[iE] = anisoRatioOpt[iE]; 
		    galaxy.IC_aniso_skymap[1].d2[i_long][i_lat].s[iE] = anisoRatioIR[iE]; 
		    galaxy.IC_aniso_skymap[2].d2[i_long][i_lat].s[iE] = anisoRatioCMB[iE]; 
		  }
		}
	      }
	      
	      if (2 == galdef.gamma_rays) { // split according to rings

		for ( int iR(0); iR < galaxy.n_Ring; ++iR) {
		  
                  //This could be done more efficiently, but this should be removed anyway
		  for ( int i_long=0; i_long<galaxy.n_long; i_long++ ) {
		    for ( int i_lat =0; i_lat <galaxy.n_lat; i_lat++ ) {
		      
		      const double l=(galaxy.long_min + i_long*galaxy.d_long)*utl::kConvertDegreesToRadians;
		      const double b=(galaxy. lat_min + i_lat *galaxy.d_lat )*utl::kConvertDegreesToRadians;
		      const SM::Coordinate co(l,b, SM::CoordSys::GAL);
		      
                      //Get the interpolated spectra
                      auto anisoRatioOpt = anisoSkymapOpt->GetSkyInterpolatedSpectra(co);
                      auto anisoRatioIR = anisoSkymapIR->GetSkyInterpolatedSpectra(co);
                      auto anisoRatioCMB = anisoSkymapCMB->GetSkyInterpolatedSpectra(co);

		      for ( int iE=0; iE < galaxy.n_E_gammagrid; ++iE ) {
			galaxy.IC_aniso_rings_skymap[0].d3[i_long][i_lat][iR].s[iE] = anisoRatioOpt[iE]; 
			galaxy.IC_aniso_rings_skymap[1].d3[i_long][i_lat][iR].s[iE] = anisoRatioIR[iE]; 
			galaxy.IC_aniso_rings_skymap[2].d3[i_long][i_lat][iR].s[iE] = anisoRatioCMB[iE]; 
		      }
		    }
		  }
		  
		}
		
	      }
	      
	    }
	    
	  }
	}
      }
      
      if (3 == galdef.skymap_format) {
	
	galaxy.IC_iso_hp_skymap.resize(galaxy.n_ISRF_components);
	
	for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {

           galaxy.IC_iso_hp_skymap[icomp] = galaxy.createFullSky(galdef.healpix_order);
	  
	}
	
	if (2 == galdef.gamma_rays) { // split according to rings
	  
	  galaxy.IC_iso_rings_hp_skymap.resize(galaxy.n_ISRF_components);
	  
	  for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
	      
	    galaxy.IC_iso_rings_hp_skymap[icomp].resize(galaxy.n_Ring);
	    
	    for (int iR = 0; iR < galaxy.n_Ring; ++iR) { 
	      
               //Do Sparse for inner galaxy only
               if (galaxy.R_bins[iR+1] < Rsun) {

                  galaxy.IC_iso_rings_hp_skymap[icomp][iR] = galaxy.createSparseSky(galdef.healpix_order);

               } else {

                  galaxy.IC_iso_rings_hp_skymap[icomp][iR] = galaxy.createFullSky(galdef.healpix_order);

               }

	    }

	  }
	  
	}

	INFO("Starting IC integration");

#pragma omp parallel for schedule(dynamic) default(shared)
	for ( int ii = 0; ii < galaxy.IC_iso_hp_skymap[0]->Npix(); ++ii ) 
         {
	   if ( ii % (galaxy.IC_iso_hp_skymap[0]->Npix()/100) == 0) 
           {
	      ostringstream buf;
	      buf << "Beginning pixel " << ii << " (" << (100*ii)/galaxy.IC_iso_hp_skymap[0]->Npix() << "%)";
	      INFO(buf.str());
           }

           auto co = galaxy.IC_iso_hp_skymap[0]->GetCoordinate( ii );
           double l, b;
           co.getCoordinates(l, b, SM::CoordSys::GAL);
           l *= utl::kConvertRadiansToDegrees;
           b *= utl::kConvertRadiansToDegrees;
	   
	   vector< vector<double> > gas;
	   vector< vector< std::valarray<double> > > emiss;
	   
	   calculatePixel(l, b, 90./galaxy.IC_iso_hp_skymap[0]->Nside(), cameraTheta, cameraR, Rbins, gasInt, emissInt, gasFuncs, emissFuncs, absorption, gas, emiss);
	   
	   //Add the pixels to the maps
	   for (int iC = 0; iC < galaxy.n_ISRF_components; ++iC) 
	     {
	       
               for (size_t iR = 0; iR < emiss[iC].size(); ++iR) 
		 {
		   
                    for (size_t iE = 0; iE < emiss[iC][iR].size(); ++iE)
                       galaxy.IC_iso_hp_skymap[iC]->GetReference(ii, iE) += emiss[iC][iR][iE];
		   
                  //Split into rings
		   if (galdef.gamma_rays==2 && iR < size_t(galaxy.n_Ring) && emiss[iC][iR][0] > 0) 
                  {
#pragma omp critical (IC_iso_rings_hp_skymap)
                     for (size_t iE = 0; iE < emiss[iC][iR].size(); ++iE)
                        galaxy.IC_iso_rings_hp_skymap[iC][iR]->GetReference(ii,iE) += emiss[iC][iR][iE];

                  }

               }

            }

         }

        //Calculate the aniso maps
        if (anisoValid) {

	   for (int iC = 0; iC < galaxy.n_ISRF_components; ++iC) 
              SM::ApplyFunction<double,double>(*galaxy.IC_aniso_hp_skymap[iC], *galaxy.IC_iso_hp_skymap[iC], [](double &a, const double &b) { a *= b; } );

           if (galdef.gamma_rays==2) 
              for (int iC = 0; iC < galaxy.n_ISRF_components; ++iC) 
                 for (int iR = 0; iR < galaxy.n_Ring; ++iR) 
                    SM::ApplyFunction<double,double>(*galaxy.IC_aniso_rings_hp_skymap[iC][iR], *galaxy.IC_iso_rings_hp_skymap[iC][iR], [](double &a, const double &b) { a *= b; } );
        }
      } 
      else 
      {

         galaxy.IC_iso_skymap = new Distribution[galaxy.n_ISRF_components];

         for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) 
         {

            galaxy.IC_iso_skymap[icomp].init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);

         }

         if (2 == galdef.gamma_rays) 
         { // split according to rings

            galaxy.IC_iso_rings_skymap = new Distribution[galaxy.n_ISRF_components];

            for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) 
            {

               galaxy.IC_iso_rings_skymap[icomp].init(galaxy.n_long, galaxy.n_lat, galaxy.n_Ring, galaxy.n_E_gammagrid);

            }

         }

         INFO("Starting IC integration");
#pragma omp parallel for schedule(dynamic) default(shared)
         for ( int i_long=0; i_long<galaxy.n_long; i_long++ ) 
         {
            if ( i_long % (galaxy.n_long/100) == 0) 
            {
               ostringstream buf;
               buf << "Beginning longitude row " << i_long << " (" << (100*i_long)/galaxy.n_long << "%)";
               INFO(buf.str());
            }
            for ( int i_lat =0; i_lat <galaxy.n_lat; i_lat++ ) 
            {
               const double l=galaxy.long_min + i_long*galaxy.d_long;
               const double b=galaxy. lat_min + i_lat *galaxy.d_lat ;

               vector< vector<double> > gas;
               vector< vector< std::valarray<double> > > emiss;

               calculatePixel(l, b, 360./galaxy.n_long, cameraTheta, cameraR, Rbins, gasInt, emissInt, gasFuncs, emissFuncs, absorption, gas, emiss);

               //Add the pixels to the maps
               for (int iC = 0; iC < galaxy.n_ISRF_components; ++iC) 
               {

                  for (size_t iR = 0; iR < emiss[0].size(); ++iR) 
                  {

                     for ( int iE=0; iE < galaxy.n_E_gammagrid; ++iE ) 
                        galaxy.IC_iso_skymap[iC].d2[i_long][i_lat].s[iE] += emiss[iC][iR][iE];

                     //Split into rings
                     if (galdef.gamma_rays==2 && iR < size_t(galaxy.n_Ring)) 
                     {

                        for ( int iE=0; iE < galaxy.n_E_gammagrid; ++iE ) 
                           galaxy.IC_iso_rings_skymap[iC].d3[i_long][i_lat][iR].s[iE] += emiss[iC][iR][iE];

                     }

                  }

               }

               if (anisoValid) 
               {
                  for ( int iE=0; iE < galaxy.n_E_gammagrid; ++iE ) 
                  {
                     galaxy.IC_aniso_skymap[0].d2[i_long][i_lat].s[iE] *= galaxy.IC_iso_skymap[0].d2[i_long][i_lat].s[iE]; 
                     galaxy.IC_aniso_skymap[1].d2[i_long][i_lat].s[iE] *= galaxy.IC_iso_skymap[1].d2[i_long][i_lat].s[iE]; 
                     galaxy.IC_aniso_skymap[2].d2[i_long][i_lat].s[iE] *= galaxy.IC_iso_skymap[2].d2[i_long][i_lat].s[iE];
                  }

                  if (galdef.gamma_rays==2) 
                  {

                     for (int iR = 0; iR < galaxy.n_Ring; ++iR) 
                     {

                        for ( int iE=0; iE < galaxy.n_E_gammagrid; ++iE ) 
                        {
                           galaxy.IC_aniso_rings_skymap[0].d3[i_long][i_lat][iR].s[iE] *= galaxy.IC_iso_rings_skymap[0].d3[i_long][i_lat][iR].s[iE]; 
                           galaxy.IC_aniso_rings_skymap[1].d3[i_long][i_lat][iR].s[iE] *= galaxy.IC_iso_rings_skymap[1].d3[i_long][i_lat][iR].s[iE]; 
                           galaxy.IC_aniso_rings_skymap[2].d3[i_long][i_lat][iR].s[iE] *= galaxy.IC_iso_rings_skymap[2].d3[i_long][i_lat][iR].s[iE]; 
                        }

                     }

                  }
               }

            }
         }
      }

      //Store the maps and delete allocated memory
      string icType = "isotropic";

      TIME_SUBROUTINE(store_IC_skymap_comp, icType);
      TIME_SUBROUTINE(store_IC_skymap, icType); //AWS20090415

      if (galdef.IC_anisotropic) 
      {

         icType = "anisotropic";

         TIME_SUBROUTINE(store_IC_skymap_comp, icType);
         TIME_SUBROUTINE(store_IC_skymap, icType); 

      }

      if (3 == galdef.skymap_format) 
      {

         galaxy.IC_iso_hp_skymap.clear();
         galaxy.IC_iso_rings_hp_skymap.clear();
         galaxy.IC_aniso_hp_skymap.clear();
         galaxy.IC_aniso_rings_hp_skymap.clear();

      } 
      else 
      {

         for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) 
         {

            galaxy.IC_iso_skymap[icomp].delete_array();

         }

         if (galdef.IC_anisotropic) 
         {

            for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) 
            {

               galaxy.IC_aniso_skymap[icomp].delete_array();

            }

         }

      }

      emissFuncs.resize(0);
   }


   if (galdef.DM_gammas) 
   {                           // IMOS20050912

      std::vector<double> z(&galaxy.z[0], &galaxy.z[0]+galaxy.n_zgrid);

      if (2 == galaxy.n_spatial_dimensions) 
      {
         std::vector<double> r(&galaxy.r[0], &galaxy.r[0]+galaxy.n_rgrid);
         emissFuncs.emplace_back(new DistributionFunction(galaxy.DM_emiss, z, r));

      } 
      else if (3 == galaxy.n_spatial_dimensions) 
      {
         std::vector<double> x(&galaxy.x[0], &galaxy.x[0]+galaxy.n_xgrid);
         std::vector<double> y(&galaxy.y[0], &galaxy.y[0]+galaxy.n_ygrid);
         emissFuncs.emplace_back(new DistributionFunction(galaxy.DM_emiss, z, x, y));

      }

      if (3 == galdef.skymap_format) 
      {

         galaxy.DM_hp_skymap = galaxy.createFullSky(galdef.healpix_order);

         INFO("Starting DM integration");

#pragma omp parallel for schedule(dynamic) default(shared)
         for ( int ii = 0; ii < galaxy.DM_hp_skymap->Npix(); ++ii ) 
         {
            if ( ii % (galaxy.DM_hp_skymap->Npix()/100) == 0) 
            {
               ostringstream buf;
               buf << "Beginning pixel " << ii << " (" << (100*ii)/galaxy.DM_hp_skymap->Npix() << "%)";
               INFO(buf.str());
            }
            auto co = galaxy.DM_hp_skymap->GetCoordinate( ii );
            double l, b;
            co.getCoordinates(l, b, SM::CoordSys::GAL);
            l *= utl::kConvertRadiansToDegrees;
            b *= utl::kConvertRadiansToDegrees;

            vector< vector<double> > gas;
            vector< vector< std::valarray<double> > > emiss;

            calculatePixel(l, b, 90./galaxy.DM_hp_skymap->Nside(), cameraTheta, cameraR, Rbins, gasInt, emissInt, gasFuncs, emissFuncs, absorption, gas, emiss);

            //Add the pixels to the maps
            for (size_t iR = 0; iR < emiss[0].size(); ++iR) 
               for (size_t iE = 0; iE < emiss[0][iR].size(); ++iE)
                  galaxy.DM_hp_skymap->GetReference(ii,iE) += emiss[0][iR][iE];
         }

      } 
      else 
      {

         galaxy.DM_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid); // DM: IMOS20050912 IMOS20060420

         INFO("Starting DM integration");
#pragma omp parallel for schedule(dynamic) default(shared)
         for ( int i_long=0; i_long<galaxy.n_long; i_long++ ) 
         {
            if ( i_long % (galaxy.n_long/100) == 0) 
            {
               ostringstream buf;
               buf << "Beginning longitude row " << i_long << " (" << (100*i_long)/galaxy.n_long << "%)";
               INFO(buf.str());
            }
            for ( int i_lat =0; i_lat <galaxy.n_lat; i_lat++ ) 
            {
               const double l=galaxy.long_min + i_long*galaxy.d_long;
               const double b=galaxy. lat_min + i_lat *galaxy.d_lat ;

               vector< vector<double> > gas;
               vector< vector< std::valarray<double> > > emiss;

               calculatePixel(l, b, 360./galaxy.n_long, cameraTheta, cameraR, Rbins, gasInt, emissInt, gasFuncs, emissFuncs, absorption, gas, emiss);

               //Add the pixels to the maps
               for (size_t iR = 0; iR < emiss[0].size(); ++iR) 
                  for ( int iE=0; iE < galaxy.n_E_gammagrid; ++iE ) 
                     galaxy.DM_skymap.d2[i_long][i_lat].s[iE] += emiss[0][iR][iE];
            }
         }
      }
      TIME_SUBROUTINE(store_DM_skymap);

      if (3 == galdef.skymap_format) 
      {

         galaxy.DM_hp_skymap.reset(nullptr);

      } 
      else 
      {

         galaxy.DM_skymap.delete_array();//Gulli20070810

      }

      emissFuncs.resize(0);
   }

   absorption.reset( nullptr );

   //Now we can do synchrotron.  It has different absorption and has to be done separately
   if (galdef.synchrotron) 
   {
      std::vector<double> z(&galaxy.z[0], &galaxy.z[0]+galaxy.n_zgrid);

      if (2 == galaxy.n_spatial_dimensions) 
      {
         std::vector<double> r(&galaxy.r[0], &galaxy.r[0]+galaxy.n_rgrid);
         emissFuncs.emplace_back(new DistributionFunction(galaxy.synchrotron_emiss, z, r));
         emissFuncs.emplace_back(new DistributionFunction(galaxy.synchrotron_Q_emiss, z, r));
         emissFuncs.emplace_back(new DistributionFunction(galaxy.synchrotron_U_emiss, z, r));

      } 
      else if (3 == galaxy.n_spatial_dimensions) 
      {
         std::vector<double> x(&galaxy.x[0], &galaxy.x[0]+galaxy.n_xgrid);
         std::vector<double> y(&galaxy.y[0], &galaxy.y[0]+galaxy.n_ygrid);
         emissFuncs.emplace_back(new DistributionFunction(galaxy.synchrotron_emiss, z, x, y));
         emissFuncs.emplace_back(new DistributionFunction(galaxy.synchrotron_Q_emiss, z, x, y));
         emissFuncs.emplace_back(new DistributionFunction(galaxy.synchrotron_U_emiss, z, x, y));

      } 
      else 
      {

	// Error handling code?

      }

      if ( galdef.free_free_absorption >= 1 ) 
      {
         emissFuncs.emplace_back(new synchro::EmissFreeFree(galaxy.nu_synch, galdef.HII_clumping_factor, galdef.HII_Te));
         absorption.reset(new synchro::KappaFreeFree(galaxy.nu_synch, galdef.HII_clumping_factor, galdef.HII_Te));
      }

      if (3 == galdef.skymap_format) 
      { // HEALPix format

         //Cannot use the createFullSky, that assumes gamma-rays
         galaxy.synchrotron_hp_skymap.reset(new SM::FullSky<double> ( 
                  SM::SpectralBinning(std::vector<double>(std::begin(galaxy.nu_synch), std::end(galaxy.nu_synch))), 
                  galdef.healpix_order, RING, SM::CoordSys::GAL, 0.0));
         galaxy.synchrotron_Q_hp_skymap = galaxy.synchrotron_hp_skymap->clone();
         galaxy.synchrotron_U_hp_skymap = galaxy.synchrotron_hp_skymap->clone();
         galaxy.synchrotron_P_hp_skymap = galaxy.synchrotron_hp_skymap->clone();
         galaxy.synchrotron_polang_hp_skymap = galaxy.synchrotron_hp_skymap->clone();
         galaxy.synchrotron_polfra_hp_skymap = galaxy.synchrotron_hp_skymap->clone();
         galaxy.free_free_hp_skymap = galaxy.synchrotron_hp_skymap->clone();

      } 
      else 
      {

         galaxy.synchrotron_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_nu_synchgrid);
         galaxy.synchrotron_Q_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_nu_synchgrid);                    //AWS20100708
         galaxy.synchrotron_U_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_nu_synchgrid);                    //AWS20100708
         galaxy.synchrotron_P_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_nu_synchgrid);                    //AWS20100708
         galaxy.synchrotron_polang_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_nu_synchgrid);                    //AWS20100708
         galaxy.synchrotron_polfra_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_nu_synchgrid);                    //AWS20100708
         galaxy.free_free_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_nu_synchgrid);                    //AWS20100708

      }

      //Perform the actual integration of the skymaps
      if ( 3 == galdef.skymap_format) 
      {
   
         //Create a healpix pixelization to get the angles correct
         Healpix_Base hp(galdef.healpix_order, RING);
   
         INFO("Starting synchrotron skymap integration");
   
#pragma omp parallel for schedule(dynamic) default(shared)
         for ( int ii = 0; ii < hp.Npix(); ++ii ) 
         {
            if ( ii % (hp.Npix()/100) == 0) 
            {
               ostringstream buf;
               buf << "Beginning pixel " << ii << " (" << (100*ii)/hp.Npix() << "%)";// << anisoSkymapOpt.Npix() << " " << outsideGalaxy2D << " " << outsideGalaxy3D;
               INFO(buf.str());
            }
            
            auto co = galaxy.synchrotron_hp_skymap->GetCoordinate( ii );
            double l, b;
            co.getCoordinates(l, b, SM::CoordSys::GAL);
            l *= utl::kConvertRadiansToDegrees;
            b *= utl::kConvertRadiansToDegrees;

            vector< vector< std::valarray<double> > > emiss;

            emiss = emissInt->integrate(l, b, emissFuncs, absorption);

            //Sum over all rings
            for (size_t iR = 0; iR < emiss[0].size(); ++iR) 
            {
               for ( int iE=0; iE < galaxy.n_nu_synchgrid; ++iE ) {
                  galaxy.synchrotron_hp_skymap->GetReference(ii,iE)   += emiss[0][iR][iE];
                  galaxy.synchrotron_Q_hp_skymap->GetReference(ii,iE) += emiss[1][iR][iE];
                  galaxy.synchrotron_U_hp_skymap->GetReference(ii,iE) += emiss[2][iR][iE];
               }
            }
            
            for ( int iE=0; iE < galaxy.n_nu_synchgrid; ++iE ) 
            {
               galaxy.synchrotron_P_hp_skymap->GetReference(ii,iE) = sqrt(galaxy.synchrotron_Q_hp_skymap->GetReference(ii,iE)*galaxy.synchrotron_Q_hp_skymap->GetReference(ii,iE) + galaxy.synchrotron_U_hp_skymap->GetReference(ii,iE)*galaxy.synchrotron_U_hp_skymap->GetReference(ii,iE));
               galaxy.synchrotron_polang_hp_skymap->GetReference(ii,iE) = 0.5*atan2(galaxy.synchrotron_U_hp_skymap->GetReference(ii,iE),galaxy.synchrotron_Q_hp_skymap->GetReference(ii,iE)) * 180./M_PI;
               galaxy.synchrotron_polfra_hp_skymap->GetReference(ii,iE) = galaxy.synchrotron_P_hp_skymap->GetReference(ii,iE)/galaxy.synchrotron_hp_skymap->GetReference(ii,iE);
            }

            if ( galdef.free_free_absorption >= 1 ) 
            {
               for (size_t iR = 0; iR < emiss[0].size(); ++iR) 
               {
                  for ( int iE=0; iE < galaxy.n_nu_synchgrid; ++iE ) {
                     galaxy.free_free_hp_skymap->GetReference(ii,iE) += emiss[3][iR][iE];
                  }
               }
            }

         }
      
      } 
      else 
      { // standard (non-HEALPIX) format IMOS20080114

#pragma omp parallel for schedule(dynamic) default(shared)
         for ( int i_long=0; i_long<galaxy.n_long; i_long++ ) 
         {
            for ( int i_lat =0; i_lat <galaxy.n_lat; i_lat++ ) 
            {
               const double l=galaxy.long_min + i_long*galaxy.d_long;
               const double b=galaxy. lat_min + i_lat *galaxy.d_lat ;

               vector< vector< std::valarray<double> > > emiss;

               emiss = emissInt->integrate(l, b, emissFuncs);

               //Sum over all rings
               for (size_t iR = 0; iR < emiss[0].size(); ++iR) 
               {

                  for ( int iE=0; iE < galaxy.n_nu_synchgrid; ++iE ) 
                  {
                     galaxy.synchrotron_skymap.d2[i_long][i_lat].s[iE] += emiss[0][iR][iE];
                     galaxy.synchrotron_Q_skymap.d2[i_long][i_lat].s[iE] += emiss[1][iR][iE];
                     galaxy.synchrotron_U_skymap.d2[i_long][i_lat].s[iE] += emiss[2][iR][iE];
                  }

               }
               
               for ( int iE=0; iE < galaxy.n_nu_synchgrid; ++iE ) 
               {
                  galaxy.synchrotron_P_skymap.d2[i_long][i_lat].s[iE] = sqrt(galaxy.synchrotron_Q_skymap.d2[i_long][i_lat].s[iE]*galaxy.synchrotron_Q_skymap.d2[i_long][i_lat].s[iE] + galaxy.synchrotron_U_skymap.d2[i_long][i_lat].s[iE]*galaxy.synchrotron_U_skymap.d2[i_long][i_lat].s[iE]);
                  galaxy.synchrotron_polang_skymap.d2[i_long][i_lat].s[iE] = 0.5*atan2(galaxy.synchrotron_U_skymap.d2[i_long][i_lat].s[iE],galaxy.synchrotron_Q_skymap.d2[i_long][i_lat].s[iE]) * 180./M_PI;
                  galaxy.synchrotron_polfra_skymap.d2[i_long][i_lat].s[iE] = galaxy.synchrotron_P_skymap.d2[i_long][i_lat].s[iE]/galaxy.synchrotron_skymap.d2[i_long][i_lat].s[iE];
               }

               if ( galdef.free_free_absorption >= 1 ) 
               {
                  for (size_t iR = 0; iR < emiss[0].size(); ++iR) 
                  {
                     for ( int iE=0; iE < galaxy.n_nu_synchgrid; ++iE ) 
                     {
                        galaxy.free_free_skymap.d2[i_long][i_lat].s[iE] += emiss[3][iR][iE];
                     }
                  }
               }
            }
         }
      }

      TIME_SUBROUTINE(store_synch_skymap);

      if (3 == galdef.skymap_format) 
      {

         galaxy.synchrotron_hp_skymap.reset(nullptr);
         galaxy.synchrotron_Q_hp_skymap.reset(nullptr);
         galaxy.synchrotron_U_hp_skymap.reset(nullptr);
         galaxy.synchrotron_P_hp_skymap.reset(nullptr);
         galaxy.synchrotron_polang_hp_skymap.reset(nullptr);
         galaxy.synchrotron_polfra_hp_skymap.reset(nullptr);
         galaxy.free_free_hp_skymap.reset(nullptr);

      } 
      else 
      {  

         galaxy.synchrotron_skymap.delete_array();
         galaxy.synchrotron_Q_skymap.delete_array(); //AWS20100708
         galaxy.synchrotron_U_skymap.delete_array(); //AWS20100708
         galaxy.synchrotron_P_skymap.delete_array(); //AWS20100708
         galaxy.synchrotron_polang_skymap.delete_array(); //AWS20100708
         galaxy.synchrotron_polfra_skymap.delete_array(); //AWS20100708
         galaxy.free_free_skymap.delete_array(); //AWS20100708

      }  

   }

   emissFuncs.resize(0);
    
   INFO( "Exit" );

}


inline void calculatePixel(const double l, const double b, const double dl, 
      const double cameraTheta, const double cameraR,
      const vector<double> &Rbins,
      const std::unique_ptr< SM::LOSintegrator<double> > &gasInt, 
      const std::unique_ptr< SM::LOSintegrator< valarray<double> > > &emissInt, 
      const vector< std::unique_ptr< SM::LOSfunction<double> > > &gasFuncs, 
      const vector< std::unique_ptr< SM::LOSfunction< valarray<double> > > > &emissFuncs,
      const std::unique_ptr< SM::LOSfunction<valarray<double> > > &absorption,
      vector< vector< double > >&gas, vector< vector< valarray<double> > > &emiss)
{

   if (gasFuncs.size() > 0)
      gas = gasInt->integrate(l, b, gasFuncs);
   else
      gas.resize(0);

   emiss = emissInt->integrate(l, b, emissFuncs, absorption);

   //Check to see if the center of Pixel misses an annulus but the edge doesn't
   //Begin by calculating the angle the los makes to the line to the GC in a plane parallel to the Galactic plane
   double alpha = 360 - l - cameraTheta;
   if (alpha < 0)
      alpha += 360;

   //No need to do this if we only have a single bin in radius
   if ( ( alpha < 90 || alpha > 270 ) && dl > 0 && Rbins.size() > 2 ) {

      int sign = 1;
      if (alpha > 270) 
         sign = -1;

      const double rTang = cameraR*fabs(sin(alpha*utl::kConvertDegreesToRadians));
      const double rTangMin = cameraR*fabs(sin((alpha-sign*dl/2.)*utl::kConvertDegreesToRadians));

      //If the annuli differ we have to take action 
      const size_t iR = findAnnulus(Rbins,rTang);
      const size_t iRMin = findAnnulus(Rbins,rTangMin);

      //Never do this for the outermost bin if camera is outside the Galaxy
      if (iRMin < iR && iR < Rbins.size()-1 && emiss[0][iR].sum() > 0) {
         const double newl = l - sign*dl/2.;

         //The weight of the new line of sight, assuming square pixels
         double tanalpha = asin(Rbins[iR]/cameraR)/utl::kConvertDegreesToRadians;
         if ( sign == -1 )
            tanalpha = 360 - tanalpha;
         double tanl = 360 - tanalpha - cameraTheta;
         if (tanl < 0)
            tanl += 360;
         const double weight = (dl/2. - fabs(l-tanl)) / dl;

         //Calculate the los integral and add to the total with the correct weight
         vector< vector< double > > tmpgas; 
         if (gasFuncs.size() > 0)
            tmpgas = gasInt->integrate(newl, b, gasFuncs);

         for (size_t i = 0; i < gas.size(); ++i) {
            for (size_t j = 0; j < gas[i].size(); ++j) 
               gas[i][j] = (1-weight)*gas[i][j] + weight*tmpgas[i][j];
         }

         vector< vector< valarray<double> > > tmpemiss;
         tmpemiss = emissInt->integrate(newl, b, emissFuncs);

         for (size_t i = 0; i < emiss.size(); ++i) {
            for (size_t j = 0; j < emiss[i].size(); ++j) {
               emiss[i][j] = (1-weight)*emiss[i][j] + weight*tmpemiss[i][j];
            }
         }

      }

   }
}
