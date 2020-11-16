
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_IC_skymap.cc *                            galprop package * 4/20/2006
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// generate inverse Compton skymaps 

// list of selectable debugs:  IMOS20080114
// galdef.verbose==-457 is the return to the old method of integration to compare with older versions

#include"galprop_classes.h"
#include"galprop_internal.h"

#include <ErrorLogger.h>
#include <PhysicalConstants.h>
#include <RadiationField.h>
#include <Skymap.h>

#include <config.h>

#include <constants.h>

#include <algorithm>
#include <cassert>
#include <sstream>
#include <valarray>
#include <vector>
#include <stdexcept>

using namespace std;

static double ClosestDistance2D(const double l, const double b, 
				const double rMax, 
				const double zMin, const double zMax,
				const valarray<double>& cameraLocation) { 

  double s = 0;

  const double cameraX = cameraLocation[0];
  const double cameraY = cameraLocation[1];
  const double cameraZ = cameraLocation[2];
  
  const double lRad = l*utl::kConvertDegreesToRadians, bRad = b*utl::kConvertDegreesToRadians;

  const double sinb = sin(bRad), cosb = cos(bRad), sinl = sin(lRad), cosl = cos(lRad);  

  // Intersection with cylindrical volume 
  // Method: solve quadratic equation for cylinder with ray intersection
  // Let point X = (xo + tpx, y0 + tpy, z0 + tpz) intersect with 
  // cylinder with equation x^2 + y^2 = rMax^2, z = +/- zMax.

  const pointing co(vec3(-cosb*cosl, -cosb*sinl, sinb));

  vec3 dir = co.to_vec3();

  const double pX = dir.x, pY = dir.y, pZ = dir.z;
    
  const double cameraR2 = cameraX*cameraX + cameraY*cameraY;

  const double C = (cameraR2 - rMax*rMax);
  
  const double A = pX*pX + pY*pY;

  const double B = 2.*(pX*cameraX + pY*cameraY);
    
  const double discriminant = B*B - 4.*A*C;

  // No real roots, no intersection with cylinder
 
  if (discriminant < 0.) 
    return 0;
  
  const double t1 = 0.5/A*(-B + sqrt(discriminant));
  
  const double t2 = 0.5/A*(-B - sqrt(discriminant));
  
  if (t1 < 0 && t2 < 0)
    return 0;
  
  const double z1 = cameraZ + t1*pZ, z2 = cameraZ + t2*pZ;
  
  if (cameraZ >= 0 && z1 > zMax && z2 > zMax)
    return 0;
  
  if (cameraZ <= 0 && z1 < zMin && z2 < zMin)
    return 0;
  
  if (cameraZ >= zMin && cameraZ <= zMax) {
    
    s = std::min(t1, t2);
    
  } else if (cameraZ > 0. && pZ < 0.) {
    
    s = (zMax - cameraZ)/pZ;
    
  } else if (cameraZ < 0. && pZ > 0.) {
    
    s = (zMin - cameraZ)/pZ;
    
  } else 
    return 0;
  
  const double x = cameraX - s*cosb*cosl, y = cameraY - s*cosb*sinl, r = sqrt(x*x + y*y);
  
  s = (r > rMax ? std::min(t1, t2) : s);
  
  //cout << l << " " << b << " " << A << " " << B << " " << C << " " << t1 << " " << t2 << " " << z1 << " " << z2 << " " << s << " " << cameraX - s*cosb*cosl << " " << cameraY - s*cosb*sinl << " " << cameraZ + s*sinb << endl;
  
  return s;
    
}

static double ClosestDistance3D(const double l, const double b, 
				const double xMin, const double xMax,
				const double yMin, const double yMax,
				const double zMin, const double zMax,
				const valarray<double>& cameraLocation) { 

  double s = 0;

  const double cameraX = cameraLocation[0];
  const double cameraY = cameraLocation[1];
  const double cameraZ = cameraLocation[2];
  
  const double lRad = l*utl::kConvertDegreesToRadians, bRad = b*utl::kConvertDegreesToRadians;

  const double sinb = sin(bRad), cosb = cos(bRad), sinl = sin(lRad), cosl = cos(lRad);  

  // Intersection with box volume
  // Method: equation for plane with ray intersection for all faces of box
  // Let point X = (xo + tpx, y0 + tpy, z0 + tpz) intersect with 
  // plane with equation Ax + By + Cz + D for all faces of box.
  // Use algorithm from `An efficient and robust ray-box intersection
  // algorithm' -- Williams et al. http://www.cs.utah.edu/~awilliam/box/
     
  const pointing co(vec3(-cosb*cosl, -cosb*sinl, sinb));

  vec3 dir = co.to_vec3();
  
  valarray<double> invDir(0., 3);
  
  invDir[0] = 1./dir.x, invDir[1] = 1./dir.y, invDir[2] = 1./dir.z;
  
  valarray<int> sign(0, 3);
  
  sign[0] = invDir[0] < 0;
  sign[1] = invDir[1] < 0;
  sign[2] = invDir[2] < 0;
  
  valarray< valarray<double> > boxSize;
  
  boxSize.resize(2);
  boxSize[0].resize(3);
  boxSize[1].resize(3);
  
  boxSize[0][0] = xMin, boxSize[0][1] = yMin, boxSize[0][2] = zMin;
  
  boxSize[1][0] = xMax, boxSize[1][1] = yMax, boxSize[1][2] = zMax;
  
  valarray<double> tMin(0., 3), tMax(0., 3);
  
  tMin[0] = (boxSize[sign[0]][0] - cameraLocation[0])*invDir[0];
  tMax[0] = (boxSize[1 - sign[0]][0] - cameraLocation[0])*invDir[0];
  
  tMin[1] = (boxSize[sign[1]][1] - cameraLocation[1])*invDir[1];
  tMax[1] = (boxSize[1 - sign[1]][1] - cameraLocation[1])*invDir[1];
  
  valarray<double> tMinInit(0., 3), tMaxInit(0., 3);
  
  tMinInit = tMin; tMaxInit = tMax;
  
  if ((tMin[0] > tMax[1]) || (tMin[1] > tMax[0]))
    return 0; // No intersection
  
  if (tMin[1] > tMin[0])
    tMin[0] = tMin[1];
  
  if (tMax[1] < tMax[0])
    tMax[0] = tMax[1];
  
  tMin[2] = (boxSize[sign[2]][2] - cameraLocation[2])*invDir[2];
  tMax[2] = (boxSize[1 - sign[2]][2] - cameraLocation[2])*invDir[2];
  
  if ((tMin[0] > tMax[2]) || (tMin[2] > tMax[0]))
    return 0; // No intersection
  
  tMinInit[2] = tMin[2]; tMaxInit[2] = tMax[2];
  
  if (tMin[2] > tMin[0])
    tMin[0] = tMin[2];
  
  if (tMax[2] < tMax[0])
    tMax[0] = tMax[2];
  
  if (tMin[0] < 0 && tMax[0] < 0)
    return 0;
  
  s = tMin[0];

  //cout << l << " " << b << " " << tMin[0] << " " << tMax[0] << " " << cameraX - s*cosb*cosl << " " << cameraY - s*cosb*sinl << " " << cameraZ + s*sinb << endl;

  return s;

}

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
   fgp(gp),
   cameraX(gp.galdef.fCameraLocation[0]),
   cameraY(gp.galdef.fCameraLocation[1]),
   cameraZ(gp.galdef.fCameraLocation[2]),
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
      const unsigned int zBins = fgp.galaxy.n_zgrid;
      
      assert (fgp.galaxy.n_ISRF_components <= 3);

      rGrid.resize(rBins, 0.);

      const double rMax3D = std::max(fgp.galaxy.y_max, fgp.galaxy.x_max);

      rMax = (2 == fgp.gcr[0].n_spatial_dimensions ? fgp.galaxy.r_max : utl::kSqrtTwo*rMax3D); 
      //	   sqrt(galaxy.x_max*galaxy.x_max + galaxy.y_max*galaxy.y_max));

      for (auto iR = 0; iR < rBins; ++iR)
	rGrid[iR] = (2 == fgp.gcr[0].n_spatial_dimensions ? fgp.galaxy.r[iR] : rMax/rBins*iR);
      
      const int targetSkymapOrder = fgp.galdef.ISRF_healpixOrder;
      
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

   for (auto iElectron = 0; iElectron < electronBins; ++iElectron) {

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
         
      for (auto iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {
        
         double anisoSumOpt = 0, isoSumOpt = 0;
         double anisoSumIR = 0, isoSumIR = 0;
         double anisoSumCMB = 0, isoSumCMB = 0; // For checking
      
         for (auto iTarget = 0; iTarget < targetBins; ++iTarget) {
      
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
      
           for (auto iSkymap = 0; iSkymap < cmbSkymap.Npix(); ++iSkymap) {
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


int Galprop::gen_IC_skymap(const valarray<double>& cameraLocation) { //IMOS20080114

  INFO("Entry: gen_IC_skymap");
  
  assert(3 == cameraLocation.size());
  assert(2 == gcr[0].n_spatial_dimensions || 3 == gcr[0].n_spatial_dimensions);

  int stat = 0;

  const unsigned int nComps = galaxy.n_ISRF_components;

  Skymap<double> anisoRatioOpt, anisoRatioIR, anisoRatioCMB;

  const unsigned int nEGammaBins = galaxy.n_E_gammagrid;

  bool anisoValid = false;

  if (galdef.IC_anisotropic) {

    bool haveOMP = false;

#ifdef _OPENMP
    haveOMP = true;
#endif

    if (galdef.ISRF_filetype < 3 || nComps < 3 || !haveOMP) {

      ostringstream buf;
      buf << "Anisotropic IC calculation only available for version 3 ISRF files and later with all 3 ISRF components (optical, infrared, CMB) with openmp enabled. Only isotropic IC will be calculated.";
      INFO(buf.str());

    } else {

      anisoValid = true;

      AnisoICFunction *anisoFunc;
      try{
         anisoFunc = new AnisoICFunction(*this);
      } catch (std::invalid_argument) {
         return 1;
      }

      const unsigned int electronBins = gcr[0].n_pgrid;
      const unsigned int targetBins = galaxy.ISRF[0].n_pgrid;
      const int cosThetaBins = anisoFunc->get_cosThetaBins();

      const unsigned int rBins = anisoFunc->get_rBins();
      const double rMax = anisoFunc->get_rMax();

      const std::valarray<double> &isoCrossSection = anisoFunc->get_isoCrossSection();
      const std::valarray<double> &anisoCrossSection = anisoFunc->get_anisoCrossSection();

      const std::valarray< std::valarray<double> >& optISRF = anisoFunc->get_optISRF();
      const std::valarray< std::valarray<double> >& irISRF = anisoFunc->get_irISRF();

      const std::vector< Skymap<double> >& optAngDist = anisoFunc->get_optAngDist();
      const std::vector< Skymap<double> >& irAngDist = anisoFunc->get_irAngDist();

      const std::valarray<double> &targetE = anisoFunc->get_targetE();
      const std::valarray<double> &gammaE = anisoFunc->get_gammaE();
      const std::valarray<double> &electronE = anisoFunc->get_electronE();
      const std::valarray<double> &rGrid = anisoFunc->get_rGrid();
      const std::valarray<double> &cmbNumberDensity = anisoFunc->get_cmbNumberDensity();

      const Skymap<double> & cmbSkymap = anisoFunc->get_cmbSkymap();
      const Distribution & electrons = anisoFunc->get_electrons();

      //exit(0);

      // Easier to modify location of camera in future ...

      const double cameraX = cameraLocation[0];
      const double cameraY = cameraLocation[1];
      const double cameraZ = cameraLocation[2];

      const double cameraR = sqrt(cameraX*cameraX + cameraY*cameraY);
      const double cameraTheta = atan2(cameraY, cameraX);

      const double ds = 0.25; // kpc -- hard wired so get predictable integration time. Errors due to the coarseness of this step size are effectively cancelled since we calculate the aniso/iso ratio

      const bool outsideGalaxy2D = 
	(cameraR > galaxy.r_max || 
	 cameraZ < galaxy.z_min || cameraZ > galaxy.z_max);
      
      const bool outsideGalaxy3D = 
	(cameraX < galaxy.x_min || cameraX > galaxy.x_max || 
	 cameraY < galaxy.y_min || cameraY > galaxy.y_max || 
	 cameraZ < galaxy.z_min || cameraZ > galaxy.z_max);
    
      // We have to set this depending on the angular size of the object
      // that we are computing the interpolation skymaps for if we are 
      // viewing the object externally (naturally if internal to an object 
      // we can choose this fairly coarsely since we are interpolating 
      // over 4 pi sr). We normally adopt order 2, but if the size of the 
      // object is too small this is obviously too coarse. We use the 
      // heuristic that the interpolation skymap order is increased by 1 for 
      // every factor of 1/8 that the object is of 4 pi sr (so, if the 
      // object is pi/2 sr, increase to order 3, pi/16 sr increase to 
      // order 4, etc.). The calculation of the `solid angle' for the
      // externally viewed object is only approximate, since this is all 
      // we require ...

      const double distance2 = cameraX*cameraX + cameraY*cameraY + cameraZ*cameraZ;

      //const double dr2 = (outsideGalaxy2D || outsideGalaxy3D ? rMax*rMax : 0);

      //const double solidAngleRatio = utl::kFourPi/utl::kTwoPi/(dr2/distance2);

      //cout << "SolidAngleRatio: " << dr2 << " " << distance2 << " " << dr2/distance2 << " " << solidAngleRatio << " " << log10(solidAngleRatio)/log10(4.) << endl;

      // This is set to 2 for internal to a galaxy, and scaled as above using
      // the solid angle ratio. The maximum it can be is healpix 
      // order 8 (for now) 

      //exit(0);

      const unsigned int internalHealpixOrder = galdef.anisoHealpixOrder;//(outsideGalaxy2D || outsideGalaxy3D ? std::min(8, 2 + int(log10(solidAngleRatio)/log10(2.))) : 2);

      Skymap<double> anisoSkymapOpt(internalHealpixOrder, gammaE);
      Skymap<double> anisoSkymapIR(internalHealpixOrder, gammaE);
      Skymap<double> anisoSkymapCMB(internalHealpixOrder, gammaE);
      Skymap<double> isoSkymapOpt(internalHealpixOrder, gammaE);
      Skymap<double> isoSkymapIR(internalHealpixOrder, gammaE);
      Skymap<double> isoSkymapCMB(internalHealpixOrder, gammaE);
      
      anisoSkymapOpt = 0;
      isoSkymapOpt = 0;
      anisoSkymapIR = 0;
      isoSkymapIR = 0;
      anisoSkymapCMB = 0;
      isoSkymapCMB = 0;
      
      // Precompute the inner product for target photon and emitted 
      // photon directions
 
      valarray<vec3> dirGamma(vec3(0, 0, 0), anisoSkymapOpt.Npix());
      
      for (unsigned int i = 0; i < dirGamma.size(); ++i)
	dirGamma[i] = anisoSkymapOpt.pix2coord(i).healpixAng().to_vec3();
       
      valarray<vec3> dirTarget(vec3(0, 0, 0), cmbSkymap.Npix());
      
      for (unsigned int i = 0; i < dirTarget.size(); ++i)
	dirTarget[i] = cmbSkymap.pix2coord(i).healpixAng().to_vec3();

      valarray<double> cosZeta(0., dirTarget.size()*dirGamma.size());
      
      for (unsigned int i = 0; i < dirTarget.size(); ++i)
	for (unsigned int j = 0; j < dirGamma.size(); ++j) {
	  
	  const unsigned int index = i*dirGamma.size() + j;
	  
	  cosZeta[index] = dotprod(-dirTarget[i], dirGamma[j]);
	  
	}
      
      INFO("Calculating interpolation skymaps");

      // Find the total pixels to be integrated over

      int totalPixels = anisoSkymapOpt.Npix();

      /*if (outsideGalaxy2D || outsideGalaxy3D) {

	totalPixels = 0;
	
	for (int iPix = 0; iPix < anisoSkymapOpt.Npix(); ++iPix) {
	
	  SM::Coordinate coord(anisoSkymapOpt.pix2ang(iPix));
	  
	  const double l = coord.l(), b = coord.b();
	
	  // Test if the camera is outside the galaxy/CR region. If so, then we 
	  // set the starting s value so that it is at the boundary of 
	  // the galaxy/CR region, then the integration proceeds as before ...
	  
	  double s =
	    (2 == gcr[0].n_spatial_dimensions && outsideGalaxy2D ? 
	     ClosestDistance2D(l, b, galaxy.r_max, galaxy.z_min, galaxy.z_max, cameraLocation) : 
	     (3 == gcr[0].n_spatial_dimensions && outsideGalaxy3D ?
	      ClosestDistance3D(l, b, galaxy.x_min, galaxy.x_max, galaxy.y_min, galaxy.y_max, galaxy.z_min, galaxy.z_max, cameraLocation) : 0));

	  if (s > 0)
	    ++totalPixels;

	}

      }
      */
#pragma omp parallel for schedule(dynamic) default(shared)
      for (int iPix = 0; iPix < anisoSkymapOpt.Npix(); ++iPix) {

#pragma omp critical 
	{
	  ostringstream buf;
	  buf << "Beginning pixel " << iPix << " (total = " << totalPixels << ")";// << anisoSkymapOpt.Npix() << " " << outsideGalaxy2D << " " << outsideGalaxy3D;
	  INFO(buf.str());
	}
	
	SM::Coordinate coord(anisoSkymapOpt.pix2ang(iPix));
	
	const double l = coord.l(), b = coord.b();
	const double lRad = l*utl::kConvertDegreesToRadians, bRad = b*utl::kConvertDegreesToRadians;

	const double sinb = sin(bRad), cosb = cos(bRad), sinl = sin(lRad), cosl = cos(lRad);
	
	// Test if the camera is outside the galaxy/CR region. If so, then we 
	// set the starting s value so that it is at the boundary of 
	// the galaxy/CR region, then the integration proceeds as before ...
	
	double s = //0;
	(2 == gcr[0].n_spatial_dimensions && outsideGalaxy2D ? 
	   ClosestDistance2D(l, b, galaxy.r_max, galaxy.z_min, galaxy.z_max, cameraLocation) : 
	   (3 == gcr[0].n_spatial_dimensions && outsideGalaxy3D ?
	    ClosestDistance3D(l, b, galaxy.x_min, galaxy.x_max, galaxy.y_min, galaxy.y_max, galaxy.z_min, galaxy.z_max, cameraLocation) : 0));
	
	bool complete = false;

	while (1) {

	  s += ds;
	  
	  const double dx = s*cosb*cosl;
	  const double dy = s*cosb*sinl;
	  const double dz = s*sinb;
	  
	  const double x = cameraX - dx;
	  const double y = cameraY - dy;
	  const double z = cameraZ + dz;
	  
	  const double r = sqrt(x*x + y*y);

	  // These are the loop break criteria. Much simpler encoding of 
	  // the break condition than in the earlier pixel routines.

	  if (2 == gcr[0].n_spatial_dimensions) {

	    complete = (r > galaxy.r_max || 
			z < galaxy.z_min || z > galaxy.z_max);

	  }

	  if (3 == gcr[0].n_spatial_dimensions) {

	    complete = (x < galaxy.x_min || x > galaxy.x_max || 
			y < galaxy.y_min || y > galaxy.y_max || 
			z < galaxy.z_min || z > galaxy.z_max);

	  }

	  //cout << iPix << " " << x << " " << y << " " << z << " " << r << " " << complete << endl;

	  if (complete)
	    break;

	  // Determine angle of gamma rays toward observer from the 
	  // emission region. This is the direction we use for the inner product
	  // in the cross section calculation. Having found the direction we 
	  // just index into the precomputed cosZeta array. Recall also that
	  // we are doing it this way because the coordinate system for the 
	  // ISRF is centred toward the GC due to the azimuthal symmetry.

	  const double cosXi = (z - cameraZ)/s; // = sinb
	  const double cos2Xi = cosXi*cosXi;
	  const double sin2Xi = (cos2Xi > 1. ? 0. : 1. - cos2Xi);
	  const double sinXi = sqrt(sin2Xi); // = fabs(cosb)
	  
	  const double theta = atan2(y, x);
	  
	  const double thetaDelta = theta - cameraTheta;

	  const double rho = sqrt(dx*dx + dy*dy);

	  const double sinChi = -cameraR/rho*sin(theta);
	  const double cosChi = -(r - cameraR*cos(theta))/rho;

	  const vec3 dir(cosChi*sinXi, sinChi*sinXi, cosXi);

	  pointing point(dir);
	  
	  // This is the index into the inner product array calculated earlier

	  const int gammaPixIndex = anisoSkymapOpt.ang2pix(point);

	  valarray<double> anisoEmissivityOpt(0., nEGammaBins), isoEmissivityOpt(0., nEGammaBins);
	  valarray<double> anisoEmissivityIR(0., nEGammaBins), isoEmissivityIR(0., nEGammaBins);
	  valarray<double> anisoEmissivityCMB(0., nEGammaBins), isoEmissivityCMB(0., nEGammaBins);
	  
	  // Essentially the same algorithm for 2D and 3D: if we're at the 
	  // spatial bin boundary, use the electron spectra, target intensities
	  // and energy densities, there. Otherwise, interpolate (bilinear for
	  // 2D, trilinear for 3D) amongst the 4(8) surrounding grid points
	  // to get the electron spectra, etc. Once that is done, it's the 
	  // usual calculation for the isotropic cross section. For the 
	  // anisotropic cross section, there is an inner loop over the 
	  // number of pixels in the target skymap that sums the anisotropic
	  // cross section weighted by the target photon intensity at each
	  // pixel. The anisotropic cross section is evaluated for the 
	  // cosTheta value corresponding to the emission angle toward 
	  // the camera.

	  if (2 == gcr[0].n_spatial_dimensions) {

	    const double dR = (r - galaxy.r_min)/(galaxy.r_max - galaxy.r_min);
	    const double dZ = (z - galaxy.z_min)/(galaxy.z_max - galaxy.z_min);

	    const unsigned int iR = (dR < 1. ? dR*(galaxy.n_rgrid-1) : galaxy.n_rgrid-1);
	    const unsigned int iZ = (dZ < 1. ? dZ*(galaxy.n_zgrid-1) : galaxy.n_zgrid-1);
      
	    if (dR >=1 || dZ >= 1) {

	      const int index = iR*galaxy.n_zgrid + iZ;

	      const Skymap<double>& opt = optAngDist[index];
	      const Skymap<double>& ir = irAngDist[index];

	      for (unsigned int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

		for (unsigned int iElectron = 0; iElectron < electronBins; ++iElectron) {

		  double anisoSumOpt = 0, isoSumOpt = 0;
		  double anisoSumIR = 0, isoSumIR = 0;
		  double anisoSumCMB = 0, isoSumCMB = 0; // For checking

		  for (unsigned int iTarget = 0; iTarget < targetBins; ++iTarget) {

		    const unsigned int isoIndex = iEGamma*electronBins*targetBins + iElectron*targetBins + iTarget;

		    const double xsIso = isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0]);

		    // This is the isotropic cross section summation
		    
		    isoSumOpt += xsIso*optISRF[index][iTarget];

		    isoSumIR += xsIso*irISRF[index][iTarget];

		    isoSumCMB += xsIso*cmbNumberDensity[iTarget];

		    // Inner loop over the target skymaps
		
		    double sumSkymapOpt = 0, sumSkymapIR = 0, sumSkymapCMB = 0;
    
		    for (unsigned int iSkymap = 0; iSkymap < opt.Npix(); ++iSkymap) {

		      const unsigned int cosZetaIndex = iSkymap*dirGamma.size() + gammaPixIndex;

		      const unsigned int cosThetaIndex = 
			((1. + cosZeta[cosZetaIndex]) < 2. ? int(0.5*(1. + cosZeta[cosZetaIndex])*cosThetaBins) : cosThetaBins-1); 

		      const unsigned int anisoIndex = iEGamma*electronBins*targetBins*cosThetaBins + iElectron*targetBins*cosThetaBins + iTarget*cosThetaBins + cosThetaIndex;
		      
		      sumSkymapOpt += anisoCrossSection[anisoIndex]*
			opt[iSkymap][iTarget];
		      
		      sumSkymapIR += anisoCrossSection[anisoIndex]*
			ir[iSkymap][iTarget];
		      
		      sumSkymapCMB += anisoCrossSection[anisoIndex]*
			cmbSkymap[iSkymap][iTarget];

		    }
		     
		    anisoSumOpt += sumSkymapOpt*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumIR += sumSkymapIR*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumCMB += sumSkymapCMB*targetE[iTarget]*log(targetE[1]/targetE[0]);

		  }

		  const double elecSum = electrons.d2[iR][iZ].s[iElectron]*electronE[iElectron]*log(electronE[1]/electronE[0]);
		  
		  isoEmissivityOpt[iEGamma] += isoSumOpt*elecSum;
		  isoEmissivityIR[iEGamma] += isoSumIR*elecSum;
		  isoEmissivityCMB[iEGamma] += isoSumCMB*elecSum;

		  anisoEmissivityOpt[iEGamma] += anisoSumOpt*elecSum;
		  anisoEmissivityIR[iEGamma] += anisoSumIR*elecSum;
		  anisoEmissivityCMB[iEGamma] += anisoSumCMB*elecSum;

		}

	      }

	    } else {

	      // Bilinear interpolation

	      const double rCoeff = (r - galaxy.r[iR])/(galaxy.r[iR+1] - galaxy.r[iR]);
	      const double zCoeff = (z - galaxy.z[iZ])/(galaxy.z[iZ+1] - galaxy.z[iZ]);

	      const unsigned int index1 = iR*galaxy.n_zgrid + iZ;
	      const unsigned int index2 = (iR+1)*galaxy.n_zgrid + iZ;
	      const unsigned int index3 = iR*galaxy.n_zgrid + (iZ+1);
	      const unsigned int index4 = (iR+1)*galaxy.n_zgrid + (iZ+1);

	      const Skymap<double>& opt1 = optAngDist[index1];
	      const Skymap<double>& ir1 = irAngDist[index1];

	      const Skymap<double>& opt2 = optAngDist[index2];
	      const Skymap<double>& ir2 = irAngDist[index2];

	      const Skymap<double>& opt3 = optAngDist[index3];
	      const Skymap<double>& ir3 = irAngDist[index3];

	      const Skymap<double>& opt4 = optAngDist[index4];
	      const Skymap<double>& ir4 = irAngDist[index4];

	      for (unsigned int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

		for (unsigned int iElectron = 0; iElectron < electronBins; ++iElectron) {
		  
		  double anisoSumOpt = 0, isoSumOpt = 0;
		  double anisoSumIR = 0, isoSumIR = 0;
		  double anisoSumCMB = 0, isoSumCMB = 0; // For checking

		  for (unsigned int iTarget = 0; iTarget < targetBins; ++iTarget) {

		    const unsigned int isoIndex = iEGamma*electronBins*targetBins + iElectron*targetBins + iTarget;

		    // This is the isotropic cross section summation
		    
		    const double xsIso = isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0]);

		    isoSumOpt += xsIso*
		      (optISRF[index1][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
		       optISRF[index2][iTarget]*rCoeff*(1. - zCoeff) +
		       optISRF[index3][iTarget]*(1. - rCoeff)*zCoeff +
		       optISRF[index4][iTarget]*rCoeff*zCoeff);
		   
		    isoSumIR += xsIso*
		      (irISRF[index1][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
		       irISRF[index2][iTarget]*rCoeff*(1. - zCoeff) +
		       irISRF[index3][iTarget]*(1. - rCoeff)*zCoeff +
		       irISRF[index4][iTarget]*rCoeff*zCoeff);
						
		    isoSumCMB += xsIso*
		      cmbNumberDensity[iTarget];
		
		    // Inner loop over the target skymaps

		    double sumSkymapOpt = 0, sumSkymapIR = 0, sumSkymapCMB = 0;
    
		    for (unsigned int iSkymap = 0; iSkymap < opt1.Npix(); ++iSkymap) {

		      const unsigned int cosZetaIndex = iSkymap*dirGamma.size() + gammaPixIndex;

		      const unsigned int cosThetaIndex = 
			((1. + cosZeta[cosZetaIndex]) < 2. ? int(0.5*(1. + cosZeta[cosZetaIndex])*cosThetaBins) : cosThetaBins-1); 

		      const unsigned int anisoIndex = iEGamma*electronBins*targetBins*cosThetaBins + iElectron*targetBins*cosThetaBins + iTarget*cosThetaBins + cosThetaIndex;
		    
		      sumSkymapOpt += anisoCrossSection[anisoIndex]*
			(opt1[iSkymap][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
			 opt2[iSkymap][iTarget]*rCoeff*(1. - zCoeff) +
			 opt3[iSkymap][iTarget]*(1. - rCoeff)*zCoeff +
			 opt4[iSkymap][iTarget]*rCoeff*zCoeff);

		      sumSkymapIR += anisoCrossSection[anisoIndex]*
			(ir1[iSkymap][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
			 ir2[iSkymap][iTarget]*rCoeff*(1. - zCoeff) +
			 ir3[iSkymap][iTarget]*(1. - rCoeff)*zCoeff +
			 ir4[iSkymap][iTarget]*rCoeff*zCoeff);
		      
		      sumSkymapCMB += anisoCrossSection[anisoIndex]*
			cmbSkymap[iSkymap][iTarget];

		    }
		     
		    anisoSumOpt += sumSkymapOpt*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumIR += sumSkymapIR*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumCMB += sumSkymapCMB*targetE[iTarget]*log(targetE[1]/targetE[0]);

		  }

		  // All 4 grid points for the bilinear interpolation
		  
		  const double elecSum =
		    (electrons.d2[iR][iZ].s[iElectron]*(1. - rCoeff)*(1. - zCoeff) +
		     electrons.d2[iR+1][iZ].s[iElectron]*rCoeff*(1. - zCoeff) +
		     electrons.d2[iR][iZ+1].s[iElectron]*(1. - rCoeff)*zCoeff +
		     electrons.d2[iR+1][iZ+1].s[iElectron]*rCoeff*zCoeff)*
		    electronE[iElectron]*log(electronE[1]/electronE[0]);
		  
		  isoEmissivityOpt[iEGamma] += isoSumOpt*elecSum;
		  isoEmissivityIR[iEGamma] += isoSumIR*elecSum;
		  isoEmissivityCMB[iEGamma] += isoSumCMB*elecSum;

		  anisoEmissivityOpt[iEGamma] += anisoSumOpt*elecSum;
		  anisoEmissivityIR[iEGamma] += anisoSumIR*elecSum;
		  anisoEmissivityCMB[iEGamma] += anisoSumCMB*elecSum;

		  //cout << iEGamma << " " << iElectron << " " << isoSumCMB << endl;

		}

	      }

	    }

	  }

	  if (3 == gcr[0].n_spatial_dimensions) {

	    const double dX = (x - galaxy.x_min)/(galaxy.x_max - galaxy.x_min);
	    const double dY = (y - galaxy.y_min)/(galaxy.y_max - galaxy.y_min);
	    const double dZ = (z - galaxy.z_min)/(galaxy.z_max - galaxy.z_min);

	    const unsigned int iX = (dX < 1. ? dX*(galaxy.n_xgrid-1) : galaxy.n_xgrid-1);
	    const unsigned int iY = (dY < 1. ? dY*(galaxy.n_ygrid-1) : galaxy.n_ygrid-1);
	    const unsigned int iZ = (dZ < 1. ? dZ*(galaxy.n_zgrid-1) : galaxy.n_zgrid-1);

	    // Radiation field is sampled only on a R,z grid due to 
	    // computational limitations (for now).

	    const double dR = r/rMax;
	    
	    const unsigned int iR = (dR < 1. ? dR*(rBins-1) : rBins-1);	    
   
	    if (dX >= 1 || dY >= 1 || dZ >= 1 || iR >= 1) {

	      const int index = iR*galaxy.n_zgrid + iZ;

	      const Skymap<double>& opt = optAngDist[index];
	      const Skymap<double>& ir = irAngDist[index];

	      for (unsigned int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

		for (unsigned int iElectron = 0; iElectron < electronBins; ++iElectron) {

		  double anisoSumOpt = 0, isoSumOpt = 0;
		  double anisoSumIR = 0, isoSumIR = 0;
		  double anisoSumCMB = 0, isoSumCMB = 0; // For checking

		  for (unsigned int iTarget = 0; iTarget < targetBins; ++iTarget) {

		    const unsigned int isoIndex = iEGamma*electronBins*targetBins + iElectron*targetBins + iTarget;

		    // This is the isotropic cross section summation
		    
		    const double xsIso = isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0]);

		    isoSumOpt += xsIso*optISRF[index][iTarget];

		    isoSumIR += xsIso*irISRF[index][iTarget];

		    isoSumCMB += xsIso*cmbNumberDensity[iTarget];

		    // Inner loop over target skymaps
		
		    double sumSkymapOpt = 0, sumSkymapIR = 0, sumSkymapCMB = 0;
    
		    for (unsigned int iSkymap = 0; iSkymap < opt.Npix(); ++iSkymap) {

		      const unsigned int cosZetaIndex = iSkymap*dirGamma.size() + gammaPixIndex;

		      const unsigned int cosThetaIndex = 
			((1. + cosZeta[cosZetaIndex]) < 2. ? int(0.5*(1. + cosZeta[cosZetaIndex])*cosThetaBins) : cosThetaBins-1); 

		      const unsigned int anisoIndex = iEGamma*electronBins*targetBins*cosThetaBins + iElectron*targetBins*cosThetaBins + iTarget*cosThetaBins + cosThetaIndex;
		      
		      sumSkymapOpt += anisoCrossSection[anisoIndex]*
			opt[iSkymap][iTarget];
		      
		      sumSkymapIR += anisoCrossSection[anisoIndex]*
			ir[iSkymap][iTarget];
		      
		      sumSkymapCMB += anisoCrossSection[anisoIndex]*
			cmbSkymap[iSkymap][iTarget];

		    }
		     
		    anisoSumOpt += sumSkymapOpt*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumIR += sumSkymapIR*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumCMB += sumSkymapCMB*targetE[iTarget]*log(targetE[1]/targetE[0]);

		  }

		  const double elecSum = electrons.d3[iX][iY][iZ].s[iElectron]*electronE[iElectron]*log(electronE[1]/electronE[0]);
		  
		  isoEmissivityOpt[iEGamma] += isoSumOpt*elecSum;
		  isoEmissivityIR[iEGamma] += isoSumIR*elecSum;
		  isoEmissivityCMB[iEGamma] += isoSumCMB*elecSum;

		  anisoEmissivityOpt[iEGamma] += anisoSumOpt*elecSum;
		  anisoEmissivityIR[iEGamma] += anisoSumIR*elecSum;
		  anisoEmissivityCMB[iEGamma] += anisoSumCMB*elecSum;

		}

	      }

	    } else {

	      // Trilinear interpolation

	      const double xCoeff = (x - galaxy.x[iX])/(galaxy.x[iX+1] - galaxy.x[iX]);

	      const double yCoeff = (y - galaxy.y[iY])/(galaxy.y[iY+1] - galaxy.y[iY]);

	      const double zCoeff = (z - galaxy.z[iZ])/(galaxy.z[iZ+1] - galaxy.z[iZ]);

	      const double rCoeff = (r - rGrid[iR])/(rGrid[iR+1] - rGrid[iR]);
 
	      const unsigned int index1 = iR*galaxy.n_zgrid + iZ;
	      const unsigned int index2 = (iR+1)*galaxy.n_zgrid + iZ;
	      const unsigned int index3 = iR*galaxy.n_zgrid + (iZ+1);
	      const unsigned int index4 = (iR+1)*galaxy.n_zgrid + (iZ+1);

	      const Skymap<double>& opt1 = optAngDist[index1];
	      const Skymap<double>& ir1 = irAngDist[index1];

	      const Skymap<double>& opt2 = optAngDist[index2];
	      const Skymap<double>& ir2 = irAngDist[index2];

	      const Skymap<double>& opt3 = optAngDist[index3];
	      const Skymap<double>& ir3 = irAngDist[index3];

	      const Skymap<double>& opt4 = optAngDist[index4];
	      const Skymap<double>& ir4 = irAngDist[index4];

	      for (unsigned int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

		for (unsigned int iElectron = 0; iElectron < electronBins; ++iElectron) {
		  
		  double anisoSumOpt = 0, isoSumOpt = 0;
		  double anisoSumIR = 0, isoSumIR = 0;
		  double anisoSumCMB = 0, isoSumCMB = 0; // For checking

		  for (unsigned int iTarget = 0; iTarget < targetBins; ++iTarget) {

		    const unsigned int isoIndex = iEGamma*electronBins*targetBins + iElectron*targetBins + iTarget;

		    // This is the isotropic cross section summation

		    const double xsIso = isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0]);

		    isoSumOpt += xsIso*
		      (optISRF[index1][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
		       optISRF[index2][iTarget]*rCoeff*(1. - zCoeff) +
		       optISRF[index3][iTarget]*(1. - rCoeff)*zCoeff +
		       optISRF[index4][iTarget]*rCoeff*zCoeff);
		   
		    isoSumIR += xsIso*
		      (irISRF[index1][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
		       irISRF[index2][iTarget]*rCoeff*(1. - zCoeff) +
		       irISRF[index3][iTarget]*(1. - rCoeff)*zCoeff +
		       irISRF[index4][iTarget]*rCoeff*zCoeff);
						
		    isoSumCMB += xsIso*cmbNumberDensity[iTarget];

		    // Inner loop over target skymaps
		
		    double sumSkymapOpt = 0, sumSkymapIR = 0, sumSkymapCMB = 0;
    
		    for (unsigned int iSkymap = 0; iSkymap < opt1.Npix(); ++iSkymap) {

		      const unsigned int cosZetaIndex = iSkymap*dirGamma.size() + gammaPixIndex;

		      const unsigned int cosThetaIndex = 
			((1. + cosZeta[cosZetaIndex]) < 2. ? int(0.5*(1. + cosZeta[cosZetaIndex])*cosThetaBins) : cosThetaBins-1); 

		      const unsigned int anisoIndex = iEGamma*electronBins*targetBins*cosThetaBins + iElectron*targetBins*cosThetaBins + iTarget*cosThetaBins + cosThetaIndex;
		    
		      sumSkymapOpt += anisoCrossSection[anisoIndex]*
			(opt1[iSkymap][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
			 opt2[iSkymap][iTarget]*rCoeff*(1. - zCoeff) +
			 opt3[iSkymap][iTarget]*(1. - rCoeff)*zCoeff +
			 opt4[iSkymap][iTarget]*rCoeff*zCoeff);

		      sumSkymapIR += anisoCrossSection[anisoIndex]*
			(ir1[iSkymap][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
			 ir2[iSkymap][iTarget]*rCoeff*(1. - zCoeff) +
			 ir3[iSkymap][iTarget]*(1. - rCoeff)*zCoeff +
			 ir4[iSkymap][iTarget]*rCoeff*zCoeff);
		      
		      sumSkymapCMB += anisoCrossSection[anisoIndex]*
			cmbSkymap[iSkymap][iTarget];

		    }
		     
		    anisoSumOpt += sumSkymapOpt*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumIR += sumSkymapIR*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumCMB += sumSkymapCMB*targetE[iTarget]*log(targetE[1]/targetE[0]);

		  }
		  
		  // All 8 grid points for the trilinear interpolation

		  const double elecSum =
		    (electrons.d3[iX][iY][iZ].s[iElectron]*(1. - xCoeff)*(1. - yCoeff)*(1. - zCoeff) +
		     electrons.d3[iX+1][iY][iZ].s[iElectron]*xCoeff*(1. - yCoeff)*(1. - zCoeff) +
		     electrons.d3[iX][iY+1][iZ].s[iElectron]*(1. - xCoeff)*yCoeff*(1. - zCoeff) +
		     electrons.d3[iX+1][iY+1][iZ].s[iElectron]*xCoeff*yCoeff*(1. - zCoeff) +
		     electrons.d3[iX][iY][iZ+1].s[iElectron]*(1. - xCoeff)*(1. - yCoeff)*zCoeff +
		     electrons.d3[iX+1][iY][iZ+1].s[iElectron]*xCoeff*(1. - yCoeff)*zCoeff +
		     electrons.d3[iX][iY+1][iZ+1].s[iElectron]*(1. - xCoeff)*yCoeff*zCoeff +
		     electrons.d3[iX+1][iY+1][iZ+1].s[iElectron]*xCoeff*yCoeff*zCoeff)*
		    electronE[iElectron]*log(electronE[1]/electronE[0]);
		  
		  isoEmissivityOpt[iEGamma] += isoSumOpt*elecSum;
		  isoEmissivityIR[iEGamma] += isoSumIR*elecSum;
		  isoEmissivityCMB[iEGamma] += isoSumCMB*elecSum;

		  anisoEmissivityOpt[iEGamma] += anisoSumOpt*elecSum;
		  anisoEmissivityIR[iEGamma] += anisoSumIR*elecSum;
		  anisoEmissivityCMB[iEGamma] += anisoSumCMB*elecSum;

		  //cout << iEGamma << " " << iElectron << " " << isoSumCMB << endl;

		}

	      }

	    }

	  }

	  for (int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

	    anisoSkymapOpt[iPix][iEGamma] += anisoEmissivityOpt[iEGamma]*ds*kpc2cm;
	    anisoSkymapIR[iPix][iEGamma] += anisoEmissivityIR[iEGamma]*ds*kpc2cm;
	    anisoSkymapCMB[iPix][iEGamma] += anisoEmissivityCMB[iEGamma]*ds*kpc2cm;
	    
	    isoSkymapOpt[iPix][iEGamma] += isoEmissivityOpt[iEGamma]*ds*kpc2cm;
	    isoSkymapIR[iPix][iEGamma] += isoEmissivityIR[iEGamma]*ds*kpc2cm;
	    isoSkymapCMB[iPix][iEGamma] += isoEmissivityCMB[iEGamma]*ds*kpc2cm;
	    
	  }

	}

	//#pragma omp critical 
	//{
	//ostringstream buf;
	//buf << "Finished pixel " << iPix << " (total = " << totalPixels << ") " << anisoSkymapOpt.Npix() << " " << outsideGalaxy2D << " " << outsideGalaxy3D;
	//INFO(buf.str());
	//}

	/*for (unsigned int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

	  SM::Coordinate coord(anisoSkymapOpt.pix2ang(iPix));
	
	  const double l = coord.l(), b = coord.b();
	  
	  cout << "EGamma: " 
	       << iPix << " " 
	       << l << " " << b << " " 
	       << iEGamma << " " 
	       << isoSkymapCMB[iPix][iEGamma] << " "
	    //<< anisoSkymapCMB[iPix][iEGamma] << " "
	       << (isoSkymapCMB[iPix][iEGamma] > 0. ? anisoSkymapCMB[iPix][iEGamma]/isoSkymapCMB[iPix][iEGamma] : 0.) << " "
	       << isoSkymapIR[iPix][iEGamma] << " "
	    // << anisoSkymapIR[iPix][iEGamma] << " "
	       << (isoSkymapIR[iPix][iEGamma] > 0. ? anisoSkymapIR[iPix][iEGamma]/isoSkymapIR[iPix][iEGamma] : 0.) << " "
	       << isoSkymapOpt[iPix][iEGamma] << " "
	    // << anisoSkymapOpt[iPix][iEGamma] << " "
	       << (isoSkymapOpt[iPix][iEGamma] > 0. ? anisoSkymapOpt[iPix][iEGamma]/isoSkymapOpt[iPix][iEGamma] : 0.) << endl;

	       }*/

      }
 
      Skymap<double> optRatio(internalHealpixOrder, gammaE);
      Skymap<double> irRatio(internalHealpixOrder, gammaE);
      Skymap<double> cmbRatio(internalHealpixOrder, gammaE);

      for (int iPix = 0; iPix < cmbRatio.Npix(); ++iPix) {
	
	for (int iSpec = 0; iSpec < gammaE.size(); ++iSpec) {

	  optRatio[iPix][iSpec] = 
	    (isoSkymapOpt[iPix][iSpec] > 0 && anisoSkymapOpt[iPix][iSpec] > 0 ? anisoSkymapOpt[iPix][iSpec]/isoSkymapOpt[iPix][iSpec] : 0);

	  irRatio[iPix][iSpec] = 
	    (isoSkymapIR[iPix][iSpec] > 0 && anisoSkymapIR[iPix][iSpec] > 0 ? anisoSkymapIR[iPix][iSpec]/isoSkymapIR[iPix][iSpec] : 0);

	  cmbRatio[iPix][iSpec] = 
	    (isoSkymapCMB[iPix][iSpec] > 0 && anisoSkymapCMB[iPix][iSpec] > 0 ? anisoSkymapCMB[iPix][iSpec]/isoSkymapCMB[iPix][iSpec] : 0);
	  
	}

      }
 
      // The skymaps calculated above are most likely done for a healpix 
      // pixelisation that is considerably coarser than the pixelisation 
      // used for the skymaps to be output to disk. Here we interpolate 
      // the skymaps calculated above to the higher order pixelisation 
      // and obtain a ratio map as a function of energy. This can be applied 
      // to the isotropic IC skymap calculated below to obtain the full 
      // anisotropic skymap for the higher pixelisations (order 6, 7, ..) 
      // typically used in calculations.

      const unsigned int healpixOrder = galdef.healpix_order;
 
      anisoRatioOpt = optRatio.interpolate(healpixOrder);
      anisoRatioIR = irRatio.interpolate(healpixOrder);
      anisoRatioCMB = cmbRatio.interpolate(healpixOrder);

      //const string fN = configure.fOutputDirectory + configure.fOutputPrefix;
    
      //isoSkymapCMB.write(fN + "iso_cmb_hp.fits.gz");
      //anisoSkymapCMB.write(fN + "aniso_cmb_hp.fits.gz");
      //cmbRatio.write(fN + "ratio_cmb_hp.fits.gz");

      //isoSkymapIR.write(fN + "iso_ir_hp.fits.gz");
      //anisoSkymapIR.write(fN + "aniso_ir_hp.fits.gz");
      //irRatio.write(fN + "ratio_ir_hp.fits.gz");

      //isoSkymapOpt.write(fN + "iso_opt_hp.fits.gz");
      //anisoSkymapOpt.write(fN + "aniso_opt_hp.fits.gz");
      //optRatio.write(fN + "ratio_opt_hp.fits.gz");

      //anisoRatioOpt.write(fN + "aniso_ratio_opt_inter_hp.fits.gz");
      //anisoRatioIR.write(fN + "aniso_ratio_ir_inter_hp.fits.gz");
      //anisoRatioCMB.write(fN + "aniso_ratio_cmb_inter_hp.fits.gz");

      //electrons.delete_array();
      delete anisoFunc;

    }
    
  } 

  if (3 == galdef.skymap_format || 4 == galdef.skymap_format) {

    //cout << "iPix: ";

#pragma omp parallel for schedule(dynamic) default(shared) 
    for (int iPix = 0; iPix < galaxy.IC_iso_hp_skymap[0].Npix(); ++iPix) {

      //cout << iPix << " ";

      SM::Coordinate co(galaxy.IC_iso_hp_skymap[0].pix2ang(iPix));
      const double l = co.l();
      const double b = co.b();
      vector< vector<double> > iso_IC;
            
      gen_IC_skymap_pixel(l, b, iso_IC, cameraLocation);

      for (int iComp = 0; iComp < nComps; ++iComp)
	for (int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

	  galaxy.IC_iso_hp_skymap[iComp][co][iEGamma] = iso_IC[iComp][iEGamma];// Already applied*galdef.ISRF_factors[iComp];

	  //cout << iComp << " " << iEGamma << " " << galaxy.IC_iso_hp_skymap[iComp][co][iEGamma];

	  //if(galdef.IC_anisotropic) galaxy.IC_aniso_hp_skymap[i_comp][co][iEgamma] = aniso_IC[i_comp][iEgamma]*galdef.ISRF_factors[i_comp];
	
	  //cout << endl;

	}
    }

    //cout << endl;

    if (galdef.IC_anisotropic && anisoValid) {

      for (int iPix = 0; iPix < galaxy.IC_iso_hp_skymap[0].Npix(); ++iPix) {

	for (int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

	  galaxy.IC_aniso_hp_skymap[0][iPix][iEGamma] = double(anisoRatioOpt[iPix][iEGamma])*galaxy.IC_iso_hp_skymap[0][iPix][iEGamma];
	  
	  galaxy.IC_aniso_hp_skymap[1][iPix][iEGamma] = double(anisoRatioIR[iPix][iEGamma])*galaxy.IC_iso_hp_skymap[1][iPix][iEGamma];
	  
	  galaxy.IC_aniso_hp_skymap[2][iPix][iEGamma] = double(anisoRatioCMB[iPix][iEGamma])*galaxy.IC_iso_hp_skymap[2][iPix][iEGamma];

	}

      }
    
    }
	
    for (int iComp = 0; iComp < nComps; ++iComp) {

      galaxy.IC_iso_hp_skymap[iComp].setSpectra(&galaxy.E_gamma[0], nEGammaBins);
      
      if (galdef.IC_anisotropic && anisoValid) 
	galaxy.IC_aniso_hp_skymap[iComp].setSpectra(&galaxy.E_gamma[0], nEGammaBins);

    }

    if (galdef.verbose >= 2) {

      for (unsigned int iComp = 0; iComp < nComps; ++iComp) {

	ostringstream isoBuf;
	isoBuf << "Isotropic inverse Compton skymap for ISRF component " << iComp;
	INFO(isoBuf.str());

	galaxy.IC_iso_hp_skymap[iComp].print(cout);
	
	if (galdef.IC_anisotropic && anisoValid) {

	  ostringstream anisoBuf;
	  anisoBuf << "Anisotropic inverse Compton skymap for component " << iComp;
	  INFO(anisoBuf.str());

	  galaxy.IC_aniso_hp_skymap[iComp].print(cout);
      
	}
      
      }

    }
      
  } else {

    cout << "iLong: ";

#pragma omp parallel for schedule(dynamic) default(shared)
    for (int iLong = 0; iLong < galaxy.n_long; ++iLong) {

      cout << iLong << " ";
      
      for (int iLat = 0; iLat < galaxy.n_lat; ++iLat) {

	const double l = galaxy.long_min + iLong*galaxy.d_long;
	const double b = galaxy.lat_min + iLat*galaxy.d_lat;
	vector< vector<double> > iso_IC;//, aniso_IC;
	
	//gen_IC_skymap_pixel(l, b, electrons, ielectrons, iso_IC, aniso_IC, Etarget, factor);
	
	gen_IC_skymap_pixel(l, b, iso_IC, cameraLocation);

	//Store and apply user defined factors to compontents

	for (int iComp = 0; iComp < nComps; ++iComp)
	  for (int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {
	    
	    galaxy.IC_iso_skymap[iComp].d2[iLong][iLat].s[iEGamma] = iso_IC[iComp][iEGamma]*galdef.ISRF_factors[iComp];

	    //if(galdef.IC_anisotropic) galaxy.IC_aniso_skymap[i_comp].d2[i_long][i_lat].s[iEgamma] = aniso_IC[i_comp][iEgamma]*galdef.ISRF_factors[i_comp];
	  }
    
      }//lat
    
    }//long

    cout << endl;

    if (galdef.IC_anisotropic && anisoValid) {

      for (int iLong = 0; iLong < galaxy.n_long; ++iLong) {
	
	for (int iLat = 0; iLat < galaxy.n_lat; ++iLat) {
	  
	  const double l = galaxy.long_min + iLong*galaxy.d_long;
	  const double b = galaxy.lat_min + iLat*galaxy.d_lat;
	  
	  SM::Coordinate coord(l, b);

	  const int iPix = anisoRatioOpt.ang2pix(coord.healpixAng());
	
	  for (int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {
	  
	    galaxy.IC_aniso_skymap[0].d2[iLong][iLat].s[iEGamma] = galaxy.IC_iso_skymap[0].d2[iLong][iLat].s[iEGamma]*double(anisoRatioOpt[iPix][iEGamma]);

	    galaxy.IC_aniso_skymap[1].d2[iLong][iLat].s[iEGamma] = galaxy.IC_iso_skymap[1].d2[iLong][iLat].s[iEGamma]*double(anisoRatioIR[iPix][iEGamma]);

	    galaxy.IC_aniso_skymap[2].d2[iLong][iLat].s[iEGamma] = galaxy.IC_iso_skymap[2].d2[iLong][iLat].s[iEGamma]*double(anisoRatioCMB[iPix][iEGamma]);

	  }

	}
    
      }

    }

    if (galdef.verbose >= 2) {

      for (int iComp = 0; iComp < nComps; ++iComp) {

	ostringstream isoBuf;
	isoBuf << "Isotropic inverse Compton skymap for ISRF component " << iComp;
	INFO(isoBuf.str());

	galaxy.IC_iso_skymap[iComp].print();
	
	if (galdef.IC_anisotropic && anisoValid) {

	  ostringstream anisoBuf;
	  anisoBuf << "Anisotropic inverse Compton skymap for component " << iComp;
	  INFO(anisoBuf.str());

	  galaxy.IC_aniso_skymap[iComp].print();
      
	}
      
      }
    
    } // galdef.verbose>=2
  
  }

  //if(galdef.IC_anisotropic==1) {

  //electrons.delete_array();  // IMOS20060420
  //}

  INFO("Exit: gen_IC_skymap");

  //exit(0);

  return stat;
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// calc. of a 2D array of g-ray emission for the given E_gammagrid for a
// particular pixel (l,b)   IMOS20080114
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::gen_IC_skymap_pixel(const double l, const double b, 
				 vector< vector<double> >& iso_IC,
				 const valarray<double>& cameraLocation) { 
  
  // Resizing and zeroing vectors

  iso_IC.resize(galaxy.n_ISRF_components);

  for (int i = 0; i < galaxy.n_ISRF_components; ++i) 
    iso_IC[i].resize(galaxy.n_E_gammagrid, 0);

  // Easier to modify location of camera in future ...

  const double cameraX = cameraLocation[0];
  const double cameraY = cameraLocation[1];
  const double cameraZ = cameraLocation[2];

  //cout << cameraX << " " << cameraY << " " << cameraZ << endl;

  const double cameraR = sqrt(cameraX*cameraX + cameraY*cameraY);
  
  const double lRad = l*utl::kConvertDegreesToRadians, bRad = b*utl::kConvertDegreesToRadians;

  const double sinb = sin(bRad), cosb = cos(bRad), sinl = sin(lRad), cosl = cos(lRad);
  
  // Integration along the line of sight
    
  const double ds = galdef.LoS_step/galdef.LoS_substep_number;

  // Test if the camera is outside the galaxy/CR region. If so, then we 
  // set the starting s value so that it is at the boundary of 
  // the galaxy/CR region, then the integration proceeds as before ...

  const bool outsideGalaxy2D = ( 2 == gcr[0].n_spatial_dimensions ) ?
    (cameraR > galaxy.r_max || 
     cameraZ < galaxy.z_min || cameraZ > galaxy.z_max)
    : false;

  const bool outsideGalaxy3D = ( 3 == gcr[0].n_spatial_dimensions ) ? 
    (cameraX < galaxy.x_min || cameraX > galaxy.x_max || 
     cameraY < galaxy.y_min || cameraY > galaxy.y_max || 
     cameraZ < galaxy.z_min || cameraZ > galaxy.z_max) 
    : false;

  double s =
    (outsideGalaxy2D ? 
     ClosestDistance2D(l, b, galaxy.r_max, galaxy.z_min, galaxy.z_max, cameraLocation) : 
     (outsideGalaxy3D ?
      ClosestDistance3D(l, b, galaxy.x_min, galaxy.x_max, galaxy.y_min, galaxy.y_max, galaxy.z_min, galaxy.z_max, cameraLocation) : 0));

  bool complete = false;

  while (1) {

    s += ds;
	  
    const double dx = s*cosb*cosl;
    const double dy = s*cosb*sinl;
    const double dz = s*sinb;
    
    const double x = cameraX - dx;
    const double y = cameraY - dy;
    const double z = cameraZ + dz;
    
    const double r = sqrt(x*x + y*y);

    // These are the loop break criteria. Much simpler encoding of 
    // the break condition than in the earlier pixel routines.

    if (2 == gcr[0].n_spatial_dimensions) {
      
      complete = (r > galaxy.r_max || 
		  z < galaxy.z_min || z > galaxy.z_max);
      
    }
    
    if (3 == gcr[0].n_spatial_dimensions) {
      
      complete = (x < galaxy.x_min || x > galaxy.x_max || 
		  y < galaxy.y_min || y > galaxy.y_max || 
		  z < galaxy.z_min || z > galaxy.z_max);
      
    }
    
    if (complete)
      break;

    if (2 == gcr[0].n_spatial_dimensions) {
      
      const double dR = (r - galaxy.r_min)/(galaxy.r_max - galaxy.r_min);
      const double dZ = (z - galaxy.z_min)/(galaxy.z_max - galaxy.z_min);
      
      const unsigned int iR = (dR < 1. ? dR*(galaxy.n_rgrid-1) : galaxy.n_rgrid-1);
      const unsigned int iZ = (dZ < 1. ? dZ*(galaxy.n_zgrid-1) : galaxy.n_zgrid-1);
      
      for (int i_comp = 0; i_comp < galaxy.n_ISRF_components; ++i_comp) {//IMOS20060420
	
	for (int iEgamma = 0; iEgamma < galaxy.n_E_gammagrid; ++iEgamma) {//IMOS20060420

	  double delta, xx[8][3], f[8], yy[7];
	  	
	  if (iR == galaxy.n_rgrid-1 || 
	      iZ == galaxy.n_zgrid-1 || 
	      galdef.verbose == -457) { //-457 -old method
	    
	    delta = ds*kpc2cm*galaxy.IC_iso_emiss[i_comp].d2[iR][iZ].s[iEgamma];
	
	  } else {

	    xx[0][0] = galaxy.r[iR]; 
	    xx[0][1] = galaxy.z[iZ + 1];  
	    f[0] = galaxy.IC_iso_emiss[i_comp].d2[iR][iZ + 1].s[iEgamma];
	    
	    xx[1][0] = galaxy.r[iR + 1]; 
	    xx[1][1] = galaxy.z[iZ + 1];  
	    f[1] = galaxy.IC_iso_emiss[i_comp].d2[iR + 1][iZ + 1].s[iEgamma];
	    
	    xx[2][0] = galaxy.r[iR]; 
	    xx[2][1] = galaxy.z[iZ];  
	    f[2] = galaxy.IC_iso_emiss[i_comp].d2[iR][iZ].s[iEgamma];
	    
	    xx[3][0] = galaxy.r[iR + 1]; 
	    xx[3][1] = galaxy.z[iZ];  
	    f[3] = galaxy.IC_iso_emiss[i_comp].d2[iR + 1][iZ].s[iEgamma];
	    
	    yy[0] = (f[0] - f[1])/(xx[0][0] - xx[1][0])*(r - xx[0][0]) + f[0]; // interpolation in R
	    yy[1] = (f[2] - f[3])/(xx[2][0] - xx[3][0])*(r - xx[2][0]) + f[2];
	    
	    yy[2] = (yy[0] - yy[1])/(xx[0][1] - xx[2][1])*(z - xx[0][1]) + yy[0];   // interpolation in z
	    
	    delta = ds*kpc2cm*yy[2];

	  }

	  iso_IC[i_comp][iEgamma] += delta;

	}

      }

    }

    if (3 == gcr[0].n_spatial_dimensions) {
      
      const double dX = (x - galaxy.x_min)/(galaxy.x_max - galaxy.x_min);
      const double dY = (y - galaxy.y_min)/(galaxy.y_max - galaxy.y_min);
      const double dZ = (z - galaxy.z_min)/(galaxy.z_max - galaxy.z_min);
      
      const unsigned int iX = (dX < 1. ? dX*(galaxy.n_xgrid-1) : galaxy.n_xgrid-1);
      const unsigned int iY = (dY < 1. ? dY*(galaxy.n_ygrid-1) : galaxy.n_ygrid-1);
      const unsigned int iZ = (dZ < 1. ? dZ*(galaxy.n_zgrid-1) : galaxy.n_zgrid-1);

      for (int i_comp = 0; i_comp < galaxy.n_ISRF_components; ++i_comp) {//IMOS20060420
	
	for (int iEgamma = 0; iEgamma < galaxy.n_E_gammagrid; ++iEgamma) {//IMOS20060420

	  double delta, xx[8][3], f[8], yy[7];

	  if (iX == galaxy.n_xgrid-1 || 
	      iY == galaxy.n_ygrid-1 || 
	      iZ == galaxy.n_zgrid-1) {

	    delta = ds*kpc2cm*galaxy.IC_iso_emiss[i_comp].d3[iX][iY][iZ].s[iEgamma];
	  
	  } else {  // linear interpolation
	    
	    //  x[0]=(x0,z1,y0), x[1]=(x1,z1,y0);   x[4]=(x0,z1,y1), x[5]=(x1,z1,y1);
	    //  x[2]=(x0,z0,y0), x[3]=(x1,z0,y0);   x[6]=(x0,z0,y1), x[7]=(x1,z0,y1);
	    xx[0][0] = galaxy.x[iX]; 
	    xx[0][1] = galaxy.z[iZ + 1]; 
	    xx[0][2] = galaxy.y[iY];  
	    f[0] = galaxy.IC_iso_emiss[i_comp].d3[iX][iY][iZ + 1].s[iEgamma];

	    xx[1][0] = galaxy.x[iX + 1]; 
	    xx[1][1] = galaxy.z[iZ + 1]; 
	    xx[1][2] = galaxy.y[iY];  
	    f[1] = galaxy.IC_iso_emiss[i_comp].d3[iX + 1][iY][iZ + 1].s[iEgamma];

	    xx[2][0] = galaxy.x[iX]; 
	    xx[2][1] = galaxy.z[iZ]; 
	    xx[2][2] = galaxy.y[iY];  
	    f[2] = galaxy.IC_iso_emiss[i_comp].d3[iX][iY][iZ].s[iEgamma];

	    xx[3][0] = galaxy.x[iX + 1]; 
	    xx[3][1] = galaxy.z[iZ]; 
	    xx[3][2] = galaxy.y[iY];  
	    f[3] = galaxy.IC_iso_emiss[i_comp].d3[iX + 1][iY][iZ].s[iEgamma];

	    xx[4][0] = galaxy.x[iX]; 
	    xx[4][1] = galaxy.z[iZ + 1]; 
	    xx[4][2] = galaxy.y[iY + 1];  
	    f[4] = galaxy.IC_iso_emiss[i_comp].d3[iX][iY + 1][iZ + 1].s[iEgamma];

	    xx[5][0] = galaxy.x[iX + 1]; 
	    xx[5][1] = galaxy.z[iZ + 1]; 
	    xx[5][2] = galaxy.y[iY + 1];  
	    f[5] = galaxy.IC_iso_emiss[i_comp].d3[iX + 1][iY + 1][iZ + 1].s[iEgamma];

	    xx[6][0] = galaxy.x[iX]; 
	    xx[6][1] = galaxy.z[iZ]; 
	    xx[6][2] = galaxy.y[iY + 1];  
	    f[6] = galaxy.IC_iso_emiss[i_comp].d3[iX][iY + 1][iZ].s[iEgamma];

	    xx[7][0] = galaxy.x[iX + 1]; 
	    xx[7][1] = galaxy.z[iZ]; 
	    xx[7][2] = galaxy.y[iY + 1];  
	    f[7] = galaxy.IC_iso_emiss[i_comp].d3[iX + 1][iY + 1][iZ].s[iEgamma];
	    
	    //f[0] = f[1] = f[2] = f[3] = f[4] = f[5] = f[6] = f[7] = 1;

	    yy[0] = (f[0] - f[1])/(xx[0][0] - xx[1][0])*(x - xx[0][0]) + f[0]; // interpolation in x
	    yy[1] = (f[2] - f[3])/(xx[2][0] - xx[3][0])*(x - xx[2][0]) + f[2];
	    yy[2] = (f[4] - f[5])/(xx[4][0] - xx[5][0])*(x - xx[4][0]) + f[4];
	    yy[3] = (f[6] - f[7])/(xx[6][0] - xx[7][0])*(x - xx[6][0]) + f[6];
	    
	    yy[4] = (yy[0] - yy[1])/(xx[0][1] - xx[2][1])*(z - xx[0][1]) + yy[0];   // interpolation in z
	    yy[5] = (yy[2] - yy[3])/(xx[4][1] - xx[6][1])*(z - xx[4][1]) + yy[2];
	    
	    yy[6] = (yy[4] - yy[5])/(xx[0][2] - xx[4][2])*(y - xx[0][2]) + yy[4];   // interpolation in y
	    	    
	    delta = ds*kpc2cm*yy[6];

	  }

	  iso_IC[i_comp][iEgamma] += delta;

	}

      }
   
    }

  } 
  
  return 0;

}
