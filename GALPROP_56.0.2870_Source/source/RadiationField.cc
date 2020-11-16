#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <dirent.h>

#include <PhysicalConstants.h>
#include <ErrorLogger.h>
#include <StatusIndicator.h>

#include <RadiationField.h>

#include <CCfits/CCfits>

//#include <CLHEP/Vector/ThreeVector.h>

//using namespace std;
using namespace rf;
//using namespace utl;
//using namespace CLHEP;

bool DirectoryExists(const char* path) {

  if (path == nullptr)
    return false;

  bool exists = false;

  DIR* pDir = opendir(path);

  if (pDir != nullptr) {

    exists = true;
    (void) closedir(pDir);

  }

  return exists;

}

RadiationField::RadiationField() : 
  fRebinnedSkymapOrder(1), 
  fNumberOfComponents(0), 
  fLuminosity(0), 
  fDustMass(0), 
  fCacheBuilt(false),
  fUse2D(false), fUse3D(false) {

}

RadiationField::RadiationField(const std::string& filename,
			       const std::valarray<double>& freq,
			       int rebinnedSkymapOrder) :
  fRebinnedSkymapOrder(rebinnedSkymapOrder),
  fNumberOfComponents(0),
  fLuminosity(0),
  fDustMass(0),
  fCacheBuilt(false),
  fUse2D(false), fUse3D(false) {

  fFrequency.resize(freq.size());
  fFrequency = freq;
  
  const std::string prefix = filename.substr(0, filename.find_last_of("/"));
  
  //const std::string dimensionality = prefix.substr(prefix.find_last_of("/")+1)

  if (DirectoryExists(prefix.c_str())) {

    if (filename.empty()) {

      INFO("Empty radiation field filename passed to ReadRadiationField.");
      exit(-1);

    }
  
    std::ostringstream buf;
    buf << "Reading from " << filename;
    INFO(buf.str());

    const std::string prefix = filename.substr(0, filename.find_last_of("/")+1);
    
    //cout << prefix << endl;
    
    fPrefix = prefix;
    
    std::ifstream rfFile(filename.c_str());

    // Read first line of ISRF file -- if it is `!!3D' then we're reading a 
    // 3D ISRF run. Otherwise we go with 2D ...

    std::string line;
    std::getline(rfFile, line);
    rfFile.close();

    auto pos = line.find_first_of("!!3D");

    const std::string dimensionality = (pos > line.size() ? "2D" : "3D");

    std::cout << prefix << " " << dimensionality << " " << filename << " " << DirectoryExists(prefix.c_str()) << std::endl;

    //exit(0);

    if (dimensionality == "3D") {
      
      ReadRadiationField3D(filename);
      fUse3D = true;
      
    } else if (dimensionality == "2D") {
      
      ReadRadiationField2D(filename);
      fUse2D = true;
      
    } else {
      
      std::ostringstream buf;
      buf.str("");
      buf << "Error reading from file " << filename << ": cannot find dimensionality";
      ERROR(buf.str());
      throw(std::invalid_argument(buf.str()));
      
    }

  } else {

    std::ostringstream buf;
    buf.str("");
    buf << "Error reading from file " << filename << ": directory does not exist";
    ERROR(buf.str());
    throw(std::invalid_argument(buf.str()));
      
  }

}

RadiationField::~RadiationField() {

  ClearData();

}

const Skymap<double>
RadiationField::GetSkymap(const ThreeVector& pos,
			  const STELLARCOMPONENT component,
			  const int healpixOrder) {

  assert (component >= TOTAL && component <= INFRARED);

  if (fRebinnedSkymapOrder == healpixOrder)
    return GetSkymap(pos, component);
  else {

    //FlushSkymapCache(component);

    //const int rebinnedSkymapOrder = fRebinnedSkymapOrder;

    //fRebinnedSkymapOrder = healpixOrder;

    Skymap<double> skymapTmp = GetSkymap(pos, component);

    Skymap<double> skymap = skymapTmp.rebin(healpixOrder);

    //FlushSkymapCache(component);

    //fRebinnedSkymapOrder = rebinnedSkymapOrder;

    return skymap;

  }

}

const Skymap<double> 
RadiationField::GetSkymap(const ThreeVector& pos, 
			  const STELLARCOMPONENT component) {

  assert(fUse2D || fUse3D);

  return (fUse2D ? GetSkymap2D(pos, component) : GetSkymap3D(pos, component));

  //cout << "Cache built: " << fCacheBuilt[component] << endl;

  /*if (!fCacheBuilt)//[component])
    BuildSkymapCache2D();//component);

  const double r = sqrt(pos.x()*pos.x() + pos.y()*pos.y()), z = pos.z(), absZ = fabs(z);

  if (r > fRData[fRData.size()-1] || absZ > fZData[fZData.size()-1]) {

    const Skymap<double>& skymap = *fSkymapOrderedData2D[0][0][component];

    return ((Skymap<double>() = skymap) = 0); 

  }

  // Bilinear interpolation to form skymap

  const double* rVal = std::lower_bound(&fRData[0], &fRData[fRData.size()-1], r);

  const int rIndex = rVal - &fRData[0];

  double rCoeff = 0, dR = 0;

  if (r <= fRData[rIndex] && r >= fRData[rIndex])
    rCoeff = 1;
  else if (r < fRData[fRData.size()-1]) {
    
    dR = fRData[rIndex] - fRData[rIndex-1];
    
    rCoeff = (fRData[rIndex] - r)/dR;
    
  }

  const double* zVal = std::lower_bound(&fZData[0], &fZData[fZData.size()-1], absZ);

  const int zIndex = zVal - &fZData[0];
 
  double zCoeff = 0, dZ = 0;

  if (absZ <= fZData[zIndex] && absZ >= fZData[zIndex]) 
    zCoeff = 1;
  else if (absZ < fZData[fZData.size()-1]) {
    
    dZ = fZData[zIndex] - fZData[zIndex-1];
    
    zCoeff = (fZData[zIndex] - absZ)/dZ;
    
  }
  
  // Form the bilinear interpolated skymap 

  Skymap<double> skymap;
  
  if (rCoeff >= 1. && zCoeff >= 1.) {

    if (z >= 0.)
      return *fSkymapOrderedData2D[rIndex][zIndex][component];
    else {

      const Skymap<double>& skymap1 = *fSkymapOrderedData2D[rIndex][zIndex][component];

      Skymap<double> skymap1Mirror = skymap1;

      for (unsigned long i = 0; i < skymap1.Npix(); ++i) {

	SM::Coordinate coord = skymap1.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap1Mirror[mirrorCoord] = skymap1[coord];

      }

      return skymap1Mirror;

    }

  } else if (rCoeff >= 1.) {

    const Skymap<double>& skymap1 = *fSkymapOrderedData2D[rIndex][zIndex-1][component];
    
    const Skymap<double>& skymap2 = *fSkymapOrderedData2D[rIndex][zIndex][component];

    if (z >= 0) 
      skymap = skymap1*zCoeff + skymap2*(1. - zCoeff);
    else {

      Skymap<double> skymap1Mirror = skymap1;
      Skymap<double> skymap2Mirror = skymap2;

      for (unsigned long i = 0; i < skymap1.Npix(); ++i) {

	SM::Coordinate coord = skymap1.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap1Mirror[mirrorCoord] = skymap1[coord];

      }

      for (unsigned long i = 0; i < skymap2.Npix(); ++i) {

	SM::Coordinate coord = skymap2.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap2Mirror[mirrorCoord] = skymap2[coord];

      }

      skymap = skymap1Mirror*zCoeff + skymap2Mirror*(1. - zCoeff); 

    }

  } else if (zCoeff >= 1.) {

    const Skymap<double>& skymap1 = *fSkymapOrderedData2D[rIndex-1][zIndex][component];

    const Skymap<double>& skymap2 = *fSkymapOrderedData2D[rIndex][zIndex][component];

    if (z >= 0.)
      skymap = skymap1*rCoeff + skymap2*(1. - rCoeff);
    else {

      Skymap<double> skymap1Mirror = skymap1;
      Skymap<double> skymap2Mirror = skymap2;

      for (unsigned long i = 0; i < skymap1.Npix(); ++i) {

	SM::Coordinate coord = skymap1.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap1Mirror[mirrorCoord] = skymap1[coord];

      }

      for (unsigned long i = 0; i < skymap2.Npix(); ++i) {

	SM::Coordinate coord = skymap2.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap2Mirror[mirrorCoord] = skymap2[coord];

      }

      skymap = skymap1Mirror*rCoeff + skymap2Mirror*(1. - rCoeff);

    }

  } else {

    const Skymap<double>& skymap1 = *fSkymapOrderedData2D[rIndex-1][zIndex-1][component];

    const Skymap<double>& skymap2 = *fSkymapOrderedData2D[rIndex-1][zIndex][component];

    const Skymap<double>& skymap3 = *fSkymapOrderedData2D[rIndex][zIndex-1][component];

    const Skymap<double>& skymap4 = *fSkymapOrderedData2D[rIndex][zIndex][component];

    if (z >= 0.) 
      skymap = skymap1*rCoeff*zCoeff + skymap2*(1. - zCoeff)*rCoeff + skymap3*zCoeff*(1. - rCoeff) + skymap4*(1. - rCoeff)*(1. - zCoeff);
    else {

      Skymap<double> skymap1Mirror = skymap1;
      Skymap<double> skymap2Mirror = skymap2;
      Skymap<double> skymap3Mirror = skymap3;
      Skymap<double> skymap4Mirror = skymap4;

      for (unsigned long i = 0; i < skymap1.Npix(); ++i) {

	SM::Coordinate coord = skymap1.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap1Mirror[mirrorCoord] = skymap1[coord];

      }

      for (unsigned long i = 0; i < skymap2.Npix(); ++i) {

	SM::Coordinate coord = skymap2.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap2Mirror[mirrorCoord] = skymap2[coord];

      }

      for (unsigned long i = 0; i < skymap3.Npix(); ++i) {

	SM::Coordinate coord = skymap3.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap3Mirror[mirrorCoord] = skymap3[coord];

      }

      for (unsigned long i = 0; i < skymap4.Npix(); ++i) {

	SM::Coordinate coord = skymap4.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap4Mirror[mirrorCoord] = skymap4[coord];

      }
      
      skymap = skymap1Mirror*rCoeff*zCoeff + skymap2Mirror*(1. - zCoeff)*rCoeff + skymap3Mirror*zCoeff*(1. - rCoeff) + skymap4Mirror*(1. - rCoeff)*(1. - zCoeff); 

    }

  }

  return skymap;
  */
}

const Skymap<double> 
RadiationField::GetSkymap2D(const ThreeVector& pos, 
			    const STELLARCOMPONENT component) {

  //cout << "Cache built: " << fCacheBuilt[component] << endl;

  if (!fCacheBuilt)//[component])
    BuildSkymapCache2D();//component);

  const double r = sqrt(pos.x()*pos.x() + pos.y()*pos.y()), z = pos.z(), absZ = fabs(z);

  if (r > fRData[fRData.size()-1] || absZ > fZData[fZData.size()-1]) {

    const Skymap<double>& skymap = *fSkymapOrderedData2D[0][0][component];

    return ((Skymap<double>() = skymap) = 0); 

  }

  // Bilinear interpolation to form skymap

  const double* rVal = std::lower_bound(&fRData[0], &fRData[fRData.size()-1], r);

  const int rIndex = rVal - &fRData[0];

  double rCoeff = 0, dR = 0;

  if (r <= fRData[rIndex] && r >= fRData[rIndex])
    rCoeff = 1;
  else if (r < fRData[fRData.size()-1]) {
    
    dR = fRData[rIndex] - fRData[rIndex-1];
    
    rCoeff = (fRData[rIndex] - r)/dR;
    
  }

  const double* zVal = std::lower_bound(&fZData[0], &fZData[fZData.size()-1], absZ);

  const int zIndex = zVal - &fZData[0];
 
  double zCoeff = 0, dZ = 0;

  if (absZ <= fZData[zIndex] && absZ >= fZData[zIndex]) 
    zCoeff = 1;
  else if (absZ < fZData[fZData.size()-1]) {
    
    dZ = fZData[zIndex] - fZData[zIndex-1];
    
    zCoeff = (fZData[zIndex] - absZ)/dZ;
    
  }
  
  // Form the bilinear interpolated skymap 

  Skymap<double> skymap;
  
  if (rCoeff >= 1. && zCoeff >= 1.) {

    if (z >= 0.)
      return *fSkymapOrderedData2D[rIndex][zIndex][component];
    else {

      const Skymap<double>& skymap1 = *fSkymapOrderedData2D[rIndex][zIndex][component];

      Skymap<double> skymap1Mirror = skymap1;

      for (unsigned long i = 0; i < skymap1.Npix(); ++i) {

	SM::Coordinate coord = skymap1.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap1Mirror[mirrorCoord] = skymap1[coord];

      }

      return skymap1Mirror;

    }

  } else if (rCoeff >= 1.) {

    const Skymap<double>& skymap1 = *fSkymapOrderedData2D[rIndex][zIndex-1][component];
    
    const Skymap<double>& skymap2 = *fSkymapOrderedData2D[rIndex][zIndex][component];

    if (z >= 0) 
      skymap = skymap1*zCoeff + skymap2*(1. - zCoeff);
    else {

      Skymap<double> skymap1Mirror = skymap1;
      Skymap<double> skymap2Mirror = skymap2;

      for (unsigned long i = 0; i < skymap1.Npix(); ++i) {

	SM::Coordinate coord = skymap1.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap1Mirror[mirrorCoord] = skymap1[coord];

      }

      for (unsigned long i = 0; i < skymap2.Npix(); ++i) {

	SM::Coordinate coord = skymap2.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap2Mirror[mirrorCoord] = skymap2[coord];

      }

      skymap = skymap1Mirror*zCoeff + skymap2Mirror*(1. - zCoeff); 

    }

  } else if (zCoeff >= 1.) {

    const Skymap<double>& skymap1 = *fSkymapOrderedData2D[rIndex-1][zIndex][component];

    const Skymap<double>& skymap2 = *fSkymapOrderedData2D[rIndex][zIndex][component];

    if (z >= 0.)
      skymap = skymap1*rCoeff + skymap2*(1. - rCoeff);
    else {

      Skymap<double> skymap1Mirror = skymap1;
      Skymap<double> skymap2Mirror = skymap2;

      for (unsigned long i = 0; i < skymap1.Npix(); ++i) {

	SM::Coordinate coord = skymap1.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap1Mirror[mirrorCoord] = skymap1[coord];

      }

      for (unsigned long i = 0; i < skymap2.Npix(); ++i) {

	SM::Coordinate coord = skymap2.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap2Mirror[mirrorCoord] = skymap2[coord];

      }

      skymap = skymap1Mirror*rCoeff + skymap2Mirror*(1. - rCoeff);

    }

  } else {

    const Skymap<double>& skymap1 = *fSkymapOrderedData2D[rIndex-1][zIndex-1][component];

    const Skymap<double>& skymap2 = *fSkymapOrderedData2D[rIndex-1][zIndex][component];

    const Skymap<double>& skymap3 = *fSkymapOrderedData2D[rIndex][zIndex-1][component];

    const Skymap<double>& skymap4 = *fSkymapOrderedData2D[rIndex][zIndex][component];

    if (z >= 0.) 
      skymap = skymap1*rCoeff*zCoeff + skymap2*(1. - zCoeff)*rCoeff + skymap3*zCoeff*(1. - rCoeff) + skymap4*(1. - rCoeff)*(1. - zCoeff);
    else {

      Skymap<double> skymap1Mirror = skymap1;
      Skymap<double> skymap2Mirror = skymap2;
      Skymap<double> skymap3Mirror = skymap3;
      Skymap<double> skymap4Mirror = skymap4;

      for (unsigned long i = 0; i < skymap1.Npix(); ++i) {

	SM::Coordinate coord = skymap1.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap1Mirror[mirrorCoord] = skymap1[coord];

      }

      for (unsigned long i = 0; i < skymap2.Npix(); ++i) {

	SM::Coordinate coord = skymap2.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap2Mirror[mirrorCoord] = skymap2[coord];

      }

      for (unsigned long i = 0; i < skymap3.Npix(); ++i) {

	SM::Coordinate coord = skymap3.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap3Mirror[mirrorCoord] = skymap3[coord];

      }

      for (unsigned long i = 0; i < skymap4.Npix(); ++i) {

	SM::Coordinate coord = skymap4.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap4Mirror[mirrorCoord] = skymap4[coord];

      }
      
      skymap = skymap1Mirror*rCoeff*zCoeff + skymap2Mirror*(1. - zCoeff)*rCoeff + skymap3Mirror*zCoeff*(1. - rCoeff) + skymap4Mirror*(1. - rCoeff)*(1. - zCoeff); 

    }

  }

  assert(skymap.Order() ==  (*fSkymapOrderedData2D[0][0][component]).Order());
  
  return skymap;

}

const Skymap<double> 
RadiationField::GetSkymap3D(const ThreeVector& pos, 
			    const STELLARCOMPONENT component) {

  //cout << "Cache built: " << fCacheBuilt[component] << endl;

  if (!fCacheBuilt)//[component])
    BuildSkymapCache3D();//component);

  auto cidx(-1);

  if (component == RadiationField::TOTAL)
    cidx = 0;

  if (component == RadiationField::OPTICAL)
    cidx = 1; // Explicit mapping because we only store total, optical, infrared for 3D models

  if (component == RadiationField::INFRARED)
    cidx = 2; // Explicit mapping because we only store total, optical, infrared for 3D models

  assert(cidx > -1 && cidx < 3);

  auto r = sqrt(pos.x()*pos.x() + pos.y()*pos.y()), phi = atan2(pos.y(), pos.x()), phiDeg = (phi < 0 ? phi + 2.*utl::kPi : (phi > 2.*utl::kPi ? phi - 2.*utl::kPi : phi))*180./utl::kPi, z = pos.z();

  // Want to check for this early because we rely on the corrected angle and it should be corrected!
  assert(phiDeg >= 0.);
  assert(phiDeg <= 360.); 
  
  //std::cout << pos.x() << " " << pos.y() << " " << r << " " << phiDeg << " " << z << std::endl;

  if (r > fRData[fRData.size()-1] || z > fZData[fZData.size()-1] || z < fZData[0]) {

    const auto& skymap = *fSkymapData3D[0][cidx];

    return ((Skymap<double>() = skymap) = 0); 
    
  }
  
  // Deal with the case where the rarray stops short of R = 0. We assume that 
  // the supplied radiation field from R = 0 to fRData[0] is constant and 
  // taken as the value at fRData[0], so just set r to fRData[0] ...

  if (r >= 0. && r <= fRData[0])
    r = fRData[0];

  // Trilinear interpolation 

  Skymap<double> skymap;
  
  const auto ridx = std::upper_bound(&fRData[0], &fRData[fRData.size()-1], r) - &fRData[0];
  
  const auto dr = (fRData[ridx] - fRData[ridx-1]);
  assert(dr > 0.); // We want to fail otherwise

  const auto rcoeff = (fRData[ridx] - r)/dr;
  
  const auto phiidx = std::upper_bound(&fPhiData[0], &fPhiData[fPhiData.size()-1], phiDeg) - &fPhiData[0];
  
  const auto dphi = (phiDeg < fPhiData[fPhiData.size()-1]) ? (fPhiData[phiidx] - fPhiData[phiidx-1]) : (360. - fPhiData[fPhiData.size()-1]);
  assert(dphi > 0.); // We want to fail otherwise

  const auto phicoeff = (phiDeg < fPhiData[fPhiData.size()-1] && dphi > 0.) ? (fPhiData[phiidx] - phiDeg)/dphi : (360. - phiDeg)/dphi;
  
  const auto zidx = std::upper_bound(&fZData[0], &fZData[fZData.size()-1], z) - &fZData[0];
  
  const auto dz = (fZData[zidx] - fZData[zidx-1]);
  assert(dz > 0.); // We want to fail otherwise

  const auto zcoeff = (fZData[zidx] - z)/dz;
 
  if (phi > fPhiData[fPhiData.size()-1]) { // This deals with the wrap around in phi coordinate between the last bin in the grid and 360./0. deg

    const size_t idx1 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (phiidx-1);

    const auto& skymap1 = *fSkymapData3D[idx1][cidx];
     
    const size_t idx2 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (0); // This is explicit to remind how the indexing is done

    const auto& skymap2 = *fSkymapData3D[idx2][cidx];

    const size_t idx3 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (phiidx-1);

    const auto& skymap3 = *fSkymapData3D[idx3][cidx];

    const size_t idx4 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (0);

    const auto& skymap4 = *fSkymapData3D[idx4][cidx];
          
    const size_t idx5 = (zidx)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (phiidx-1);

    const auto& skymap5 = *fSkymapData3D[idx5][cidx];
    
    const size_t idx6 = (zidx)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (0);

    const auto& skymap6 = *fSkymapData3D[idx6][cidx];
    
    const size_t idx7 = (zidx)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (phiidx-1);

    const auto& skymap7 = *fSkymapData3D[idx7][cidx];
    
    const size_t idx8 = (zidx)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (0);

    const auto& skymap8 = *fSkymapData3D[idx8][cidx];
        
    skymap = skymap1*zcoeff*rcoeff*phicoeff + 
      skymap2*zcoeff*rcoeff*(1. - phicoeff) + 
      skymap3*zcoeff*(1. - rcoeff)*phicoeff + 
      skymap4*zcoeff*(1. - rcoeff)*(1. - phicoeff) + 
      skymap5*(1. - zcoeff)*rcoeff*phicoeff + 
      skymap6*(1. - zcoeff)*rcoeff*(1. - phicoeff) +
      skymap7*(1. - zcoeff)*(1. - rcoeff)*phicoeff + 
      skymap8*(1. - zcoeff)*(1. - rcoeff)*(1. - phicoeff);

  } else {

    const size_t idx1 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (phiidx-1);
    
    const auto& skymap1 = *fSkymapData3D[idx1][cidx];
 
    const size_t idx2 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (phiidx);
    
    const auto& skymap2 = *fSkymapData3D[idx2][cidx];
 
    const size_t idx3 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (phiidx-1);
    
    const auto& skymap3 = *fSkymapData3D[idx3][cidx];
 
    const size_t idx4 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (phiidx);
    
    const auto& skymap4 = *fSkymapData3D[idx4][cidx];
          
    const size_t idx5 = (zidx)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (phiidx-1);
    
    const auto& skymap5 = *fSkymapData3D[idx5][cidx];
    
    const size_t idx6 = (zidx)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (phiidx);
    
    const auto& skymap6 = *fSkymapData3D[idx6][cidx];
    
    const size_t idx7 = (zidx)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (phiidx-1);
    
    const auto& skymap7 = *fSkymapData3D[idx7][cidx];
    
    const size_t idx8 = (zidx)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (phiidx);
    
    const auto& skymap8 = *fSkymapData3D[idx8][cidx];
        
    skymap = skymap1*zcoeff*rcoeff*phicoeff + 
      skymap2*zcoeff*rcoeff*(1. - phicoeff) + 
      skymap3*zcoeff*(1. - rcoeff)*phicoeff + 
      skymap4*zcoeff*(1. - rcoeff)*(1. - phicoeff) + 
      skymap5*(1. - zcoeff)*rcoeff*phicoeff + 
      skymap6*(1. - zcoeff)*rcoeff*(1. - phicoeff) +
      skymap7*(1. - zcoeff)*(1. - rcoeff)*phicoeff + 
      skymap8*(1. - zcoeff)*(1. - rcoeff)*(1. - phicoeff);

  }
  
  assert(skymap.Order() == (*fSkymapData3D[0][cidx]).Order());
  
  return skymap;

}

const std::valarray<double>
RadiationField::GetNumberDensity(const ThreeVector& pos, 
				 const STELLARCOMPONENT component) {

  assert(fUse2D || fUse3D);

  return (fUse2D ? GetNumberDensity2D(pos, component) : GetNumberDensity3D(pos, component));

}

const std::valarray<double>
RadiationField::GetNumberDensity2D(const ThreeVector& pos, 
				   const STELLARCOMPONENT component) {

  std::valarray<double> spec(0., fFrequency.size());

  const double r = sqrt(pos.x()*pos.x() + pos.y()*pos.y()), z = pos.z(), absZ = fabs(z);
  
  if (r > fRData[fRData.size()-1] || z > fZData[fZData.size()-1]) {

    return spec;

  }

  // Bilinear interpolation 

  const double* rVal = std::lower_bound(&fRData[0], &fRData[fRData.size()-1], r);

  const int rIndex = rVal - &fRData[0];

  double rCoeff = 0, dR = 0;

  if (r <= fRData[rIndex] && r >= fRData[rIndex])
    rCoeff = 1;
  else if (r < fRData[fRData.size()-1]) {
    
    dR = fRData[rIndex] - fRData[rIndex-1];
    
    rCoeff = (fRData[rIndex] - r)/dR;
    
  }

  const double* zVal = std::lower_bound(&fZData[0], &fZData[fZData.size()-1], absZ);

  const int zIndex = zVal - &fZData[0];
 
  double zCoeff = 0, dZ = 0;

  if (absZ <= fZData[zIndex] && absZ >= fZData[zIndex]) 
    zCoeff = 1;
  else if (absZ < fZData[fZData.size()-1]) {
    
    dZ = fZData[zIndex] - fZData[zIndex-1];
    
    zCoeff = (fZData[zIndex] - absZ)/dZ;
    
  }
  
  //cout << z << " " << zIndex << " " << dZ << " " << zCoeff << " " << fZData[zIndex] << endl;

  //exit(0);

  auto& wl = fWavelength;

  std::valarray<double> freq(0., wl.size());
      
  for (size_t iWl = 0; iWl < wl.size(); ++iWl)
    freq[iWl] = utl::kSpeedOfLight_SI*1./(wl[wl.size()-1-iWl]*utl::micron/utl::m);
  
  std::valarray<double> energyDensity(0., wl.size());  

  if (rCoeff >= 1. && zCoeff >= 1.) {

    energyDensity = *fEnergyDensity2D[rIndex][zIndex][component];

  } else if (rCoeff >= 1.) {

    auto& energyDensity1 = *fEnergyDensity2D[rIndex][zIndex-1][component];
    
    auto& energyDensity2 = *fEnergyDensity2D[rIndex][zIndex][component];

    energyDensity = energyDensity1*zCoeff + energyDensity2*(1. - zCoeff);

  } else if (zCoeff >= 1.) {

    auto& energyDensity1 = *fEnergyDensity2D[rIndex-1][zIndex][component];

    auto& energyDensity2 = *fEnergyDensity2D[rIndex][zIndex][component];

    energyDensity = energyDensity1*rCoeff + energyDensity2*(1. - rCoeff);

  } else {

    auto& energyDensity1 = *fEnergyDensity2D[rIndex-1][zIndex-1][component];

    auto& energyDensity2 = *fEnergyDensity2D[rIndex-1][zIndex][component];

    auto& energyDensity3 = *fEnergyDensity2D[rIndex][zIndex-1][component];

    auto& energyDensity4 = *fEnergyDensity2D[rIndex][zIndex][component];

    energyDensity = energyDensity1*rCoeff*zCoeff + energyDensity2*(1. - zCoeff)*rCoeff + energyDensity3*zCoeff*(1. - rCoeff) + energyDensity4*(1. - rCoeff)*(1. - zCoeff);

  }

  if (fFrequency.size() < freq.size()) { 

    // For binning of supplied ISRF finer than required
    
    std::valarray<double> binCount(0., fFrequency.size()), energy(fFrequency.size());
    
    energy = utl::kPlanck_SI/utl::e_SI*fFrequency;
    
    for (auto iFreq = 0; iFreq < freq.size(); ++iFreq) {
      
      auto index = (log10(freq[iFreq]) - log10(fFrequency[0]))/(log10(fFrequency[fFrequency.size()-1]) - log10(fFrequency[0]))*fFrequency.size();
      
      if (index >= 0 && index < fFrequency.size()) {
	
	auto specVal = energyDensity[freq.size()-1-iFreq];
      
	spec[index] += specVal;
	binCount[index] += 1.;

	//cout << iFreq << " " << index << " " << freq[iFreq] << " " << fFrequency[index] << " " << specVal << " " << spec[index] << " " << binCount[index] << " " << spec[index]/binCount[index] << endl;
	
      }
      
    }
    
    spec *= 1./binCount*1./energy*1./energy;

  } else if (fFrequency.size() >= freq.size()) {

    // For binning of supplied ISRF coarser than required

    std::valarray<double> energy(fFrequency.size());
    
    energy = utl::kPlanck_SI/utl::e_SI*fFrequency;

    for (auto iFreq = 0; iFreq < fFrequency.size(); ++iFreq) {

      const auto index = (log(fFrequency[iFreq]) - log(freq[0]))/(log(freq[freq.size()-1]) - log(freq[0]))*freq.size();

      if (index >= 0 && index < freq.size()-1) {

	const auto interpolant = (log(fFrequency[iFreq]) - log(freq[index]))/(log(freq[index+1]) - log(freq[index]));

	spec[iFreq] = exp(log(energyDensity[freq.size()-1-index])*(1. - interpolant) + log(energyDensity[freq.size()-1-(index+1)])*interpolant);

	//#pragma omp critical
	//{
	//cout << "target gtr base: " <<  iFreq << " " << index << " " << fFrequency[iFreq] << " " << freq[index] << " " << freq[index+1] << " " << interpolant << " " << energyDensity[freq.size()-1-index] << " " << energyDensity[freq.size()-1-(index+1)] << " " << spec[iFreq] << endl;
	//}
	
      }
      
    }

    spec *= 1./energy*1./energy;

  } 

  return spec;

}

const std::valarray<double>
RadiationField::GetNumberDensity3D(const ThreeVector& pos, 
				   const STELLARCOMPONENT component) {

  auto cidx(-1);

  if (component == RadiationField::TOTAL)
    cidx = 0;

  if (component == RadiationField::OPTICAL)
    cidx = 1; // Explicit mapping because we only store total, optical, infrared for 3D models

  if (component == RadiationField::INFRARED)
    cidx = 2; // Explicit mapping because we only store total, optical, infrared for 3D models

  assert(cidx > -1 && cidx < 3);

  std::valarray<double> spec(0., fFrequency.size());
  
  auto r = sqrt(pos.x()*pos.x() + pos.y()*pos.y()), phi = atan2(pos.y(), pos.x()), phiDeg = (phi < 0 ? phi + 2.*utl::kPi : (phi > 2.*utl::kPi ? phi - 2.*utl::kPi : phi))*180./utl::kPi, z = pos.z();

  // Want to check for this early because we rely on the corrected angle and it should be within this range!
  assert(phiDeg >= 0.);
  assert(phiDeg <= 360.); 
  
  //std::cout << pos.x() << " " << pos.y() << " " << r << " " << phiDeg << " " << z << std::endl;

  if (r > fRData[fRData.size()-1] || z > fZData[fZData.size()-1] || z < fZData[0]) {
    
    return spec;
    
  }
  
  // Deal with the case where the rarray stops short of R = 0. We assume that 
  // the supplied radiation field from R = 0 to fRData[0] is constant and 
  // taken as the value at fRData[0], so just set r to fRData[0] ...

  if (r >= 0. && r <= fRData[0])
    r = fRData[0];

  const auto& wl = fWavelength;
  
  std::valarray<double> freq(0., wl.size());
  
  for (size_t iWl = 0; iWl < wl.size(); ++iWl) {
    freq[iWl] = utl::kSpeedOfLight_SI*1./(wl[wl.size()-1-iWl]*utl::micron/utl::m);
  }
  
  std::valarray<double> energydensity(0., wl.size());  

  // Trilinear interpolation 
   
  const auto ridx = std::upper_bound(&fRData[0], &fRData[fRData.size()-1], r) - &fRData[0];
  
  const auto dr = (fRData[ridx] - fRData[ridx-1]);
  assert(dr > 0.);

  const auto rcoeff = (fRData[ridx] - r)/dr;
  
  const auto phiidx = std::upper_bound(&fPhiData[0], &fPhiData[fPhiData.size()-1], phiDeg) - &fPhiData[0];
  
  const auto dphi = (phiDeg < fPhiData[fPhiData.size()-1]) ? (fPhiData[phiidx] - fPhiData[phiidx-1]) : (360. - fPhiData[fPhiData.size()-1]);
  assert(dphi > 0.);

  const auto phicoeff = (phiDeg < fPhiData[fPhiData.size()-1] && dphi > 0.) ? (fPhiData[phiidx] - phiDeg)/dphi : (360. - phiDeg)/dphi;
  
  const auto zidx = std::upper_bound(&fZData[0], &fZData[fZData.size()-1], z) - &fZData[0];
  
  const auto dz = (fZData[zidx] - fZData[zidx-1]);
  assert(dz > 0.);

  const auto zcoeff = (fZData[zidx] - z)/dz;
 
  if (phi > fPhiData[fPhiData.size()-1]) { // This deals with the wrap around in phi coordinate between the last bin in the grid and 360./0. deg

    const size_t idx1 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (phiidx-1);
    
    const auto& energydensity1 = *fEnergyDensity3D[idx1][cidx];

    const size_t idx2 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (0);
    
    const auto& energydensity2 = *fEnergyDensity3D[idx2][cidx];

    const size_t idx3 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (phiidx-1);
    
    const auto& energydensity3 = *fEnergyDensity3D[idx3][cidx];

    const size_t idx4 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (0);
    
    const auto& energydensity4 = *fEnergyDensity3D[idx4][cidx];
          
    const size_t idx5 = (zidx)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (phiidx-1);
    
    const auto& energydensity5 = *fEnergyDensity3D[idx5][cidx];
    
    const size_t idx6 = (zidx)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (0);
    
    const auto& energydensity6 = *fEnergyDensity3D[idx6][cidx];
    
    const size_t idx7 = (zidx)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (phiidx-1);
    
    const auto& energydensity7 = *fEnergyDensity3D[idx7][cidx];
    
    const size_t idx8 = (zidx)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (0);
    
    const auto& energydensity8 = *fEnergyDensity3D[idx8][cidx];
        
    energydensity = energydensity1*zcoeff*rcoeff*phicoeff + 
      energydensity2*zcoeff*rcoeff*(1. - phicoeff) + 
      energydensity3*zcoeff*(1. - rcoeff)*phicoeff + 
      energydensity4*zcoeff*(1. - rcoeff)*(1. - phicoeff) + 
      energydensity5*(1. - zcoeff)*rcoeff*phicoeff + 
      energydensity6*(1. - zcoeff)*rcoeff*(1. - phicoeff) +
      energydensity7*(1. - zcoeff)*(1. - rcoeff)*phicoeff + 
      energydensity8*(1. - zcoeff)*(1. - rcoeff)*(1. - phicoeff);

  } else {

    const size_t idx1 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (phiidx-1);
    
    const auto& energydensity1 = *fEnergyDensity3D[idx1][cidx];

    const size_t idx2 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (phiidx);
    
    const auto& energydensity2 = *fEnergyDensity3D[idx2][cidx];

    const size_t idx3 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (phiidx-1);
    
    const auto& energydensity3 = *fEnergyDensity3D[idx3][cidx];

    const size_t idx4 = (zidx-1)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (phiidx);
    
    const auto& energydensity4 = *fEnergyDensity3D[idx4][cidx];
          
    const size_t idx5 = (zidx)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (phiidx-1);
    
    const auto& energydensity5 = *fEnergyDensity3D[idx5][cidx];
    
    const size_t idx6 = (zidx)*fRData.size()*fPhiData.size() + (ridx-1)*fPhiData.size() + (phiidx);
    
    const auto& energydensity6 = *fEnergyDensity3D[idx6][cidx];
    
    const size_t idx7 = (zidx)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (phiidx-1);
    
    const auto& energydensity7 = *fEnergyDensity3D[idx7][cidx];
    
    const size_t idx8 = (zidx)*fRData.size()*fPhiData.size() + (ridx)*fPhiData.size() + (phiidx);
    
    const auto& energydensity8 = *fEnergyDensity3D[idx8][cidx];
        
    energydensity = energydensity1*zcoeff*rcoeff*phicoeff + 
      energydensity2*zcoeff*rcoeff*(1. - phicoeff) + 
      energydensity3*zcoeff*(1. - rcoeff)*phicoeff + 
      energydensity4*zcoeff*(1. - rcoeff)*(1. - phicoeff) + 
      energydensity5*(1. - zcoeff)*rcoeff*phicoeff + 
      energydensity6*(1. - zcoeff)*rcoeff*(1. - phicoeff) +
      energydensity7*(1. - zcoeff)*(1. - rcoeff)*phicoeff + 
      energydensity8*(1. - zcoeff)*(1. - rcoeff)*(1. - phicoeff);

  }
  
  /*const double* rVal = std::lower_bound(&fRData[0], &fRData[fRData.size()-1], r);
  
  const int rIndex = rVal - &fRData[0];
  
  auto rCoeff(0.), dR(0.);
  
  if (r <= fRData[rIndex] && r >= fRData[rIndex])
    rCoeff = 1.;
  else if (r < fRData[fRData.size()-1]) {
    
    dR = fRData[rIndex] - fRData[rIndex-1];
    
    rCoeff = (fRData[rIndex] - r)/dR;
    
  }
  
  const double* phiVal = std::lower_bound(&fPhiData[0], &fPhiData[fPhiData.size()-1], phiDeg);
  
  const int phiIndex = phiVal - &fPhiData[0];
  
  auto phiCoeff(0.), dPhi(0.);
  
  if (phiDeg <= fPhiData[phiIndex] && phiDeg >= fPhiData[phiIndex])
    phiCoeff = 1.;
  else if (phiDeg < fPhiData[fPhiData.size()-1]) {
    
    dPhi = fPhiData[phiIndex] - fPhiData[phiIndex-1];
    
    phiCoeff = (fPhiData[phiIndex] - phiDeg)/dPhi;
    
  }
  
  const double* zVal = std::lower_bound(&fZData[0], &fZData[fZData.size()-1], z);
  
  const int zIndex = zVal - &fZData[0];
  
  auto zCoeff(0.), dZ(0.);
  
  if (z <= fZData[zIndex] && z >= fZData[zIndex]) 
    zCoeff = 1.;
  else if (z < fZData[fZData.size()-1] && z > fZData[0]) {
    
    dZ = fZData[zIndex] - fZData[zIndex-1];
    
    zCoeff = (fZData[zIndex] - z)/dZ;
    
  }
  
  //auto index = zIndex*fRData.size()*fPhiData.size() + rIndex*fPhiData.size() + phiIndex;
  //assert(index >= 0 && index < fEnergyDensity3D.size());
  
  //cout << z << " " << zIndex << " " << dZ << " " << zCoeff << " " << fZData[zIndex] << endl;
  
  //exit(0);
  
  auto& wl = fWavelength;
  
  std::valarray<double> freq(0., wl.size());
  
  for (auto iWl = 0; iWl < wl.size(); ++iWl) {
    freq[iWl] = utl::kSpeedOfLight_SI*1./(wl[wl.size()-1-iWl]*utl::micron/utl::m);
    //std::cout << iWl << " " << freq[iWl] << " " << pos << " " << component << std::endl;
  }
  
  //std::cout << "Pos : " << pos << " " << component << std::endl;

  //std::cout << "Indices: " << fZData[zIndex] << " " << zIndex << " " << zCoeff << " " << z << " " << fRData[rIndex] << " " << rIndex << " " << rCoeff << " " << r << " " << fPhiData[phiIndex] << " " << phiIndex << " " << phiCoeff << " " << phiDeg << std::endl;

  std::valarray<double> energyDensity(0., wl.size());  
  
  if (rCoeff >= 1. && zCoeff >= 1. && phiCoeff >= 1.) {
    
    auto index = zIndex*fRData.size()*fPhiData.size() + rIndex*fPhiData.size() + phiIndex;

    energyDensity = *fEnergyDensity3D[index][component];

  } else if (rCoeff >= 1.) {

    if (phiCoeff >= 1.) {

      auto index1 = (zIndex-1)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
      auto& energyDensity1 = *fEnergyDensity3D[index1][component];
      
      auto index2 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
      auto& energyDensity2 = *fEnergyDensity3D[index2][component];
      
      energyDensity = energyDensity1*zCoeff + energyDensity2*(1. - zCoeff);
      
    } else if (zCoeff >= 1.) {
      
      auto index1 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex-1);
      auto& energyDensity1 = *fEnergyDensity3D[index1][component];
      
      auto index2 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
      auto& energyDensity2 = *fEnergyDensity3D[index2][component];
      
      energyDensity = energyDensity1*phiCoeff + energyDensity2*(1. - phiCoeff);
      
    } else {
      
      auto index1 = (zIndex-1)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex-1);
      auto& energyDensity1 = *fEnergyDensity3D[index1][component];
      
      auto index2 = (zIndex-1)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
      auto& energyDensity2 = *fEnergyDensity3D[index2][component];
     
      auto index3 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex-1);
      auto& energyDensity3 = *fEnergyDensity3D[index3][component];
      
      auto index4 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
      auto& energyDensity4 = *fEnergyDensity3D[index4][component];
            
      energyDensity = energyDensity1*phiCoeff*zCoeff + energyDensity2*(1. - phiCoeff)*zCoeff + energyDensity3*phiCoeff*(1. - zCoeff) + energyDensity4*(1. - phiCoeff)*(1. - zCoeff);

    }

  } else if (phiCoeff >= 1.) {

    if (rCoeff >= 1.) {

      auto index1 = (zIndex-1)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
      auto& energyDensity1 = *fEnergyDensity3D[index1][component];
      
      auto index2 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
      auto& energyDensity2 = *fEnergyDensity3D[index2][component];
      
      energyDensity = energyDensity1*zCoeff + energyDensity2*(1. - zCoeff);
      
    } else if (zCoeff >= 1.) {
      
      auto index1 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex-1)*fPhiData.size() + (phiIndex);
      auto& energyDensity1 = *fEnergyDensity3D[index1][component];
      
      auto index2 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
      auto& energyDensity2 = *fEnergyDensity3D[index2][component];
      
      energyDensity = energyDensity1*rCoeff + energyDensity2*(1. - rCoeff);
      
    } else {
      
      auto index1 = (zIndex-1)*fRData.size()*fPhiData.size() + (rIndex-1)*fPhiData.size() + (phiIndex);
      auto& energyDensity1 = *fEnergyDensity3D[index1][component];
      
      auto index2 = (zIndex-1)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
      auto& energyDensity2 = *fEnergyDensity3D[index2][component];
     
      auto index3 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex-1)*fPhiData.size() + (phiIndex);
      auto& energyDensity3 = *fEnergyDensity3D[index3][component];
      
      auto index4 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
      auto& energyDensity4 = *fEnergyDensity3D[index4][component];
            
      energyDensity = energyDensity1*rCoeff*zCoeff + energyDensity2*(1. - rCoeff)*zCoeff + energyDensity3*rCoeff*(1. - zCoeff) + energyDensity4*(1. - rCoeff)*(1. - zCoeff);

    }

  } else if (zCoeff >= 1.) {

    if (phiCoeff >= 1.) {

      auto index1 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex-1)*fPhiData.size() + (phiIndex);
      auto& energyDensity1 = *fEnergyDensity3D[index1][component];
      
      auto index2 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
      auto& energyDensity2 = *fEnergyDensity3D[index2][component];
      
      energyDensity = energyDensity1*rCoeff + energyDensity2*(1. - rCoeff);
      
    } else if (rCoeff >= 1.) {
      
      auto index1 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex-1);
      auto& energyDensity1 = *fEnergyDensity3D[index1][component];
      
      auto index2 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
      auto& energyDensity2 = *fEnergyDensity3D[index2][component];
      
      energyDensity = energyDensity1*phiCoeff + energyDensity2*(1. - phiCoeff);
      
    } else {
      
      auto index1 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex-1)*fPhiData.size() + (phiIndex-1);
      auto& energyDensity1 = *fEnergyDensity3D[index1][component];
      
      auto index2 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex-1)*fPhiData.size() + (phiIndex);
      auto& energyDensity2 = *fEnergyDensity3D[index2][component];
     
      auto index3 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex-1);
      auto& energyDensity3 = *fEnergyDensity3D[index3][component];
      
      auto index4 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
      auto& energyDensity4 = *fEnergyDensity3D[index4][component];
            
      energyDensity = energyDensity1*rCoeff*phiCoeff + energyDensity2*rCoeff*(1. - phiCoeff) + energyDensity3*(1. - rCoeff)*phiCoeff + energyDensity4*(1. - rCoeff)*(1. - phiCoeff);

    }

  } else {
  
    //std::cout << "Indices: " << fZData[zIndex] << " " << zIndex << " " << zCoeff << " " << z << " " << fRData[rIndex] << " " << rIndex << " " << rCoeff << " " << r << " " << fPhiData[phiIndex] << " " << phiIndex << " " << phiCoeff << " " << phiDeg << std::endl;
  
    auto index1 = (zIndex-1)*fRData.size()*fPhiData.size() + (rIndex-1)*fPhiData.size() + (phiIndex-1);
    //std::cout << "Index1: " << index1 << std::endl;
    auto& energyDensity1 = *fEnergyDensity3D[index1][component];

    auto index2 = (zIndex-1)*fRData.size()*fPhiData.size() + (rIndex-1)*fPhiData.size() + (phiIndex);
    //std::cout << "Index2: " << index2 << std::endl;
    auto& energyDensity2 = *fEnergyDensity3D[index2][component];

    auto index3 = (zIndex-1)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex-1);
    //std::cout << "Index3: " << index3 << std::endl;
    auto& energyDensity3 = *fEnergyDensity3D[index3][component];

    auto index4 = (zIndex-1)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
    //std::cout << "Index4: " << index4 << std::endl;
    auto& energyDensity4 = *fEnergyDensity3D[index4][component];

    auto index5 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex-1)*fPhiData.size() + (phiIndex-1);
    //std::cout << "Index5: " << index5 << std::endl;
    auto& energyDensity5 = *fEnergyDensity3D[index5][component];

    auto index6 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex-1)*fPhiData.size() + (phiIndex);
    //std::cout << "Index6: " << index6 << std::endl;
    auto& energyDensity6 = *fEnergyDensity3D[index6][component];

    auto index7 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex-1);
    //std::cout << "Index7: " << index7 << std::endl;
    auto& energyDensity7 = *fEnergyDensity3D[index7][component];

    auto index8 = (zIndex)*fRData.size()*fPhiData.size() + (rIndex)*fPhiData.size() + (phiIndex);
    //std::cout << "Index8: " << index8 << std::endl;
    auto& energyDensity8 = *fEnergyDensity3D[index8][component];

    energyDensity = energyDensity1*zCoeff*rCoeff*phiCoeff + energyDensity2*zCoeff*rCoeff*(1. - phiCoeff) + energyDensity3*zCoeff*(1. - rCoeff)*phiCoeff + energyDensity4*zCoeff*(1. - rCoeff)*(1. - phiCoeff) + energyDensity5*(1. - zCoeff)*rCoeff*phiCoeff + energyDensity6*(1. - zCoeff)*rCoeff*(1. - phiCoeff) + energyDensity7*(1. - zCoeff)*(1. - rCoeff)*phiCoeff + energyDensity8*(1. - zCoeff)*(1. - rCoeff)*(1. - phiCoeff);

  }
  */
  if (fFrequency.size() < freq.size()) { 

    // For binning of supplied ISRF finer than required
    
    std::valarray<double> binCount(0., fFrequency.size()), energy(fFrequency.size());
    
    energy = utl::kPlanck_SI/utl::e_SI*fFrequency;
    
    for (auto iFreq = 0; iFreq < freq.size(); ++iFreq) {
      
      const auto index = (log10(freq[iFreq]) - log10(fFrequency[0]))/(log10(fFrequency[fFrequency.size()-1]) - log10(fFrequency[0]))*fFrequency.size();
      
      if (index >= 0 && index < fFrequency.size()) {
	
	const auto specVal = energydensity[freq.size()-1-iFreq];
      
	spec[index] += specVal;
	binCount[index] += 1.;

	//cout << iFreq << " " << index << " " << freq[iFreq] << " " << fFrequency[index] << " " << specVal << " " << spec[index] << " " << binCount[index] << " " << spec[index]/binCount[index] << endl;
	
      }
      
    }
    
    spec *= 1./binCount*1./energy*1./energy;

  } else if (fFrequency.size() >= freq.size()) {

    // For binning of supplied ISRF coarser than required

    std::valarray<double> energy(fFrequency.size());
    
    energy = utl::kPlanck_SI/utl::e_SI*fFrequency;

    for (auto iFreq = 0; iFreq < fFrequency.size(); ++iFreq) {

      const auto index = (log(fFrequency[iFreq]) - log(freq[0]))/(log(freq[freq.size()-1]) - log(freq[0]))*freq.size();

      if (index >= 0 && index < freq.size()-1) {

	const auto interpolant = (log(fFrequency[iFreq]) - log(freq[index]))/(log(freq[index+1]) - log(freq[index]));

	spec[iFreq] = exp(log(energydensity[freq.size()-1-index])*(1. - interpolant) + log(energydensity[freq.size()-1-(index+1)])*interpolant);

	//#pragma omp critical
	//{
	//cout << "target gtr base: " <<  iFreq << " " << index << " " << fFrequency[iFreq] << " " << freq[index] << " " << freq[index+1] << " " << interpolant << " " << energyDensity[freq.size()-1-index] << " " << energyDensity[freq.size()-1-(index+1)] << " " << spec[iFreq] << endl;
	//}
	
      }
      
    }

    spec *= 1./energy*1./energy;

  } 

  return spec;

}

void RadiationField::ReadRadiationField2D(const std::string& filename) {

  //#pragma omp critical 
  {

    //if (filename.empty()) {

    //INFO("Empty radiation field filename passed to ReadRadiationField.");
    //exit(-1);

    //}
  
    //ostringstream buf;
    //buf << "Reading from " << filename;
    //INFO(buf.str());
    
    ClearData();
    
    const std::string prefix = filename.substr(0, filename.find_last_of("/")+1);
    
    //cout << prefix << endl;
    
    //fPrefix = prefix;
    
    std::ifstream rfFile(filename.c_str());

    if (!rfFile) {
      std::ostringstream buf;
      buf.str("");
      buf << "Error reading from file " << filename;
      ERROR(buf.str());
      throw(std::invalid_argument(buf.str()));
    }

    unsigned long volumeElements, wlBins, numFilters;
    double luminosity, dustMass;
    
    rfFile >> volumeElements >> wlBins >> numFilters;
    rfFile >> luminosity >> dustMass;
        
    fLuminosity = luminosity;
    fDustMass = dustMass;
    
    fWavelength.resize(wlBins);
    fWavelength = 0;
    
    unsigned long stellarComponents;
    
    rfFile >> stellarComponents;
    
    fStellarComponentLuminosity.resize(stellarComponents);
    fStellarComponentName.resize(stellarComponents);
    
    for (unsigned long i = 0; i < stellarComponents; ++i) {
      
      double componentLuminosity;
      std::string componentName;
      
      rfFile >> componentName >> componentLuminosity;
      
      fStellarComponentLuminosity[i] = componentLuminosity;
      fStellarComponentName[i] = componentName;
      
    }
    
    unsigned long geometry;
    double modelRegionData[6];
    
    rfFile >> geometry;
    rfFile >> modelRegionData[0] >> modelRegionData[1] >> modelRegionData[2] >> modelRegionData[3] >> modelRegionData[4] >> modelRegionData[5];
    
    std::vector<ThreeVector> posVec, rangeVec;
    
    // This is assuming that the individual elements are ordered correctly 
    // and that there are no duplicates

    for (unsigned long i = 0; i < volumeElements; ++i) {
      
      double index;
      double x, y, z, dX, dY, dZ;
      
      rfFile >> index;
      rfFile >> x >> y >> z;
      rfFile >> dX >> dY >> dZ;
      
      ThreeVector pos(x, y, z), range(dX, dY, dZ);
      
      //cout << pos << " " << range << endl;
      
      if (!posVec.size()) {
	
	posVec.push_back(pos);
	rangeVec.push_back(range);
	
      } else if (fabs(posVec[0].x() - pos.x()) <= dX) {
	
	posVec.push_back(pos);
	rangeVec.push_back(range);
	
      }
      
      if (fabs(posVec[0].x() - pos.x()) > dX || i == volumeElements - 1) {
	
	fPositionData2D.push_back(posVec);
	fRangeData2D.push_back(rangeVec);
	
	posVec.clear();
	posVec.push_back(pos);
	
	rangeVec.clear();
	rangeVec.push_back(range);
	
      }
      
      std::string filterFilename, filterCountFilename, directFilename, scatteredFilename, transientFilename, thermalFilename, totalFilename, filterFluxFilename, fluxFilename;
      
      rfFile >> filterFilename;
      rfFile >> filterCountFilename;
      rfFile >> directFilename;
      rfFile >> scatteredFilename;
      rfFile >> transientFilename;
      rfFile >> thermalFilename;
      rfFile >> totalFilename;
      rfFile >> filterFluxFilename;
      rfFile >> fluxFilename;

      //cout << filterFilename << endl;
      
      fNumberOfComponents = 7;
      
      //fCacheBuilt.resize(fNumberOfComponents);
      //fCacheBuilt = false;
      
      fFilenameData[pos].resize(fNumberOfComponents+1);
      
      fFilenameData[pos][0] = totalFilename;
      fFilenameData[pos][1] = directFilename;
      fFilenameData[pos][2] = scatteredFilename;
      fFilenameData[pos][3] = transientFilename;
      fFilenameData[pos][4] = thermalFilename;
      fFilenameData[pos][5] = filterFilename;
      fFilenameData[pos][6] = filterCountFilename;
      fFilenameData[pos][7] = fluxFilename;
      
    }
    
    for (unsigned int i = 0; i < fPositionData2D.size(); ++i) {
      
      std::vector< std::vector<std::string>* > pStrData;
      
      std::vector< std::vector<Skymap<double>* > > pSkymapData;
      
      std::vector< std::vector< std::valarray<double>* > > pEnergyDensityData;
      
      for (unsigned int j = 0; j < fPositionData2D[i].size(); ++j) {
	
	const ThreeVector& pos = fPositionData2D[i][j];
	
	//std::map<ThreeVector, std::vector<std::string> >::iterator fIt = fFilenameData.find(pos);
	
	auto fIt = fFilenameData.find(pos);

	pStrData.push_back(&fIt->second);
	
	std::vector< Skymap<double>* > skymaps;
	skymaps.resize(pStrData.back()->size()); 
	skymaps.assign(pStrData.back()->size(), 0);
	
	pSkymapData.push_back(skymaps);
	
	// Read in the flux file at this stage -- we're ready to go at exit from
	// this method if this is all that the user will be requesting.
	
	std::vector< std::valarray<double>* > energyDensity;
	
	for (int iC = 0; iC < 5; ++iC) {
	  
	  std::valarray<double>* vec = new std::valarray<double>(0., wlBins);
	  energyDensity.push_back(vec);
	  
	}
	
	const std::string fname = prefix + "/" + pStrData.back()->back();
	
	CCfits::FITS fluxFile(fname, CCfits::Read);
	
	CCfits::ExtHDU& table = fluxFile.currentExtension(); // It should be the only one in this file ...
	
	std::vector<double> wl, total, direct, scattered, transient, thermal;
	
	table.column(1).read(wl, 1, wlBins);
	table.column(2).read(total, 1, wlBins);
	table.column(3).read(direct, 1, wlBins);
	table.column(4).read(scattered, 1, wlBins);
	table.column(5).read(transient, 1, wlBins);
	table.column(6).read(thermal, 1, wlBins);
	
	for (int iWl = 0; iWl < wlBins; ++iWl) {
	  
	  fWavelength[iWl] = wl[iWl];
	  
	  (*(energyDensity[0]))[iWl] = total[iWl];
	  (*(energyDensity[1]))[iWl] = direct[iWl];
	  (*(energyDensity[2]))[iWl] = scattered[iWl];
	  (*(energyDensity[3]))[iWl] = transient[iWl];
	  (*(energyDensity[4]))[iWl] = thermal[iWl];
	  
	  //cout << i << " " << j << " " << iWl << " " << wl[iWl] << " " << (*(energyDensity[0]))[iWl] << " " << (*(energyDensity[4]))[iWl] << endl;
	  
	}
	
	//exit(0);
	
	pEnergyDensityData.push_back(energyDensity);
	
      }
      
      fFilenameOrderedData2D.push_back(pStrData);
      fSkymapOrderedData2D.push_back(pSkymapData);
      fEnergyDensity2D.push_back(pEnergyDensityData);
      
    }
    
    //exit(0);
    
    fRData.resize(fPositionData2D.size());
    
    for (unsigned int i = 0; i < fPositionData2D.size(); ++i) {
      
      fRData[i] = fPositionData2D[i][0].x();

      //cout << i << " " << fRData[i] << " " << fPositionData[i].size() << endl;
            
    }
    
    fZData.resize(fPositionData2D[0].size());
    
    for (unsigned int i = 0; i < fPositionData2D[0].size(); ++i) {
      fZData[i] = fPositionData2D[0][i].z();
      //cout << i << " " << fZData[i] << endl;
    }

    //exit(0);
    
    fRRangeData.resize(fRangeData2D.size());
        
    for (unsigned int i = 0; i < fRangeData2D.size(); ++i) {
      
      fRRangeData[i] = fRangeData2D[i][0].x();
            
    }
    
    fZRangeData.resize(fRangeData2D[0].size());
    
    for (unsigned int i = 0; i < fRangeData2D[0].size(); ++i)
      fZRangeData[i] = fRangeData2D[0][i].z();
    
  }
    
}

void RadiationField::ReadRadiationField3D(const std::string& filename) {

  {
    
    ClearData();
    
    const std::string prefix = filename.substr(0, filename.find_last_of("/")+1);

    //cout << prefix << endl;

    std::ifstream rfFile(filename.c_str());

    if (!rfFile) {
      std::ostringstream buf;
      buf.str("");
      buf << "Error reading from file " << filename;
      ERROR(buf.str());
      throw(std::invalid_argument(buf.str()));
    }
    
    std::string line;
    std::getline(rfFile, line); // eat the !!3D for the first line
    
    std::getline(rfFile, line);
    std::vector<uint64_t> totals;
    { std::stringstream ss(line); for (;;) { uint64_t val; ss >> val; if (ss.fail()) break; totals.push_back(val); }} 
    auto volumeElements(totals[0]), wlBins(totals[1]), numFilters(totals[2]), numRVals(totals[3]), numPhiVals(totals[4]), numZVals(totals[5]);
    
    line.clear();
    std::getline(rfFile, line);
    std::vector<double> rData;
    { std::stringstream ss(line); for (;;) { double val; ss >> val; if (ss.fail()) break; rData.push_back(val); }}
    fRData.resize(rData.size());
    std::copy(rData.begin(), rData.end(), &fRData[0]);

    line.clear();
    std::getline(rfFile, line);
    std::vector<double> phiData;
    { std::stringstream ss(line); for (;;) { double val; ss >> val; if (ss.fail()) break; phiData.push_back(val); }}
    fPhiData.resize(phiData.size());
    std::copy(phiData.begin(), phiData.end(), &fPhiData[0]);
    
    line.clear();
    std::getline(rfFile, line);
    std::vector<double> zData;
    { std::stringstream ss(line); for (;;) { double val; ss >> val; if (ss.fail()) break; zData.push_back(val); }}
    fZData.resize(zData.size());
    std::copy(zData.begin(), zData.end(), &fZData[0]);

    double luminosity, dustMass;
    rfFile >> luminosity >> dustMass;
        
    fLuminosity = luminosity;
    fDustMass = dustMass;
    
    fWavelength.resize(wlBins);
    fWavelength = 0;
    
    uint64_t stellarComponents;
    
    rfFile >> stellarComponents;
    
    fStellarComponentLuminosity.resize(stellarComponents);
    fStellarComponentName.resize(stellarComponents);
    
    for (auto i = 0; i < stellarComponents; ++i) {
      
      double componentLuminosity;
      std::string componentName;
      
      rfFile >> componentName >> componentLuminosity;
      
      fStellarComponentLuminosity[i] = componentLuminosity;
      fStellarComponentName[i] = componentName;
      
    }
    
    uint64_t geometry;
    double modelRegionData[6];
    
    rfFile >> geometry;
    rfFile >> modelRegionData[0] >> modelRegionData[1] >> modelRegionData[2] >> modelRegionData[3] >> modelRegionData[4] >> modelRegionData[5];
    
    // Elements are ordered as z, r, phi so the indices are easily computed.

    for (auto i = 0; i < volumeElements; ++i) {
      
      double index;
      double x, y, z, dX, dY, dZ;
      
      rfFile >> index;
      rfFile >> x >> y >> z;
      rfFile >> dX >> dY >> dZ;
      
      ThreeVector pos(x, y, z), range(dX, dY, dZ);

      std::string filterFilename, filterCountFilename, directFilename, scatteredFilename, transientFilename, thermalFilename, totalFilename, opticalfilename, infraredfilename, filterFluxFilename, fluxFilename;
      
      rfFile >> filterFilename;
      rfFile >> filterCountFilename;
      rfFile >> directFilename;
      rfFile >> scatteredFilename;
      rfFile >> transientFilename;
      rfFile >> thermalFilename;
      rfFile >> totalFilename;
      rfFile >> opticalfilename;
      rfFile >> infraredfilename;
      rfFile >> filterFluxFilename;
      rfFile >> fluxFilename;

      //std::cout << i << " " << totalFilename << std::endl;
      
      fNumberOfComponents = 3;
      
      //fCacheBuilt.resize(fNumberOfComponents);
      //fCacheBuilt = false;
      
      fFilenameData[pos].resize(fNumberOfComponents+1);
      
      fFilenameData[pos][0] = totalFilename;
      fFilenameData[pos][1] = opticalfilename;//directFilename;
      fFilenameData[pos][2] = infraredfilename;//scatteredFilename;
      fFilenameData[pos][3] = fluxFilename;//transientFilename;
      //fFilenameData[pos][4] = thermalFilename;
      //fFilenameData[pos][5] = filterFilename;
      //fFilenameData[pos][6] = filterCountFilename;
      //fFilenameData[pos][7] = fluxFilename;

      const auto fIt = fFilenameData.find(pos);

      fFilenameData3D.push_back(&fIt->second);
            
      //cout << i << " " << pos << " " << posVecArr.size() << " " << posVecArr.back().size() << " " << posVecArr.back()[ii] << endl;

      std::vector< std::valarray<double>* > energyDensity;

      // Read in the flux file at this stage -- we're ready to go at exit from
      // this method if this is all that the user will be requesting.
	  	      
      for (auto iC = 0; iC < 3; ++iC) {
		
	std::valarray<double>* vec = new std::valarray<double>(0., wlBins);
	energyDensity.push_back(vec);
			
      }
	      
      const std::string fname = prefix + "/" + fluxFilename;
	  	      
      //std::cout << fname << std::endl;

      CCfits::FITS fluxFile(fname, CCfits::Read);
	      
      CCfits::ExtHDU& table = fluxFile.currentExtension(); // It should be the only one in this file ...
      
      std::vector<double> wl, total, direct, scattered, transient, thermal;
	      
      table.column(1).read(wl, 1, wlBins);
      table.column(2).read(total, 1, wlBins);
      table.column(3).read(direct, 1, wlBins);
      table.column(4).read(scattered, 1, wlBins);
      table.column(5).read(transient, 1, wlBins);
      table.column(6).read(thermal, 1, wlBins);
      
      for (auto iWl = 0; iWl < wlBins; ++iWl) {
		
	fWavelength[iWl] = wl[iWl];
	
	(*(energyDensity[0]))[iWl] = total[iWl];
	(*(energyDensity[1]))[iWl] = direct[iWl] + scattered[iWl];
	(*(energyDensity[2]))[iWl] = transient[iWl] + thermal[iWl];//scattered[iWl];
	//(*(energyDensity[3]))[iWl] = transient[iWl];
	//(*(energyDensity[4]))[iWl] = thermal[iWl];
	
	//std::cout << i << " " << pos << " " << iWl << " " << wl[iWl] << " " << (*(energyDensity[0]))[iWl] << " " << (*(energyDensity[4]))[iWl] << std::endl;
	
      }

      /*if (fabs(z) >= 0. && fabs(z) <= 0.) {
      auto sum(0.);
      for (auto i = 0; i < energyDensity[0]->size(); ++i)
	sum += (*(energyDensity[0]))[i];
      std::cout << x << " " << y << " " << sum*log(wl[1]/wl[0]) << std::endl;
      }
      */

      fEnergyDensity3D.push_back(energyDensity);

      std::vector< Skymap<double>* > skymaps(energyDensity.size(),  nullptr);
      fSkymapData3D.push_back(skymaps); // Built in BuildSkymapCache3D() as needed

    }

  }

}

void RadiationField::BuildSkymapCache2D() { //const unsigned int component) {

#pragma omp critical (RadiationFieldCache2D)
  {

    INFO("Entry");

    if (!fCacheBuilt) {
        
      std::ostringstream buf;
      buf << "Building skymap cache";// for component" << component;
      INFO(buf.str());
      
      //assert(component < fNumberOfComponents);
      
      std::valarray<double> energy(fFrequency.size());
      
      energy = utl::kPlanck_SI/utl::e_SI*fFrequency;
      
      for (auto i = 0; i < fSkymapOrderedData2D.size(); ++i) {
	
	for (auto j = 0; j < fSkymapOrderedData2D[i].size(); ++j) {
	  
	  for (size_t k = RadiationField::TOTAL; k <= RadiationField::THERMAL; ++k) {
	    
	    const std::string& filename = (*fFilenameOrderedData2D[i][j])[k];
	    
	    //cout << filename << endl;
	    
	    Skymap<double> skymap(fPrefix + filename);
	    
	    const std::valarray<double>& wl = skymap.getSpectra();
	    
	    std::valarray<double> freq(0., wl.size()), en(0., wl.size());
	    
	    for (auto iWl = 0; iWl < wl.size(); ++iWl)
	      freq[iWl] = utl::kSpeedOfLight_SI*1./(wl[wl.size()-1-iWl]*utl::micron/utl::m);
	    
	    en = utl::kPlanck_SI/utl::e_SI*freq;
	    
	    Skymap<double> skymapRebinned(fRebinnedSkymapOrder, freq);
	    
	    skymapRebinned = skymap.rebin(fRebinnedSkymapOrder);
	    
	    fSkymapOrderedData2D[i][j][k] = new Skymap<double>(fRebinnedSkymapOrder, fFrequency);	
	    
	    for (auto nPix = 0; nPix < skymapRebinned.Npix(); ++nPix) {

	      std::valarray<double> spec(0., fFrequency.size());

	      if (fFrequency.size() < freq.size()) { 

		// For binning of supplied ISRF finer than required
    
		std::valarray<double> binCount(0., fFrequency.size()), energy(fFrequency.size());
    
		energy = utl::kPlanck_SI/utl::e_SI*fFrequency;
    
		for (auto iFreq = 0; iFreq < freq.size(); ++iFreq) {
      
		  auto index = (log10(freq[iFreq]) - log10(fFrequency[0]))/(log10(fFrequency[fFrequency.size()-1]) - log10(fFrequency[0]))*fFrequency.size();
      
		  if (index >= 0 && index < fFrequency.size()) {
	
		    auto specVal = skymapRebinned[nPix][freq.size()-1-iFreq];
      
		    spec[index] += specVal;
		    binCount[index] += 1.;
      
		  }
      
		}
    
		spec *= 1./binCount*1./energy*1./energy*1./(utl::kSpeedOfLight_SI*utl::m/utl::cm);;

	      } else if (fFrequency.size() > freq.size()) {

		// For binning of supplied ISRF coarser than required
		
		std::valarray<double> energy(fFrequency.size());
    
		energy = utl::kPlanck_SI/utl::e_SI*fFrequency;
		
		for (auto iFreq = 0; iFreq < fFrequency.size(); ++iFreq) {
		  
		  auto index = (log10(fFrequency[iFreq]) - log10(freq[0]))/(log10(freq[freq.size()-1]) - log10(freq[0]))*freq.size();
		  
		  if (index >= 0 && index < freq.size()-1) {
		    
		    auto interpolant = (log(fFrequency[iFreq]) - log(freq[index]))/(log(freq[index+1]) - log(freq[index]));
		    
		    spec[iFreq] = exp(log(skymapRebinned[nPix][freq.size()-1-index])*(1. - interpolant) + log(skymapRebinned[nPix][freq.size()-1-(index+1)])*interpolant);
		    
		  }

		  spec *= 1./energy*1./energy*1./(utl::kSpeedOfLight_SI*utl::m/utl::cm);

		}
		
	      } else {
		
		// For supplied ISRF binning the same as required
		
		for (auto iFreq = 0; iFreq < freq.size(); ++iFreq)
		  spec[iFreq] = skymapRebinned[nPix][freq.size()-1-iFreq];
		
		spec *= 1./(utl::kPlanck_SI/utl::e_SI*fFrequency)*1./(utl::kPlanck_SI/utl::e_SI*fFrequency)*1./(utl::kSpeedOfLight_SI*utl::m/utl::cm);

	      }
	      
	      /*std::valarray<double> spec(0., fFrequency.size()), binCount(0., fFrequency.size());
	      
	      for (unsigned int iFreq = 0; iFreq < freq.size(); ++iFreq) {
		
		const int index = (log10(freq[iFreq]) - log10(fFrequency[0]))/(log10(fFrequency[fFrequency.size()-1]) - log10(fFrequency[0]))*fFrequency.size();
		
		if (index >= 0 && index < fFrequency.size()) {
		  
		  const double specVal = skymapRebinned[nPix][freq.size()-1-iFreq];//energyRaw[iFreq]/energyRaw[iFreq];
		  
		  spec[index] += specVal;
		  binCount[index] += 1.;
		  
		}
		
	      }
	      
	      spec *= 1./binCount*1./energy*1./energy*1./(kSpeedOfLight_SI*m/cm);
	      */

	      (*fSkymapOrderedData2D[i][j][k])[nPix] = spec; // eV^-1 cm^-3 -- per pixel
	      
	    }
	    
	    /*for (size_t iFreq = 0; iFreq < freq.size(); ++iFreq) {
	      
	      size_t index = (log10(freq[iFreq]) - log10(fFrequency[0]))/(log10(fFrequency[fFrequency.size()-1]) - log10(fFrequency[0]))*fFrequency.size();
	      
	      cout << iFreq << " " 
	      << freq[iFreq] << " " 
	      << en[iFreq] << " " 
	      << index << " " 
	      << (index >= 0 && index < fFrequency.size() ? fFrequency[index] : 0) << " "
	      << (index >= 0 && index < fFrequency.size() ? energy[index] : 0) << " " 
	      << skymap.sum(freq.size()-1-iFreq)*1./(kSpeedOfLight_SI*m/cm) << " " 
	      << skymapRebinned.sum(freq.size()-1-iFreq)*1./(kSpeedOfLight_SI*m/cm) << " "
	      << (index >= 0 && index < fFrequency.size() ? fSkymapOrderedData[i][j][component]->sum(index)*energy[index]*energy[index] : 0) << " "
	      << fFrequency.size() << endl;
	      
	      }
	      
	      exit(0);
	    */
	  }
	  
	}     
	
      }
      
      fCacheBuilt = true;//[component] = true;

    } else {

      INFO("Cache already built.");

    }

    INFO("Exit");

  }

}

void RadiationField::BuildSkymapCache3D() { 

#pragma omp critical (RadiationFieldCache3D)
  {

    INFO("Entry");
    
    if (!fCacheBuilt) {
      
      //std::ostringstream buf;
      //buf << "Building skymap cache";
      //INFO(buf.str());
      
      std::valarray<double> energy(fFrequency.size());
      
      energy = utl::kPlanck_SI/utl::e_SI*fFrequency;
      
      utl::StatusIndicator status("Build 3D skymap cache", fSkymapData3D.size());
      
      //#pragma omp parallel for schedule (dynamic)
      for (size_t i = 0; i < fSkymapData3D.size(); ++i) {
	
	for (size_t j(0); j < 3; ++j) {//uint32_t j = RadiationField::TOTAL; j <= RadiationField::THERMAL; ++j) {
	  
	  Skymap<double> skymap;
	  
	  //#pragma omp critical
	  //{	    
	  
	  // Protect FITS io with critical section because segfaults occur otherwise
	  const auto& filename = (*fFilenameData3D[i])[j];
	  
	  //std::cout << fPrefix << " " << filename << " " << fPrefix+filename << std::endl;
	  
	  skymap.load(fPrefix + filename);
	  
	  //}
	  
	  const auto& wl = skymap.getSpectra();
	  
	  std::valarray<double> freq(0., wl.size()), en(0., wl.size());
	    
	  for (auto iWl = 0; iWl < wl.size(); ++iWl)
	    freq[iWl] = utl::kSpeedOfLight_SI*1./(wl[wl.size()-1-iWl]*utl::micron/utl::m);
	    
	  en = utl::kPlanck_SI/utl::e_SI*freq;
	    
	  Skymap<double> skymapRebinned(fRebinnedSkymapOrder, freq);
	    
	  skymapRebinned = skymap.rebin(fRebinnedSkymapOrder);

	  //#pragma omp critical
	  //{
	  fSkymapData3D[i][j] = new Skymap<double>(skymapRebinned.Order(), fFrequency);//fRebinnedSkymapOrder, fFrequency);	
	      //}
	  
	  std::valarray<double> energy(fFrequency.size()), conv(fFrequency.size());
	  energy = utl::kPlanck_SI/utl::e_SI*fFrequency;
	  conv = 1./energy*1./energy*1./(utl::kSpeedOfLight_SI*utl::m/utl::cm);

#pragma omp parallel for schedule(dynamic)
	  for (auto nPix = 0; nPix < skymapRebinned.Npix(); ++nPix) {

	    std::valarray<double> spec(0., fFrequency.size());

	    if (fFrequency.size() < freq.size()) { 

	      // For binning of supplied ISRF finer than required
	      
	      std::valarray<double> binCount(0., fFrequency.size());//, energy(fFrequency.size());
	      
	      //energy = utl::kPlanck_SI/utl::e_SI*fFrequency;
	      
	      for (size_t iFreq = 0; iFreq < freq.size(); ++iFreq) {
		
		const auto index = (log10(freq[iFreq]) - log10(fFrequency[0]))/(log10(fFrequency[fFrequency.size()-1]) - log10(fFrequency[0]))*fFrequency.size();
		
		if (index >= 0 && index < fFrequency.size()) {
		  
		  auto specVal = skymapRebinned[nPix][freq.size()-1-iFreq];
		  
		  spec[index] += specVal;
		  binCount[index] += 1.;
		  
		}
		
	      }
	      
	      spec *= 1./binCount*conv;//1./energy*1./energy*1./(utl::kSpeedOfLight_SI*utl::m/utl::cm);;
	      
	    } else if (fFrequency.size() > freq.size()) {
	      
	      // For binning of supplied ISRF coarser than required
	      
	      //std::valarray<double> energy(fFrequency.size());
	      
	      //energy = utl::kPlanck_SI/utl::e_SI*fFrequency;
	      
	      for (size_t iFreq = 0; iFreq < fFrequency.size(); ++iFreq) {
		
		const auto index = (log10(fFrequency[iFreq]) - log10(freq[0]))/(log10(freq[freq.size()-1]) - log10(freq[0]))*freq.size();
		
		if (index >= 0 && index < freq.size()-1) {
		  
		  auto interpolant = (log(fFrequency[iFreq]) - log(freq[index]))/(log(freq[index+1]) - log(freq[index]));
		  
		  spec[iFreq] = exp(log(skymapRebinned[nPix][freq.size()-1-index])*(1. - interpolant) + log(skymapRebinned[nPix][freq.size()-1-(index+1)])*interpolant);
		  
		}
		
		spec *= conv;//1./energy*1./energy*1./(utl::kSpeedOfLight_SI*utl::m/utl::cm);
		
	      }
	      
	    } else {
	      
	      // For supplied ISRF binning the same as required
	      
	      for (size_t iFreq = 0; iFreq < freq.size(); ++iFreq)
		spec[iFreq] = skymapRebinned[nPix][freq.size()-1-iFreq];
	      
	      spec *= conv;//1./energy*1./energy*1./(utl::kSpeedOfLight_SI*utl::m/utl::cm);
	      
	      //1./(utl::kPlanck_SI/utl::e_SI*fFrequency)*1./(utl::kPlanck_SI/utl::e_SI*fFrequency)*1./(utl::kSpeedOfLight_SI*utl::m/utl::cm);
	      
	    }
	    
	    (*fSkymapData3D[i][j])[nPix] = spec; // eV^-1 cm^-3 -- per pixel
	      
	  }

	}

	status.refresh();
	  
      }     
      
      fCacheBuilt = true;//[component] = true;
      
    } else {
      
      INFO("Cache already built.");
      
    }
    
    INFO("Exit");
    
  }
  
}
 
void RadiationField::FlushSkymapCache2D() {//const unsigned int component) {

  INFO("Entry");

  //assert(component < fNumberOfComponents);
  
  for (auto i = 0; i < fSkymapOrderedData2D.size(); ++i) {
    
    for (auto j = 0; j < fSkymapOrderedData2D[i].size(); ++j) {
      
      for (auto k = 0; k < fNumberOfComponents; ++k) {

	delete fSkymapOrderedData2D[i][j][k];
	fSkymapOrderedData2D[i][j][k] = 0;
      
      }

    }
    
  }

  fCacheBuilt = false;//[component] = false;
  
  INFO("Exit");

}

void RadiationField::FlushSkymapCache3D() {//const unsigned int component) {

  INFO("Entry");
  
  for (auto i = 0; i < fSkymapData3D.size(); ++i) {
    
    for (auto j = 0; j < fSkymapData3D[i].size(); ++j) {

      delete fSkymapData3D[i][j];
      fSkymapData3D[i][j] = nullptr;
	  	
    }

  }

  fCacheBuilt = false;//[component] = false;

  INFO("Exit");

}

void RadiationField::ClearData() {

  //#pragma omp critical
  {

    //fEnergyDensity is created without fCacheBuilt being set
    for (auto i = 0; i < fEnergyDensity2D.size(); ++i) {
      
      for (auto j = 0; j < fEnergyDensity2D[i].size(); ++j) {
	
        for (auto k = 0; k < fEnergyDensity2D[i][j].size(); ++k) 
	  delete fEnergyDensity2D[i][j][k];
	
        fEnergyDensity2D[i][j].clear();
	
      }
      
      fEnergyDensity2D[i].clear();
      
    }
    
    fEnergyDensity2D.clear();
    
    for (auto i = 0; i < fPositionData2D.size(); ++i) {
     
      fPositionData2D[i].clear();
      fRangeData2D[i].clear();
      
    }
    
    for (auto i = 0; i < fFilenameOrderedData2D.size(); ++i)
      fFilenameOrderedData2D[i].clear(); // Don't own the pointers
   
    fFilenameOrderedData2D.clear();
   
    for (auto i = 0; i < fSkymapOrderedData2D.size(); ++i) {
      
      for (auto j = 0; j < fSkymapOrderedData2D[i].size(); ++j) {
	
	if (fCacheBuilt) {
	  for (auto k = 0; k < fSkymapOrderedData2D[i][j].size(); ++k) 
	    delete fSkymapOrderedData2D[i][j][k];
	}
	
	fSkymapOrderedData2D[i][j].clear();
	
      }
      
      fSkymapOrderedData2D[i].clear();
      
    }
    
    fSkymapOrderedData2D.clear();
    
    for (auto i = 0; i < fEnergyDensity3D.size(); ++i) {

      for (auto j = 0; j < fEnergyDensity3D[i].size(); ++j)
	    delete fEnergyDensity3D[i][j];
	
      fEnergyDensity3D[i].clear();
	
    }
   
    if (fCacheBuilt)
      for (auto i = 0; i < fSkymapData3D.size(); ++i) {
      
	for (auto j = 0; j < fSkymapData3D[i].size(); ++i) {
	  delete fSkymapData3D[i][j];
	  fSkymapData3D[i][j] = nullptr;
	}

	fSkymapData3D[i].clear();

      }
    
    fCacheBuilt = false;
    
  }
  
}
