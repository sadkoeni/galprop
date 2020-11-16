#ifndef _rf_RadiationField_h_
#define _rf_RadiationField_h_

#include <map>
#include <string>
#include <valarray>
#include <vector>

#include <CLHEP/Vector/ThreeVector.h>

#include <Skymap.h>

namespace rf {

  // Explicit assumption made is that the supplied radiation field
  // is generated on a regular spatial grid. That is, there are n*m entries 
  // in the grid. For each entry in the grid, there is a range of validity 
  // (given by the cell size). Outside this range the radiation field is 
  // determined by interpolation within the grid. For points outside the grid
  // the radiation field is considered zero

  class RadiationField {

  public:

    enum STELLARCOMPONENT { TOTAL, DIRECT, SCATTERED, TRANSIENT, THERMAL, OPTICAL, INFRARED };

#ifdef CLHEP_V1_8
    typedef Hep3Vector ThreeVector;
#else
    typedef CLHEP::Hep3Vector ThreeVector;
#endif

    RadiationField();  
    
    RadiationField(const std::string& filename,
		   const std::valarray<double>& freq,
		   int rebinnedSkymapOrder = 0);

    ~RadiationField();

    // Note: The method below is *very* expensive if the desired 
    // healpix order is different from the rebinned order supplied
    // in the constructor above. Use with caution!

    const Skymap<double> GetSkymap(const ThreeVector& pos, 
				   const STELLARCOMPONENT component,
				   const int healpixOrder);

    // This method is *much* faster if you just want the number density

    const std::valarray<double> GetNumberDensity(const ThreeVector& pos,
						 const STELLARCOMPONENT component); 

    // Use these to get information about the boundaries in the binning 
    // for the supplied ISRF data files

    const std::valarray<double>& GetBoundaryR() const { return fRData; }
    const std::valarray<double>& GetBoundaryZ() const { return fZData; }
    const std::valarray<double>& GetBoundaryX() const { return fXData; }
    const std::valarray<double>& GetBoundaryY() const { return fYData; }

    bool Is3D() { return (fUse3D == true); }
    
  private:

    const Skymap<double> GetSkymap(const ThreeVector& pos,
				   const STELLARCOMPONENT component);

    const Skymap<double> GetSkymap2D(const ThreeVector& pos,
				     const STELLARCOMPONENT component);
    const Skymap<double> GetSkymap3D(const ThreeVector& pos,
				     const STELLARCOMPONENT component);

    const std::valarray<double> GetNumberDensity2D(const ThreeVector& pos,
						   const STELLARCOMPONENT component);
    const std::valarray<double> GetNumberDensity3D(const ThreeVector& pos,
						   const STELLARCOMPONENT component);

    bool fCacheBuilt = false, fUse2D = false, fUse3D = false;

    std::string fPrefix;

    int fRebinnedSkymapOrder = -1;

    unsigned int fNumberOfComponents = 0;

    double fLuminosity = 0., fDustMass = 0.;

    std::valarray<double> fFrequency, fWavelength, fStellarComponentLuminosity;

    std::vector<std::string> fStellarComponentName;

    std::map<ThreeVector, std::vector<std::string> > fFilenameData;
  
    std::vector< std::vector<ThreeVector> > fPositionData2D, fRangeData2D;

    std::vector< std::vector< std::vector<ThreeVector> > > fPositionData3D, fRangeData3D;

    std::valarray<double> fRData, fPhiData, fXData, fYData, fZData, fRRangeData, fXRangeData, fYRangeData, fZRangeData;

    std::vector< std::vector< std::vector<std::string>* > > fFilenameOrderedData2D;

    std::vector< std::vector< std::vector<Skymap<double>* > > > fSkymapOrderedData2D;

    std::vector< std::vector< std::vector< valarray<double>* > > > fEnergyDensity2D;

    std::vector< std::vector<std::string>* > fFilenameData3D;
    
    std::vector< std::vector< Skymap<double>* > > fSkymapData3D;
    std::vector< std::vector< std::valarray<double>* > > fEnergyDensity3D;

    //std::vector< std::vector< std::vector< std::vector<std::string>* > > > fFilenameOrderedData3D;

    //std::vector< std::vector< std::vector< std::vector<Skymap<double>* > > > > fSkymapOrderedData3D;

    //std::vector< std::vector< std::vector< std::vector< valarray<double>* > > > > fEnergyDensity3D;

    void ReadRadiationField2D(const std::string& filename);
    void BuildSkymapCache2D();
    void FlushSkymapCache2D();
     
    void ReadRadiationField3D(const std::string& filename);
    void BuildSkymapCache3D();
    void FlushSkymapCache3D();
  
    void ClearData();
  
  };

}

#endif
