//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * source_distribution.cc *                      galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include "SourceClass_Compatibility.h"
#include "constants.h"

#include "galprop_internal.h"
#include <cylindricalprofiles.h>
#include <Parameters.h>

#include <cmath>
#include <vector>

using namespace std;

#include <StellarGeometryDistributions.h>

//static vector<rf::ArmData> gArmData;

//static utl::NE2001ArmFunction gArmFunction;

// Profiles for the galstruct library.  This is a quick and dirty hack to get it working.
static std::vector<std::unique_ptr<GalacticStructure::CylindricalProfile> > gsProfiles;

// cosmic ray source distribution: x,y,z in kpc

double SourceClass_Compatibility::source_distribution (const double x, 
      const double y, 
      const double z, 
      int srcModel, 
      int n_spatial_dimensions,
      const std::vector<double> &parameters,
      const std::vector<double> &source_values,
      const std::vector<double> &source_radius) const {
   
  const double r = sqrt(x*x + y*y), phi = atan2(y, x);
  double alpha,beta,rmax,ro,rc,zscale,Value;
  double rm,rs,rconst,roff;                               //AWS20091001,20091008

  double result = 0;
  
  if (1 == srcModel)  //arbitrary parametrization
    {         
      alpha= parameters[1];
      beta = parameters[2]; 
      rmax = parameters[3]; 
      rconst=parameters[4]; // for rconst<r<rmax, set to value at rconst.     in kpc from GC  AWS20091008
      roff  =parameters[5]; // replace r with x=r+roff and r0 with x0=r0+roff in the equations (Yusifov & Kucuc 2004).   in kpc  GJ20100223
      
      ro=8.5;
      zscale=parameters[0];
      double xr      = r      + roff;
      double xro     = ro     + roff;
      double xrconst = rconst + roff;
      
      result =pow(xr     /xro,alpha) * exp(-beta*(xr     -xro)/xro) * exp(-fabs(z)/zscale);
      if (r >= rconst)result =pow(xrconst/xro,alpha) * exp(-beta*(xrconst-xro)/xro) * exp(-fabs(z)/zscale);//AWS20091008
      
      if (r >= rmax)  result = 0.0;
    }
  
  if (2 == srcModel)  //SNR:  Case and Bhattacharya 3rd Compton p437
    {
      alpha=1.69;
      beta=3.33;
      ro=8.5;
      zscale=parameters[0];
      result =pow(r/ro,alpha) * exp(-beta*(r-ro)/ro) * exp(-fabs(z)/zscale);
    }
  
  if (3 == srcModel)  //pulsars: Taylor, Manchester and Lyne 1993 ApJSupp 88,259
    {
      rc=3.5;
      zscale=parameters[0];
      result =cosh(ro/rc)/cosh(r/rc)  * exp(-fabs(z)/zscale);
    }    
  
  
  if (5 == srcModel) //Strong and Mattox gamma-ray distribution E> 100 MeV
    {
      zscale=parameters[0];
      Value=19.3;
      if (r >  4.0)  Value=21.9;
      if (r >  8.0)  Value=15.8;
      if (r > 10.0)  Value=18.3;
      if (r > 12.0)  Value=13.3;
      if (r > 15.0)  Value= 7.4;
      result = Value * exp(-fabs(z)/zscale);
    }
  
  if (6 == srcModel)  //Strong and Mattox gamma-ray distribution E> 100 MeV R<15 kpc, zero beyond
    {
      zscale=parameters[0];
      Value=19.3;
      if (r >  4.0)  Value=21.9;
      if (r >  8.0)  Value=15.8;
      if (r > 10.0)  Value=18.3;
      if (r > 12.0)  Value=13.3;
      if (r > 15.0)  Value= 0.0;
      result =  Value  * exp(-fabs(z)/zscale);
    }         
  
  
  if (7 == srcModel)  // Gaussian                AWS20091001
    {         
      rm   = parameters[1]; // mean  of Gaussian in kpc from GC
      rs   = parameters[2]; // sigma of Gaussian in kpc
      rmax = parameters[3]; // cutoff radius     in kpc from GC
      rconst=parameters[4]; // for rconst<r<rmax, set to value at rconst.     in kpc from GC  AWS20091008
      
      
      zscale=parameters[0];
      
      result =   exp(-(r     -rm)*(r     -rm)/(2.*rs*rs)) * exp(-fabs(z)/zscale);
      if (r >= rconst)result =   exp(-(rconst-rm)*(rconst-rm)/(2.*rs*rs)) * exp(-fabs(z)/zscale); //AWS20091008
      if (r >= rmax  )result =   0.0;
    }
  
  if (8 == srcModel) // Linear interpolation from tabulated values
    {
      //Find correct bin for interpolation (extrapolation if needed)
      if (r <= source_radius[0]) {
	result = source_values[0];
      } else if ( r >= source_radius[source_values.size()-1]) {
	result = source_values[source_values.size()-1];
      } else {
	int i = 0;
	while (r > source_radius[i])
	  i++;
	result = source_values[i-1] + 
	  (source_values[i]-source_values[i-1])/(source_radius[i]-source_radius[i-1])*
	  (r-source_radius[i-1]);
      }
      zscale=parameters[0];
      result *= exp(-fabs(z)/zscale);
    }
  
  if (9 == srcModel) //Total gas distribution
    {
      //Simply sum up the analytial gas functions
      result = nHI3D(x,y,z) + 2*nH23D(x,y,z);
    }
  
  if (10 == srcModel) //Only CO
    {
      result = nH23D(x,y,z);
    }
  
  if (11 == srcModel) {
    
    result = nHII3D(x,y,z);
    
  }
  
  if (12 == srcModel) {

    // NE2001 HII distribution for arms and GC -- this needs a lot of fixing because of how the arm functions are implemented.

    if (r >= rmax)  result = 0.0;

    //initModelM2();

    const double armContribution = 0;//(*WarpArmModelFunction)(x,y,z,vec3(0,0,0))[0];

    /*const double armRScale = 4., armZScale = 0.09; // kpc -- radial and scale height of O stars
    
    if (0 == gArmData.size()) {

      rf::ArmData arm1(0.50, 4.25, 3.48, 15., 0., 6., 1.);
      rf::ArmData arm2(1.2, 4.25, 3.48, 15., 3.141, 6., 1.5);
      rf::ArmData arm3(1.3, 4.89, 4.9, 15., 2.525, 6., 1.);
      rf::ArmData arm4(1.0, 4.89, 3.76, 15., 4.24, 6., 0.8);
      rf::ArmData arm5(0.25, 4.57, 8.1, 12., 5.847, 0.55, 1.);

      gArmData.push_back(arm1);
      gArmData.push_back(arm2);
      gArmData.push_back(arm3);
      gArmData.push_back(arm4);
      gArmData.push_back(arm5);

    }

    double armContribution = 0;
    */
    /*    for (vector<rf::ArmData>::const_iterator aIt = gArmData.begin();
	 aIt != gArmData.end(); ++aIt) {

      double result = 0;
      
      if (2 == n_spatial_dimensions) {

	const double dPhi = 20.*utl::kConvertDegreesToRadians;
	
	const int phiSteps = int(utl::kTwoPi/dPhi) + 1;
	
	for (int iPhi = 0; iPhi < phiSteps; ++iPhi) {
	  
	  const double phi = iPhi*dPhi;
	  
	  result += gArmFunction(r, armRScale, Rsun, phi, z, armZScale, *aIt);
	  
	}

	result /= phiSteps;
	
      } else 
	result = gArmFunction(r, armRScale, Rsun, phi, z, armZScale, *aIt);
      	
      armContribution += result;

    }
    */

    const double nGC = 10., rGC = 0.145, rGC2 = rGC*rGC, hGC = 0.026, hGC2 = hGC*hGC, xGC = -0.01, yGC = 0., zGC = -0.020; // For the Galactic centre component

    const double rPerp2 = (x - xGC)*(x - xGC) + (y - yGC)*(y - yGC);

    const double gcContribution = nGC*exp(-(rPerp2/rGC2 + (z - zGC)*(z - zGC)/hGC2));
   
    result = gcContribution + armContribution;
    
  }

  if (13 == srcModel) {

    // Corresponds to bulge + thin disc from the ISRF model
 
    const double discDensity = parameters[0];
    const double rScale = parameters[1];
    const double zScale = parameters[2];
    const double bulgeDensity = parameters[3];
    const double bulgeA = parameters[4];
    const double bulgeB = parameters[5];
    const double bulgeScaleLength = parameters[6];
    const double bulgeIndex = parameters[7];
    const double phiOffset = parameters[8]*utl::kConvertDegreesToRadians;
    
    const double rS = 8.5;
    const double hole = 1.3, holeIndex = 1.7;
    
    const double discContribution = 
      discDensity*exp(-(r - Rsun)/rScale)* 
      (1. - exp(-pow(r/hole, holeIndex)))*exp(-fabs(z)/zScale);
  
    const double phi = atan2(y, x);

    const double bulgeContribution =
      (2 == n_spatial_dimensions ?
       utl::LopezCorredoiraBulgeAverage(r, z, phiOffset, 1., bulgeScaleLength, bulgeA, bulgeB, bulgeIndex, utl::kPi/8) : 
       utl::LopezCorredoiraBulge(r, phi, z, phiOffset, 1., bulgeScaleLength, bulgeA, bulgeB, bulgeIndex));
     
    result = discContribution + bulgeContribution;
    
    //cout << x << " " << y << " " << r << " " << z << " " << discContribution << " " << bulgeContribution << " " << result << endl;

  }

  if (14 == srcModel) {

    // Corresponds to parameterised model + ellipsoid centred on GC.

    const double bulgeDensity = parameters[0];
    const double bulgeA = parameters[1];
    const double bulgeB = parameters[2];
    const double bulgeScaleLength = parameters[3];
    const double bulgeIndex = 2.;
    const double phiOffset = 0.0*utl::kConvertDegreesToRadians;; // deg -> radians
    
    const double phi = atan2(y, x);

    const double bulge =
      bulgeDensity*(2 == n_spatial_dimensions ?
       utl::LopezCorredoiraBulgeAverage(r, z, phiOffset, 1., bulgeScaleLength, bulgeA, bulgeB, bulgeIndex, utl::kPi/8) : 
       utl::LopezCorredoiraBulge(r, phi, z, phiOffset, 1., bulgeScaleLength, bulgeA, bulgeB, bulgeIndex));
 
    const double alpha = parameters[4];
    const double beta = parameters[5];
    const double rMax = parameters[6];
    const double rConst = parameters[7];
    const double rOffset = parameters[8];

    const double zScale = parameters[9];

    const double xr = r + rOffset;
    const double xS = Rsun + rOffset;
    const double xConst = rConst + rOffset;
      
    double parameterised = pow(xr/xS, alpha)*exp(-beta*(xr - xS)/xS)*exp(-fabs(z)/zScale);

    if (r >= rConst)
      parameterised = pow(xConst/xS, alpha)*exp(-beta*(xConst - xS)/xS)*exp(-fabs(z)/zScale);

    if (r >= rMax)
      parameterised = 0;

    result = bulge + parameterised;

    //cout << r << " " << z << " " << bulge << " " << parameterised << endl;
    
  }

  if (15 == srcModel) {

     // XML profiles using libgalstruct

     if ( gsProfiles.size() == 0 ) {

        //Read in the xml file, full path should be given
        utl::Reader xmlreader(source_xmlFile);

        // Add all cylindrical profiles in the top branch to the model
        for ( utl::Branch rp = xmlreader.GetTopBranch().GetFirstChild(); rp; rp = rp.GetNextSibling()) {
           //Make sure it is the correct type before adding
           if ( rp.GetBranchNameString() == "CylindricalProfile" ) {
              gsProfiles.push_back( GalacticStructure::CylindricalProfile::createProfile(rp));
           }
        }

        if (gsProfiles.size() == 0) {
           ERROR("No source profiles in source_xmlFile.");
           throw(std::runtime_error("Fix your GALDEF file"));
        }
     }

     result = 0;
     for (size_t i(0); i < gsProfiles.size(); ++i)
        result += (*gsProfiles[i])(r,phi,z);

  }

  /*
  if (-2 == srcModel) {

     //Use model 1 as the base
      alpha= parameters[1];
      beta = parameters[2]; 
      rmax = parameters[3]; 
      rconst=parameters[4]; // for rconst<r<rmax, set to value at rconst.     in kpc from GC  AWS20091008
      roff  =parameters[5]; // replace r with x=r+roff and r0 with x0=r0+roff in the equations (Yusifov & Kucuc 2004).   in kpc  GJ20100223
      
      ro=8.5;
      zscale=parameters[0];
      double xr      = r      + roff;
      double xro     = ro     + roff;
      double xrconst = rconst + roff;
      
      result =pow(xr/xro,alpha) * exp(-beta*(xr-xro)/xro) * exp(-fabs(z)/zscale);
      if (r >= rconst)result =pow(xrconst/xro,alpha) * exp(-beta*(xrconst-xro)/xro) * exp(-fabs(z)/zscale);//AWS20091008
      
      if (r >= rmax)  result = 0.0;

      initModelM2();

      result *= (*WarpArmModelFunction)(x,y,z,vec3(0,0,0))[0];

  }
  */

  if (3 == n_spatial_dimensions) {
     
    for (int i_cr_source=0; i_cr_source<n_cr_sources; i_cr_source++) {

      const double r2 = 
	pow(cr_source_x[i_cr_source] - x, 2.) + 
	pow(cr_source_y[i_cr_source] - y, 2.) +
	pow(cr_source_z[i_cr_source] - z, 2.);

      const double src = cr_source_L[i_cr_source]*exp(-r2/(2.*pow(cr_source_w[i_cr_source], 2.)));

      result += src;
// cout<<" source_distribution: x y z cr source r2 L s:"<<x<<" "<<y<<" "<<z<<" "<<i_cr_source<<" "<<r2<<" "<<s<<endl;
    }

  }   //  n_spatial_dimensions==3

   return result;
}
