#ifndef CYLINDRICALPROFILES_H
#define CYLINDRICALPROFILES_H

#include "radialprofiles.h"
#include "planeprofiles.h"
#include "spiralarms.h"

#include <Variables.h>

#include <functional>

#ifdef HAVE_OPENCL
#include <CL/cl2.hpp>
#endif

namespace GalacticStructure {

/** Cylindrical profiles
 *
 * Returns a double value for a given radius, azimuth and z.
 */
class CylindricalProfile {
   public:
      CylindricalProfile() {}
      /** \brief XML constructor
       *
       * Reads the name of the profile and asserts that it is a cylindrical profile
       */
      CylindricalProfile( const utl::Branch &b );
      virtual ~CylindricalProfile() {}
      virtual double operator() ( double R, double theta, double z ) const { std::cerr<<"Cylindrical profile operator() called"<<std::endl; throw(0); };

      //! Set parameters of the profile from variables.  Handy for doing pre-calculations.
      virtual void updateVariableValues ( const utl::Variables &vars ) { std::cerr<<"Cylindrical profile update called"<<std::endl; throw(0); };
      
      const utl::Variables & getVariables () { return fvars; }

      const std::string & getName () { return profileName; }

      static std::unique_ptr<CylindricalProfile> createProfile( const utl::Branch &b );

      /** \brief Add the profile to the DOMNode
       */
      virtual void addToDOM( xercesc::DOMNode *node ) const = 0;
      
#ifdef HAVE_OPENCL
      /** \brief Provides an OpenCL code for the function
       *
       * The function prototype is
       * float <name> ( float R, float theta, float z, __global const float *pars)
       * The function name is combined of a class specific prefix and a user
       * specified postfix (see below).  The pars is a global buffer, a part of
       * which the class assumes ownership of.  The user must assign the start
       * position in the buffer, the class provides information on the size 
       * required.
       *
       * \param name is on in put the user specified postfix but contains on output
       * the entier function name.
       * \param parIndex is the start position in the buffer that is assigned to this
       * function.
       * \return the code in a string
       */
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const;

      /** \brief Return the number of parameters required by this class in the pars buffer
       *
       * This number should preferably not change dynamically and should be set on class 
       * construction.  There is currently no mechanism to trigger an increase in the size of
       * the buffer.
       */
      virtual size_t getOpenCLNPars() const;

      /** \brief Get the parameters that should be loaded in the buffer
       *
       * The caller is responsible for loading the parameters in the buffer.
       * The size of the output should be the same as that returned by getOpenCLNPars.
       */
      virtual std::vector<cl_float> getOpenCLPars() const;
#endif
      
   protected:
      utl::Variables fvars;
      std::string profileName;
};


/** Generic disk-like profile based 3D distribution
 *
 * center of plane is specified with a plane profile
 * density in the plane is specified with a plane profile
 * z distribution uses 2 scaling function, one for north, the 
 * other for south.  Each comes with an plane function to define
 * the scale height.
 *
 * f(r,th,z) = g(r,th) * h( |z-m(r,th)| * s(r,th) )
 *
 * where g is the plane density, h is the scale function, m defines
 * the center of the plane and s is one over scale height.  Note that h and
 * s are different for positive and negative z-m(r,th).
 *
 * TODO: Modify the code to allow many separate functions for g, i.e.
 * f(r,th,z) = h( |z-m(r,th)| * s(r,th) ) * \sum_i g_i(r,th)
 */

class GenericDiskProfile : public CylindricalProfile {

   public:
      GenericDiskProfile( const utl::Parameters &pars,
         std::unique_ptr<PlaneProfile> planeDensity,
         std::unique_ptr<PlaneProfile> planeCenter,
         std::unique_ptr<PlaneProfile> northScaleHeight,
         std::unique_ptr<PlaneProfile> southScaleHeight,
         std::unique_ptr<ScaleFunction> northScaleFunction,
         std::unique_ptr<ScaleFunction> southScaleFunction );

      /** \brief Profile defined in XML
       *
       * Requires four plane profiles with the following intent attribute added to the profile element
       * * planeDensity for the density in the plane
       * * planeCenter for the 0 point of the plane
       * * northScaleHeight for the scale height above the plane
       * * southScaleHeight for the scale height below the plane
       * In addition to these profiles, two scale functions are required
       * * <northScaleFunction> </northScaleFunction>
       * * <southScaleFunction> </southScaleFunction>
       * Available functions are Exponential, Gaussian, and Sech2
       *
       */
      GenericDiskProfile( const utl::Branch &b );

      void updateVariableValues( const utl::Variables & vars ) override;
      double operator() ( double r, double theta, double z ) const override;

      virtual void addToDOM( xercesc::DOMNode *node ) const override;
      
#ifdef HAVE_OPENCL
      //! "GenericDiskProfile_" is added at the start of name.
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

   protected:
      std::unique_ptr<PlaneProfile> density, center, northHeight, southHeight;
      std::unique_ptr<ScaleFunction> northScale, southScale;

};

/** \brief The NE2001 arm distribution
 *
 * All parameters are fixed
 */
class NE2001ArmProfile : public CylindricalProfile {
   
   public:

      /** \brief Fixed model from NE2001, no parameters allowed
       * 
       * Returns in units of cm^-3, assuming r and z are in kpc.
       */
      NE2001ArmProfile ( );

      /** \brief Profile defined in XML
       *
       * This model is fixed, nothing to be defined
       *
       */
      NE2001ArmProfile( const utl::Branch &b );

      //! No variables to update but needed
      void updateVariableValues( const utl::Variables &vars ) override;

      double operator() (double r, double theta, double z) const override;

      virtual void addToDOM( xercesc::DOMNode *node ) const override;

   protected:
      double na = 0., wa = 0., ha = 0.; //!< Possible parameters
      double wa1 = 0., wa2 = 0., wa3 = 0., wa4 = 0., wa5 = 0.; //!< Arm width
      double na1 = 0., na2 = 0., na3 = 0., na4 = 0., na5 = 0.; //!< Arm density
      double ha1 = 0., ha2 = 0., ha3 = 0., ha4 = 0., ha5 = 0.; //!< Arm scale height
      double Aa = 0.; //!< Radial scale length and cutoff
      ArmFunction arm1, arm2, arm3, arm4, arm5; //!< The 5 arms.

      void init();
};


/** \brief Wainscoat bulge
 *
 * Distribution from Frankie
 *
 * f(r,z) = norm * exp(-x^3)/x^1.8
 * 
 * with x = sqrt( rr^2 + (z/k1)^2 )/R1 and  rr = max(r,Rmin)
 */
class WainscoatBulgeProfile : public CylindricalProfile {

   public:

      /** \brief Specify parameters manually and fix them
       *
       * The variable names are given as the strings.
       */
      WainscoatBulgeProfile ( double Rmin, double R1, double k1, double norm,
           const std::string &Rminname, const std::string &R1name, const std::string &k1name, const std::string &normname );

      /** \brief Specify parameters from a parameters file
       *
       * The variable names are given as the strings.
       * The variables have to be defined in the parameter object
       */
      WainscoatBulgeProfile ( const utl::Parameters &pars, 
           const std::string &Rminname, const std::string &R1name, const std::string &k1name, const std::string &normname );

      /** \brief XML constructor
       *
       * The variable names are: Rmin, R1, k1, and norm.
       */
      WainscoatBulgeProfile ( const utl::Branch &b );

      void updateVariableValues( const utl::Variables & vars ) override;
      double operator() ( double r, double theta, double z ) const override;

      virtual void addToDOM( xercesc::DOMNode *node ) const override;

      void setRmin(double r);
      void setR1(double r);
      void setk1(double k);
      void setnorm(double n);

   private:
      double Rmin = 0., oneOver_R1 = 0., k1 = 0., norm = 0.;
      std::string fRminname, fR1name, fk1name, fnormname;
};



/** \brief Wainscoat halo
 *
 * Distribution from Frankie
 *
 * f(r,z) = norm * 10^(-3.3307*(alpha^0.25 - 1))
 *
 * with alpha = sqrt(r^2 + (z*zE)^2)/rE
 * 
 */
class WainscoatHaloProfile : public CylindricalProfile {

   public:

      /** \brief Specify parameters manually and fix them
       *
       * The variable names are given as the strings.
       */
      WainscoatHaloProfile ( double rE, double zE, double norm, 
            const std::string &rEname, const std::string &zEname, const std::string &normname );

      /** \brief Specify parameters from a parameters file
       *
       * The variable names are given as the strings.
       * The variables have to be defined in the parameter object
       */
      WainscoatHaloProfile ( const utl::Parameters &pars,
            const std::string &rEname, const std::string &zEname, const std::string &normname );

      /** \brief XML constructor
       *
       * The variable names are: rE, zE, and norm.
       */
      WainscoatHaloProfile ( const utl::Branch &b );

      virtual void addToDOM( xercesc::DOMNode *node ) const override;

      void updateVariableValues( const utl::Variables & vars ) override;
      double operator() ( double r, double theta, double z ) const override;

      void setrE(double r);
      void setzE(double z);
      void setnorm(double n);

   private:
      double oneOver_rE = 0., zE = 0., norm = 0.;
      std::string frEname, fzEname, fnormname;
};



/** \brief SDSS halo
 *
 * Distribution from Frankie
 *
 * f(r,z) = norm * ( rS/sqrt(rp^2 + zp^2) )^alpha
 *
 * with rp = max(0.1, r) and zp = z/qH
 * 
 */
class SDSSHaloProfile : public CylindricalProfile {

   public:

      /** \brief Specify parameters manually and fix them
       *
       * The variable names are given as the strings.
       */
      SDSSHaloProfile ( double rS, double qH, double alpha, double norm,
            const std::string &rSname, const std::string &qHname, const std::string &alphaname, const std::string &normname );

      /** \brief Specify parameters from a parameters file
       *
       * The variable names are given as the strings.
       * The variables have to be defined in the parameter object
       */
      SDSSHaloProfile ( const utl::Parameters &pars, 
            const std::string &rSname, const std::string &qHname, const std::string &alphaname, const std::string &normname );

      /** \brief XML constructor
       *
       * The variable names are: rS, qH, alpha, and norm.
       */
      SDSSHaloProfile ( const utl::Branch &b );

      void updateVariableValues( const utl::Variables & vars ) override;
      double operator() ( double r, double theta, double z ) const override;

      virtual void addToDOM( xercesc::DOMNode *node ) const override;

      void setrS(double r);
      void setqH(double q);
      void setalpha(double a);
      void setnorm(double n);

   private:
      double rS = 0., oneOver_qH = 0., alpha = 0., norm = 0.;
      std::string frSname, fqHname, falphaname, fnormname;
};



/** \brief Lopez Corredoira bulge 
 *
 * Distribution from Frankie
 *
 * f(r,z) = norm * exp(-t/scaleLength)
 *
 * with 
 * t = (|xp|^index + |yp/ya|^index + |z/za|^index)^(1./index)
 * and
 * xp = cos(phiOffset)*x + sin(phiOffset)*y
 * yp = -sin(phiOffset)*x + cos(phiOffset)*y
 * 
 * phiOffset should be in radians.
 */
class LopezCorredoiraBulgeProfile : public CylindricalProfile {

   public:

      /** \brief Specify parameters manually and fix them
       *
       * The variable names are given as the strings.
       */
      LopezCorredoiraBulgeProfile ( double phiOffset, double ya, double za, double index, double scaleLength, double norm, 
            const std::string &phiOffsetname, const std::string &yaname, const std::string &zaname, 
            const std::string &indexname, const std::string &scaleLengthname, const std::string &normname );

      /** \brief Specify parameters from a parameters file
       *
       * The variable names are given as the strings.
       * The variables have to be defined in the parameter object
       */
      LopezCorredoiraBulgeProfile ( const utl::Parameters &pars, 
            const std::string &phiOffsetname, const std::string &yaname, const std::string &zaname, 
            const std::string &indexname, const std::string &scaleLengthname, const std::string &normname );

      /** \brief XML constructor
       *
       * The variable names are: phiOffset, ya, za, index, scaleLength, and norm.
       */
      LopezCorredoiraBulgeProfile ( const utl::Branch &b );

      void updateVariableValues( const utl::Variables & vars ) override;
      double operator() ( double r, double theta, double z ) const override;

      virtual void addToDOM( xercesc::DOMNode *node ) const override;

      void setphiOffset(double phi);
      void setya(double q);
      void setza(double q);
      void setindex(double i);
      void setscaleLength(double l);
      void setnorm(double n);

   private:
      double cosPhiOffset = 0., sinPhiOffset = 1., oneOver_ya = 0., oneOver_za = 0., index = 0., oneOver_index = 0., oneOver_scaleLength = 0., norm = 0.;
      std::string fphiOffsetname, fyaname, fzaname, findexname, fscaleLengthname, fnormname;
};



/** \brief Freudenreich bar
 *
 * Distribution from Frankie
 *
 * f(r,z) = norm * exp( -(rs - barREnd)^2/barHEnd^2 ) * sech^2(rs)  if r > barRend
 * f(r,z) = norm * sech^2(rs)  if r <= barRend
 *
 * with 
 * rs = ( rPerp^barPara + |z/barZ|^barPara )^(1./barPara)
 * rPerp = ( |xp/barX|^barPerp + |yp/barY|^barPerp )^(1./barPerp)
 * and
 * xp = cos(phiOffset)*cos(pitchAngle)*x + sin(phiOffset)*y - cos(phiOffset)*sin(pitchAngle)*z
 * yp = -sin(phiOffset)*cos(pitchAngle)*x + cos(phiOffset)*y + sin(phiOffset)*sin(pitchAngle)*z
 * zp = sin(pitchAngle)*x + cos(pitchAngle)*z
 * phiOffset, pitchAngle should be in radians.
 */
class FreudenreichBarProfile : public CylindricalProfile {

   public:

      /** \brief Specify parameters manually and fix them
       *
       * The variable names are given as the strings.
       */
  FreudenreichBarProfile ( double phiOffset, double pitchAngle, double barX, double barY, double barZ, double barPerp, double barPara, double barREnd, double barHEnd, double norm, const std::string &phiOffsetname, const std::string& pitchAnglename, const std::string &barXname, const std::string &barYname, const std::string &barZname, const std::string &barPerpname, const std::string &barParaname, const std::string &barREndname, const std::string &barHEndname, const std::string &normname );

      /** \brief Specify parameters from a parameters file
       *
       * The variable names are given as the strings.
       * The variables have to be defined in the parameter object
       */
  FreudenreichBarProfile ( const utl::Parameters &pars, 
			   const std::string &phiOffsetname, const std::string& pitchAnglename, const std::string &barXname, const std::string &barYname, 
            const std::string &barZname, const std::string &barPerpname, const std::string &barParaname, 
            const std::string &barREndname, const std::string &barHEndname, const std::string &normname );

      /** \brief XML constructor
       *
       * The variable names are: phiOffset, barX, barY, barZ, barPerp, barPara, barREnd, barHEnd, norm
       */
      FreudenreichBarProfile ( const utl::Branch &b );

      void updateVariableValues( const utl::Variables & vars ) override;
      double operator() ( double r, double theta, double z ) const override;

      virtual void addToDOM( xercesc::DOMNode *node ) const override;

      void setphiOffset(double phi);
      void setpitchAngle(double pa);
      void setbarX(double q);
      void setbarY(double q);
      void setbarZ(double q);
      void setbarPerp(double i);
      void setbarPara(double i);
      void setbarREnd(double l);
      void setbarHEnd(double l);
      void setnorm(double n);

#ifdef HAVE_OPENCL
      //! "FreudenreichBar_" is added at the start of name.
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

   private:
      double cosPhiOffset = 0., sinPhiOffset = 1., cosPitchAngle = 0., sinPitchAngle = 1., oneOver_barX = 0., oneOver_barY = 0., oneOver_barZ = 0., 
	barPerp = 0., oneOver_barPerp = 0., barPara = 0., oneOver_barPara = 0., barREnd = 0., oneOver_barHEnd = 0., norm = 0.;
      std::string fphiOffsetname, fpitchAnglename, fbarXname, fbarYname, fbarZname, fbarPerpname, fbarParaname, fbarREndname, fbarHEndname, fnormname;
};

/** \brief Exponential disc with hole that has eccentricity and offset
 *
 * Distribution from Frankie
 *
 * 
 * xp = cos(phiOffset)*x + sin(phiOffset)*y
 * yp = -sin(phiOffset)*x + cos(phiOffset)*y
 * 
 * phiOffset should be in radians.
 */
/*class ExpDiscWithEccentricHoleProfile : public CylindricalProfile {

   public:

  ExpDiscWithEccentricHoleProfile ( double R0, double Rs, double Rh, double hi, double phiOffset, double holeEccentricity, double Zs, double Zoffset, double norm, const std::string &R0name, const std::string &Rsname, const std::string &Rhname, const std::string &hiname, const std::string &phiOffsetname, const std::string &holeEccentricityname, const std::string &Zsname, const std::string &Zoffsetname, const std::string &normname );

  ExpDiscWithEccentricHoleProfile ( const utl::Parameters &pars, 
				    const std::string &R0name, const std::string &Rsname, const std::string &Rhname, const std::string &hiname, const std::string &phiOffsetname, const std::string &holeEccentricityname, const std::string &Zsname, const std::string &Zoffsetname, const std::string &normname );

  ExpDiscWithEccentricHoleProfile ( const utl::Branch &b );

  void updateVariableValues( const utl::Variables & vars );
  double operator() ( double r, double theta, double z ) const;

  void setR0(double r);
  void setRs(double r);
  void setRh(double r);
  void sethi(double i);
  void setphiOffset(double phi);
  void setholeEccentricity(double e);
  void setZs(double z);
  void setZoffset(double z);
  void setnorm(double n);

 private:
  double oneOver_R0 = 0., Rs = 0., oneOver_Rh = 0., hi = 0., cosPhiOffset = 0., sinPhiOffset = 1., eccentricity = 1., oneOver_Zs = 0., Zoffset = 0., norm = 0.;
  std::string fR0name, fRsname, fRhname, fhiname, fphiOffsetname, feccentricityname, fZsname, fZoffsetname, fnormname;

};
*/
/** \brief Ferriere CMZ
 *
 * Distribution from Ferriére et al. 2007
 *
 * f(r,z) = norm * exp(-((sqrt(xp^2 + (Yc*yp)^2) - Xc)/Lc)^4) * exp(-(z/Hc)^2)
 *
 * with 
 * * Yc = 2.5
 * * Xc = Xmax/2
 * * Lc = Xmax/(2*ln(2)^(1/4))
 * as defulat values
 * and
 * * xp = cos(phiOffset)*(x-x0) - sin(phiOffset)*(y+y0)
 * * yp = -sin(phiOffset)*(x-x0) - cos(phiOffset)*(y+y0)
 * 
 * phiOffset should be in radians.  Note that Ferriere uses a left handed system
 * so the sign in front of y in xp and yp has been changed.
 */
class FerriereCMZProfile : public CylindricalProfile {

   public:

      /** \brief Specify parameters manually and fix them
       *
       * The variable names are given as the strings.
       * Use default values for Yc, Xc, and Lc.
       */
      FerriereCMZProfile ( double phiOffset, double Xmax, double Hc, double x0, double y0, double norm, 
            const std::string &phiOffsetname, const std::string &Xmaxname, const std::string &Hcname,
            const std::string &x0name, const std::string &y0name, const std::string &normname );

      /** \brief Specify parameters manually and fix them
       *
       * The variable names are given as the strings.
       */
      FerriereCMZProfile ( double phiOffset, double Yc, double Xc, double Lc, double Hc, double x0, double y0, double norm, 
            const std::string &phiOffsetname, const std::string &Ycname, const std::string &Xcname, const std::string &Lcname,
            const std::string &Hcname, const std::string &x0name, const std::string &y0name, const std::string &normname );

      /** \brief Specify parameters from a parameters file
       *
       * The variable names are given as the strings.
       * The variables have to be defined in the parameter object
       * Uses default values for Yc, Xc, and Lc.
       */
      FerriereCMZProfile ( const utl::Parameters &pars, 
            const std::string &phiOffsetname, const std::string &Xmaxname, const std::string &Hcname,
            const std::string &x0name, const std::string &y0name, const std::string &normname );

      /** \brief Specify parameters from a parameters file
       *
       * The variable names are given as the strings.
       * The variables have to be defined in the parameter object
       */
      FerriereCMZProfile ( const utl::Parameters &pars, 
            const std::string &phiOffsetname, const std::string &Ycname, const std::string &Xcname, const std::string &Lcname,
            const std::string &Hcname, const std::string &x0name, const std::string &y0name, const std::string &normname );

      /** \brief XML constructor
       *
       * The variable names are: phiOffset, Yc, Xc, Lc, Hc, x0, y0, and norm.
       */
      FerriereCMZProfile ( const utl::Branch &b );

      virtual void addToDOM( xercesc::DOMNode *node ) const override;

#ifdef HAVE_OPENCL
      //! "FerriereCMZ_" is added at the start of name.
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

      void updateVariableValues( const utl::Variables & vars ) override;
      double operator() ( double r, double theta, double z ) const override;

      void setphiOffset(double phi);
      void setXmax(double q);
      void setYc(double q);
      void setXc(double q);
      void setLc(double q);
      void setHc(double q);
      void setx0(double x);
      void sety0(double y);
      void setnorm(double n);

   private:
      double cosPhiOffset, sinPhiOffset, Yc, Xc, oneOver_Lc, oneOver_Hc2, x0, y0, norm;
      std::string fphiOffsetname, fXmaxname, fYcname, fXcname, fLcname, fHcname, fx0name, fy0name, fnormname;
};



/** \brief Ferriere Disk
 *
 * Distribution from Ferriére et al. 2007
 *
 * f(r,z) = norm * exp(-((sqrt(xp^2 + (Yd*yp)^2) - Xd)/Ld)^4) * exp(-(zp/Hd)^2)
 *
 * with 
 * Xd = (Xmax+Xmin)/2
 * Lc = (Xmax-Xmin)/(2*ln(2)^(1/4))
 * and
 * xp = x*cos(beta)*cos(phid) - y*(sin(alpha)*sin(beta)*cos(phid) - cos(alpha)*sin(phid)) - z*(cos(alpha)*sin(beta)*cos(phid) + sin(alpha)*sin(phid))
 * yp = -x*cos(beta)*sin(phid) + y*(sin(alpha)*sin(beta)*sin(phid) + cos(alpha)*cos(phid)) + z*(cos(alpha)*sin(beta)*sin(phid) - sin(alpha)*cos(phid))
 * zp = -x*sin(beta) + y*sin(alpha)*cos(beta) + z*cos(alpha)*cos(beta)
 * 
 * alpha, beta, and phid should be in radians.
 */
class FerriereDiskProfile : public CylindricalProfile {

   public:

      /** \brief Specify parameters manually and fix them
       *
       * The variable names are given as the strings.
       * Yd set to the default value of 3.1
       */
      FerriereDiskProfile ( double alpha, double beta, double phid, double Xmax, double Xmin, double Hd, double norm, 
            const std::string &alphaname, const std::string &betaname, const std::string &phidname, 
            const std::string &Xmaxname, const std::string &Xminname, const std::string &Hdname, const std::string &normname );

      /** \brief Specify parameters manually and fix them
       *
       * The variable names are given as the strings.
       */
      FerriereDiskProfile ( double alpha, double beta, double phid, double Xmax, double Xmin, double Yd, double Hd, double norm, 
            const std::string &alphaname, const std::string &betaname, const std::string &phidname, const std::string &Xmaxname, 
            const std::string &Xminname, const std::string &Ydname, const std::string &Hdname, const std::string &normname );

      /** \brief Specify parameters from a parameters file
       *
       * The variable names are given as the strings.
       * The variables have to be defined in the parameter object
       * Yd set to the default value of 3.1
       */
      FerriereDiskProfile ( const utl::Parameters &pars, 
            const std::string &alphaname, const std::string &betaname, const std::string &phidname, 
            const std::string &Xmaxname, const std::string &Xminname, const std::string &Hdname, const std::string &normname );

      /** \brief Specify parameters from a parameters file
       *
       * The variable names are given as the strings.
       * The variables have to be defined in the parameter object
       */
      FerriereDiskProfile ( const utl::Parameters &pars, 
            const std::string &alphaname, const std::string &betaname, const std::string &phidname, const std::string &Xmaxname, 
            const std::string &Xminname, const std::string &Ydname, const std::string &Hdname, const std::string &normname );

      /** \brief XML constructor
       *
       * The variable names are: alpha, beta, phid, Xmax, Xmin, Yd, Hd, norm
       */
      FerriereDiskProfile ( const utl::Branch &b );

      virtual void addToDOM( xercesc::DOMNode *node ) const override;

#ifdef HAVE_OPENCL
      //! "FerriereDisk_" is added at the start of name.
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

      void updateVariableValues( const utl::Variables & vars ) override;
      double operator() ( double r, double theta, double z ) const override;

      void setalpha(double phi);
      void setbeta(double phi);
      void setphid(double phi);
      void setXmax(double q);
      void setXmin(double q);
      void setYd(double q);
      void setHd(double q);
      void setnorm(double n);

   private:
      void setXdLd();
      double cosAlpha, sinAlpha, cosBeta, sinBeta, cosPhid, sinPhid, Xd, oneOver_Ld, Yd, oneOver_Hd2, Xmax, Xmin, norm;
      std::string falphaname, fbetaname, fphidname, fXmaxname, fXminname, fYdname, fHdname, fnormname;
};

/** \brief Generic bar
 *
 * Bar profile that should encompass most of the ones listed above with the right parameters.
 *
 * f(rr) = norm * exp( -rr^ei ) * rr^pi
 *
 * where
 * rr = ( (rp/rSc)^zi + ((|zp-z0|)/zSc)^zi )^(1./zi)
 * with
 * rp = ( (|yp|/ySc)^ri + (|xp-x0|)^ri )^(1./ri)
 * and
 * xp = cos(phi)*cos(alpha)*x + sin(phi)*cos(alpha)*y + sin(alpha)*z
 * yp = -sin(phi)*x + cos(phi)*y
 * zp = -cos(phi)*sin(alpha)*x - sin(phi)*sin(alpha)*y + cos(alpha)*z
 *
 * phi and alpha should be in radians, ySc, rSc, zSc, x0, and z0 in kpc and other parameters are unitless.
 *
 * with zi and ri == 2 we can interpret rr and rp as distance from center.  Larger values
 * give boxier bars, smaller values give star like bars.
 */
class GenericBarProfile : public CylindricalProfile {

   public:

      /** \brief Specify parameters manually and fix them
       *
       * The variable names are given as the strings.
       */
      GenericBarProfile ( double phi, double alpha, 
            double ySc, double rSc, double zSc, 
            double ri, double zi, double x0, double z0,
            double ei, double pi, double norm, 
            const std::string &phiname, const std::string& alphaname, 
            const std::string &yScname, const std::string &rScname, const std::string &zScname, 
            const std::string &riname, const std::string &ziname, const std::string &x0name, const std::string &z0name,
            const std::string &einame, const std::string &piname, const std::string &normname );

      /** \brief Specify parameters from a parameters file
       *
       * The variable names are given as the strings.
       * The variables have to be defined in the parameter object
       */
      GenericBarProfile ( const utl::Parameters &pars, 
            const std::string &phiname, const std::string& alphaname, 
            const std::string &yScname, const std::string &rScname, const std::string &zScname, 
            const std::string &riname, const std::string &ziname, const std::string &x0name, const std::string &z0name,
            const std::string &einame, const std::string &piname, const std::string &normname );

      /** \brief XML constructor
       *
       * The variable names are: phi, alpha, xSc, ySc, zSc, ri, zi, x0, z0, ei, pi, norm
       */
      GenericBarProfile ( const utl::Branch &b );

      void updateVariableValues( const utl::Variables & vars ) override;
      double operator() ( double r, double theta, double z ) const override;

      virtual void addToDOM( xercesc::DOMNode *node ) const override;

      void setphi(double phi);
      void setalpha(double pa);
      void setySc(double q);
      void setrSc(double q);
      void setzSc(double q);
      void setri(double i);
      void setzi(double i);
      void setx0(double l);
      void setz0(double l);
      void setei(double l);
      void setpi(double l);
      void setnorm(double n);

#ifdef HAVE_OPENCL
      //! "GenericBar_" is added at the start of name.
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

   private:
      double fphi, fcosAlpha, fsinAlpha, fOO_ySc, fOO_rSc, fOO_zSc, fri, fzi, fOO_ri, fOO_zi, fx0, fz0, fei, fpi, fnorm;
      std::string fphiname, falphaname, fyScname, frScname, fzScname, friname, fziname, fx0name, fz0name, feiname, fpiname, fnormname;
};


/** \brief nHI profile from GALPROP branch
 */
class GALPROPnHIProfile : public CylindricalProfile {

   public:

      GALPROPnHIProfile ( ) {}

      /** \brief XML constructor
       *
       * The radius of the sun can be set with the element Rsun
       */
      GALPROPnHIProfile ( const utl::Branch &b );

      void updateVariableValues( const utl::Variables & vars ) override {}
      double operator() ( double r, double theta, double z ) const override;

      virtual void addToDOM( xercesc::DOMNode *node ) const override;

   private:
      double Rsun;
};

/** \brief nH2 profile from GALPROP branch
 */
class GALPROPnH2Profile : public CylindricalProfile {

   public:

      GALPROPnH2Profile ( ) {}

      /** \brief XML constructor
       *
       * The radius of the sun can be set with the element Rsun
       */
      GALPROPnH2Profile ( const utl::Branch &b );

      void updateVariableValues( const utl::Variables & vars ) override {}
      double operator() ( double r, double theta, double z ) const override;

      virtual void addToDOM( xercesc::DOMNode *node ) const override;

   private:
      double Rsun;
};



}

#endif
