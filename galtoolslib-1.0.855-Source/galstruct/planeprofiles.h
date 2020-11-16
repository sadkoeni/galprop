#ifndef PLANEPROFILES_H
#define PLANEPROFILES_H

#include "radialprofiles.h"
#include "spiralarms.h"
#include <memory>

#ifdef HAVE_OPENCL
#include <CL/cl2.hpp>
#endif

namespace GalacticStructure {

/** Plane profiles
 *
 * Returns a double value for a given radius and angle theta.
 * theta should be given in radians.
 */
class PlaneProfile {
   public:
      //! Default constructor
      PlaneProfile() {}
      virtual ~PlaneProfile() {}
      /** \brief XML constructor
       *
       * Reads the name of the profile and asserts that it is a plane profile
       */
      PlaneProfile( const utl::Branch &b );
      virtual double operator() ( double R, double theta ) const { std::cerr<<"Plane profile operator called"<<std::endl; throw(0); };

      //! Set parameters of the profile from variables.  Handy for doing pre-calculations.
      virtual void updateVariableValues ( const utl::Variables &vars ) { std::cerr<<"Plane profile update called"<<std::endl; throw(0); };
      
      const utl::Variables & getVariables () { return fvars; }
      
      /** \brief Add the profile to the DOMNode
       *
       * Extra attributes, if any, are added to the profile element.
       */
      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const = 0;
      
      //! Throws an exception if name is not found
      static std::unique_ptr<PlaneProfile> createProfile (const utl::Branch &b);

#ifdef HAVE_OPENCL
      /** \brief Provides an OpenCL code for the function
       *
       * The function prototype is
       * float <name> ( float R, float theta, __global const float *pars)
       * The function name is combined of a class specific prefix and a user
       * specified postfix (see below).  The pars is a global buffer, a part of
       * which the class assumes ownership of.  The user must assign the start
       * position in the buffer, the class provides information on the size 
       * required.
       *
       * \param name is on input the user specified postfix but contains on output
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
      std::string profileName;
      utl::Variables fvars;
};

/** \brief Exponential disc with hole that has eccentricity and offset
 *
 * Distribution from Frankie
 *
 * f(r, theta) = norm*exp(-(r - Rs)/R0)*(1 - exp((rhole/rh)^i)) 
 * 
 * xp = cos(phiOffset)*x + sin(phiOffset)*y
 * yp = -sin(phiOffset)*x + cos(phiOffset)*y
 * rhole = sqrt(xp*xp + e*e*yp*yp)
 * e is eccentricity, rh is hole scale, i is hole index/power, R0 is radial scale
 * phiOffset should be in radians.
 */

 class ExpDiscWithEccentricHoleProfile : public PlaneProfile {
 public:
   ExpDiscWithEccentricHoleProfile(const utl::Branch& b);
   ExpDiscWithEccentricHoleProfile(const utl::Parameters& pars, const std::string& R0name, const std::string& Rsname, const std::string& Rhname, const std::string& hiname, const std::string &phiOffsetname, const std::string& pitchAnglename, const std::string &holeEccentricityname, const std::string &normname, const std::string& Rmaxname);
   ExpDiscWithEccentricHoleProfile(double R0, double Rs, double Rh, double hi, double phiOffset, double pitchAngle, double holeEccentricity, double norm, double Rmax, const std::string& R0name, const std::string& Rsname, const std::string& Rhname, const std::string& hiname, const std::string &phiOffsetname, const std::string& pitchAnglename, const std::string &holeEccentricityname, const std::string &normname, const std::string& Rmaxname); 

   void updateVariableValues ( const utl::Variables &vars ) override;
   double operator () ( double R, double theta ) const override;

   virtual void addToDOM( xercesc::DOMNode *node, 
         const std::map<std::string,std::string> &attributes ) const override;

   void setR0(double r);
   void setRs(double r);
   void setRh(double r);
   void sethi(double i);
   void setphiOffset(double phi);
   void setpitchAngle(double pa);
   void setholeEccentricity(double e);
   void setnorm(double n);
   void setRmax(double r);

 private:

   double oneOver_R0 = 0., Rs = 0., oneOver_Rh = 0., hi = 0., cosPhiOffset = 0., sinPhiOffset = 1., cosPitchAngle = 0., sinPitchAngle = 1., eccentricity = 1., norm = 0., Rmax = 0.;
   std::string fR0name, fRsname, fRhname, fhiname, fphiOffsetname, fpitchAnglename, feccentricityname, fnormname, fRmaxname;

 };


//! Constant in theta, but possible normalization and a radial profile
class PlaneProfile1Mode : public PlaneProfile {
   public:
      /** \brief Profile defined in XML
       *
       * Requires one radial profile and a normalization variable named norm.
       * There is no need for an intent in the radial profile, the first one is used.
       */
      PlaneProfile1Mode( const utl::Branch &b );
      //! Constructor giving the normalization parameter name.  It has to be defined in the pars instance.
      PlaneProfile1Mode ( const utl::Parameters &pars, const std::string &normName, std::unique_ptr<RadialProfile> profile );

      void updateVariableValues ( const utl::Variables &vars ) override;
      double operator () ( double R, double theta ) const override;
      
      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      //! "PlaneProfile1Mode_" is added at the start of name.
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

   private:
      double norm;
      std::unique_ptr<RadialProfile> rProfile;
      std::string fnormName;

};

//! Simple cosine function with radial dependence of normalization and zero point value plus a constant in theta
class PlaneProfile2Mode : public PlaneProfile {
   public:
      /** \brief Profile defined in XML
       *
       * Requires three radial profiles with the following intent attribute added to the profile element
       * * norm1 for 0th mode normalization
       * * norm2 for 1st mode normalization
       * * zero2 the 0 point for 1st mode.
       * In addition to these profiles, two normalization variables are needed for the 
       * normalization profiles, called norm1 and norm2.
       *
       */
      PlaneProfile2Mode( const utl::Branch &b );
      /** \brief Construct the Plane profile
       *
       * Specify the radial profiles for the constant and first cosine mode.
       * Separate normalization for each mode.
       * Note: it makes little sense to have anything else than a spline profile or similar for 
       * the zero point radial profile
       */
      PlaneProfile2Mode ( const utl::Parameters &pars, 
            const std::string &norm1Name, 
            const std::string &norm2Name, 
            std::unique_ptr<RadialProfile> norm1profile,
            std::unique_ptr<RadialProfile> norm2profile, 
            std::unique_ptr<RadialProfile> zero2profile 
            );

      void updateVariableValues ( const utl::Variables &vars ) override;
      double operator () ( double R, double theta ) const override;
      
      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

   private:
      double norm1, norm2;
      std::unique_ptr<RadialProfile> fnorm1Profile, fnorm2Profile, fzero2Profile;

      std::string fnorm1Name, fnorm2Name;

};

//! constant, simple cosine, and first harmonic function, all with separate radial profiles.
class PlaneProfile3Mode : public PlaneProfile {
   public:
      /** \brief Profile defined in XML
       *
       * Requires five radial profiles with the following intent attribute added to the profile element
       * * norm1 for 0th mode normalization
       * * norm2 for 1st mode normalization
       * * zero2 the 0 point for 1st mode
       * * norm3 for 2nd mode normalization
       * * zero3 the 0 point for 2nd mode
       * In addition to these profiles, three normalization variables are needed for the 
       * normalization profiles, called norm1, norm2, and norm3.
       *
       */
      PlaneProfile3Mode( const utl::Branch &b );
      /** \brief Construct the Plane profile
       *
       * Specify the radial profiles for the constant, first cosine mode, and first harmonic
       * Separate normalization for each mode.
       * Note: it makes little sense to have anything else than a spline profile or similar for 
       * the zero point radial profile
       */
      PlaneProfile3Mode ( const utl::Parameters &pars, 
            const std::string &norm1Name, 
            const std::string &norm2Name, 
            const std::string &norm3Name, 
            std::unique_ptr<RadialProfile> norm1profile,
            std::unique_ptr<RadialProfile> norm2profile,
            std::unique_ptr<RadialProfile> zero2profile,
            std::unique_ptr<RadialProfile> norm3profile,
            std::unique_ptr<RadialProfile> zero3profile
            );

      void updateVariableValues ( const utl::Variables &vars ) override;
      double operator () ( double R, double theta ) const override;
      
      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

   private:
      double norm1, norm2, norm3;
      std::unique_ptr<RadialProfile> fnorm1Profile, fnorm2Profile, fzero2Profile;
      std::unique_ptr<RadialProfile> fnorm3Profile, fzero3Profile;
      std::string fnorm1Name, fnorm2Name, fnorm3Name;

};

//! variable number of modes, all with separate profiles for normalization and zero point (except the constant)
class PlaneProfileNMode : public PlaneProfile {
   public:
      /** \brief Profile defined in XML
       *
       * The number of modes is specified in the element numberOfModes. Then we need 
       * to specify radial profiles for all modes.  There should be one
       * normalization and zero profile per mode, expect the first one that 
       * only has a normalization.  The intent parameter of the profiles are 
       * norm<i> and zero<i>, where i is the 1 based number of the mode.  We 
       * also need normalization variables for each normalization profile, 
       * named norm<i> and zero<i>.
       *
       */
      PlaneProfileNMode( const utl::Branch &b );
      /** \brief Construct the Plane profile
       *
       * Specify the radial profiles for each mode.  The size of the normalization
       * vector should be one larger than the zero point vector.
       * Parameter name specified with a prefix and are named <prefix>_norm_ii, where
       * ii is the number of the mode, starting with 1.
       * Note: it makes little sense to have anything else than a spline profile or similar for 
       * the zero point radial profiles
       */
      PlaneProfileNMode ( const utl::Parameters &pars, 
            const std::string &prefix, 
            std::vector<std::unique_ptr<RadialProfile> > &&normProfiles,
            std::vector<std::unique_ptr<RadialProfile> > &&zeroProfiles);

      /** \brief Construct the Plane profile
       *
       * Specify the radial profiles for each mode.  The size of the normalization
       * vector should be one larger than the zero point vector.
       * Parameter names are specified with string vectors
       * Note: it makes little sense to have anything else than a spline profile or similar for 
       * the zero point radial profiles
       */
      PlaneProfileNMode ( const utl::Parameters &pars, 
            const std::vector<std::string> &normNames, 
            const std::vector<std::string> &zerNames, 
            std::vector<std::unique_ptr<RadialProfile> > &&normProfiles,
            std::vector<std::unique_ptr<RadialProfile> > &&zeroProfiles);

      void updateVariableValues ( const utl::Variables &vars ) override;
      double operator () ( double R, double theta ) const override;
      
      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      //! "PlaneProfileNMode_" is added at the start of name.
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

   private:
      std::vector<double> norm, zero;
      std::vector<std::unique_ptr<RadialProfile> > fnormProfiles, fzeroProfiles;
      std::vector<std::string> fnormNames, fzeroNames;

};


//! Single spiral arm
class PlaneProfileArm : public PlaneProfile {

   public:

      /** \brief Profile defined in XML
       *
       * The type of arm function is given in the element armfunction.
       * Available values are Exponential, Sech2, and Gaussian (case sensitive)
       * The arm structure is specified with the variables
       * norm, a, rMin, phiMin, width, and either rMax or phiExtent.
       * An exception is thrown if both rMax and phiExtent are given.
       * The radial dependence is given with a radial profile.
       * No intent is needed for the profile, the first one is used.
       */
      PlaneProfileArm( const utl::Branch &b );

      /** \brief Create the model and specify variables with a prefix
       *
       * The arm structure is specified with the armFunction
       *
       * It's radial dependence is given with radProfile
       *
       * The parameters of the arm are named: <prefix>_norm, <prefix>_a, <prefix>_rMin, <prefix>_phiMin, <prefix>_rMax,
       * <prefix>_width.  If the boolean use_rMax is set to false, phiExtent is used instead of rMax in the variables.
       */
      PlaneProfileArm ( const utl::Parameters& pars, 
            std::unique_ptr<ArmFunction> armFunction,
            std::unique_ptr<RadialProfile> radProfile, 
            const std::string &prefix,
            bool use_rMax=true);

      /** \brief Create the model and specify variables independently
       *
       * The arm structure is specified with the armFunction
       *
       * It's radial dependence is given with radProfile
       *
       * The parameters of the arm are given in varNames in the order: norm, a, rMin, phiMin,
       * width, rMax.  If the boolean use_rMax is set to false, phiExtent is used instead of rMax in the variables
       */
      PlaneProfileArm ( const utl::Parameters& pars, 
            std::unique_ptr<ArmFunction> armFunction,
            std::unique_ptr<RadialProfile> radProfile, 
            const std::vector<std::string> &varNames,
            bool use_rMax=true);

      void updateVariableValues ( const utl::Variables &vars ) override;
      double operator () ( double R, double theta ) const override;
      
      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      //! "PlaneProfileArm_" is added at the start of name.
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

   protected:
      double norm;
      std::unique_ptr<ArmFunction> arm;
      std::unique_ptr<RadialProfile> fradProfile;
      std::vector<std::string> fvarNames;
      bool use_rMax;

};

/** \brief Symmetric spiral arm structure
 *
 * N spiral arms that all have the same parameters except for the starting point.
 * phiMin for the arms is set such that the arms are distributed evenly over 2 Pi 
 * with the first at phiMin.
 */
class PlaneProfileSymmetricArms : public PlaneProfile {

   public:

      /** \brief Profile defined in XML
       *
       * The number of arms is given with the element numberOfArms
       * The type of arm function is given in the element armfunction.
       * Available values are Exponential, Sech2, and Gaussian (case sensitive)
       * The arm structure is specified with the variables
       * norm, a, rMin, phiMin, width, and either rMax or phiExtent.
       * An exception is thrown if both rMax and phiExtent are given.
       * The radial dependence is given with a radial profile.
       * No intent is needed for the profile, the first one is used.
       */
      PlaneProfileSymmetricArms( const utl::Branch &b );

      /** \brief Create the model and specify variables with a prefix
       *
       * The arm structure is specified with the armFunctions vector
       *
       * It's radial dependence is given with radProfile
       *
       * The parameters of the arms are named: <prefix>_norm, <prefix>_a, <prefix>_rMin, 
       * <prefix>_phiMin, <prefix>_rMax, <prefix>_width.  If the boolean use_rMax is set to false, 
       * phiExtent is used instead of rMax in the variables.
       */
      PlaneProfileSymmetricArms ( const utl::Parameters& pars, 
            std::vector<std::unique_ptr<ArmFunction> > &&armFunctions,
            std::unique_ptr<RadialProfile> radProfile, 
            const std::string &prefix,
            bool use_rMax=true);

      /** \brief Create the model and specify variables independently
       *
       * The arm structure is specified with the armFunctions vector
       *
       * It's radial dependence is given with radProfile
       *
       * The parameters of the arm are given in varNames in the order: norm, a, rMin, 
       * phiMin, width, rMax.  If the boolean use_rMax is set to false, phiExtent is 
       * used instead of rMax in the variables.
       */
      PlaneProfileSymmetricArms ( const utl::Parameters& pars, 
            std::vector<std::unique_ptr<ArmFunction> > &&armFunctions,
            std::unique_ptr<RadialProfile> radProfile, 
            const std::vector<std::string> &varNames,
            bool use_rMax=true);

      void updateVariableValues ( const utl::Variables &vars ) override;
      double operator () ( double R, double theta ) const override;
      
      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      //! "PlaneProfileSymmetricArms_" is added at the start of name.
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif
   protected:
      double norm;
      std::vector<std::unique_ptr<ArmFunction> > arms;
      std::unique_ptr<RadialProfile> fradProfile;
      std::vector<std::string> fvarNames;
      bool use_rMax;

};

}
#endif
