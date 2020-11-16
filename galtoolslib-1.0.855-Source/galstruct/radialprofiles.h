#ifndef RADIALPROFILES_H
#define RADIALPROFILES_H

#include <Variables.h>
#include <Reader.h>
#include <Singleton.h>

#include <vector>
#include <memory>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#ifdef HAVE_OPENCL
#include <CL/cl2.hpp>
#endif


namespace GalacticStructure {

double sech(double x);
double sech2(double x);

/** \brief Base class for defining funcional dependence in radius
 *
 * Simply a function that takes a positive double and returns a double.
 * It should always be possible to set the parameters of the profile
 * with a utl::Variables instance.
 *
 * Each profile should have a constructor taking an xml Reader Branch.
 * Variables are declared with a
 * <variable name=string value=double step=double min=double max=double varName=string/>
 * element.  varName is optional and if not given the name of the varible in 
 * the Variables class should be <profileName>_<variableName> where 
 * variableName is from the variable name attribute and profileName is the 
 * data in a <name> element that is required.
 *
 * Profiles should document the available variable names.
 */
class RadialProfile
{
   public:
      RadialProfile() {}
      //! Asserts that the branch name is RadialProfile and reads the name
      RadialProfile( const utl::Branch &b );
      virtual ~RadialProfile() {}

      virtual double operator() ( double R ) const = 0;

      //! Set parameters of the profile from variables.  Handy for doing pre-calculations.
      virtual void updateVariableValues ( const utl::Variables &vars ) = 0;

      const utl::Variables & getVariables () { return fvars; }

      /** \brief Add the profile to the DOMNode
       *
       * Extra attributes, if any, are added to the profile element.
       */
      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const = 0;
      
      //! Throws an exception if name is not found
      static std::unique_ptr<RadialProfile> createProfile (const utl::Branch &b);

#ifdef HAVE_OPENCL
      /** \brief Provides an OpenCL code for the function
       *
       * The function prototype is
       * float <name> ( float R, __global const float *pars)
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

      utl::Variables fvars;
      std::string profileName;

      /** \brief Helper function to create XML profiles.
       *
       * Adds the elements with the name and text given in the map.
       * Also adds the variables with the names and ids.
       */
      void addToDOMhelper( xercesc::DOMNode *node, 
            const std::string &type,
            const std::map<std::string,std::string> &attributes,
            const std::map<std::string,std::string> &elements, 
            const std::vector<std::string> &varNames,
            const std::vector<std::string> &varIds
            ) const;
};

//! Constant in radius
class ConstantRadialProfile : public RadialProfile
{
   public:
      ConstantRadialProfile() {}
      //! Really does nothing except asserting the type.
      ConstantRadialProfile( const utl::Branch &b );
      //virtual ~ConstantRadialProfile() {}
      virtual double operator () ( double R ) const override {return 1;};
      virtual void updateVariableValues ( const utl::Variables &vars ) override{}

      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif
};


/** \brief Implements the canonical radial pulsar shape function exponential multiplied with a powerlaw.
 *
 * The function is
 *  f(R) = norm * ((R+Roff)/R0)**alpha * exp( -beta * (R-R0)/(R0+Roff) )
 * where R0, Roff, alpha, and beta are parameters.  norm = (alpha/beta)**alpha * exp(beta-alpha) is to 
 * ensure the peak value is consistently at 1.
 */
class PulsarRadialProfile : public RadialProfile {
   public:
      /** \brief XML constructor
       *
       * The variable names are: R0, Roff, alpha, and beta.
       * Their meaning should be clear.
       */
      PulsarRadialProfile ( const utl::Branch &b );

      /** \brief Constructor with only a prefix to variable names
       *
       * The variable names are <prefix>_[name] where name is one of 
       * alpha, beta, Roff, and R0.  Those variables have to be set in the 
       * pars instance.
       */ 
      PulsarRadialProfile ( const utl::Parameters &pars, const std::string &prefix );

      /** \brief Constructor with only a prefix to variable names
       *
       * The variable names are <prefix>_[name] where name is one of 
       * alpha, beta, Roff, and R0.  The variables will be fixed.
       */ 
      PulsarRadialProfile ( double alpha, double beta, double Roff, double R0, const std::string &prefix );

      /** \brief Constructor with each variable name specified
       *
       * The variable names are given in the string vector in the order
       * alpha, beta, Roff, and R0.  Those variables have to be set in
       * the pars instance.
       */
      PulsarRadialProfile ( const utl::Parameters &pars, const std::vector<std::string> &varNames );

      virtual double operator () ( double R ) const override;
      virtual void updateVariableValues ( const utl::Variables &vars ) override;

      void setAlpha(double a);
      void setBeta(double b);
      void setRoff(double r);
      void setR0(double r);

      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif
   private:
      //Only need to store x0 because R0 is never used without Roff.
      double alpha, beta, Roff, x0, oneOver_x0, norm;
      std::vector<std::string> fvarNames;

};

/** \brief Implements the exponential radial function
 *
 * The function is
 *  f(R) = exp(-(R-Rs)/R0)
 * where R0 is the scaling parameter and Rs the radial shift or normalization
 */
class ExponentialRadialProfile : public RadialProfile {
   public:
      /** \brief XML constructor
       *
       * The variable names are: R0
       */
      ExponentialRadialProfile ( const utl::Branch &b );

      //! \brief Constructor giving the scaling parameter name.  It has to be defined in the pars instance.
      ExponentialRadialProfile ( const utl::Parameters &pars, const std::string &R0name, const std::string& RSname );

      //! \brief Constructor giving the scaling parameter name.  
      ExponentialRadialProfile ( double R0, double rs, const std::string &R0name, const std::string& RSname );

      virtual double operator () ( double R ) const override;
      virtual void updateVariableValues ( const utl::Variables &vars ) override;

      void setR0(double r);
      void setRS(double r);

      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif
   private:
      double oneOver_R0, rs;
      std::string fR0name, fRSname;

};

/** \brief Implements the Freudenreich (1998) Eq. 2 radial profile
 *
 * The function is
 *  f(R) = c1*u + c2*u^2 + c3*u^3 with u = R - Rw
 * where constants c1, c2, c3, and warp radius Rw are specified
 */
class FreudenreichWarpRadialProfile : public RadialProfile {
   public:
      /** \brief XML constructor
       *
       * The variable names are: R0
       */
      FreudenreichWarpRadialProfile ( const utl::Branch &b );

      //! \brief Constructor giving the scaling parameter name.  It has to be defined in the pars instance.
      FreudenreichWarpRadialProfile ( const utl::Parameters &pars, const std::string &Rwname, const std::string& c1name, const std::string& c2name, const std::string& c3name );

      //! \brief Constructor giving the scaling parameter name.  
      FreudenreichWarpRadialProfile ( double Rw, double c1, double c2, double c3, const std::string &Rwname, const std::string& c1name, const std::string& c2name, const std::string& c3name );

      virtual double operator () ( double R ) const override;
      virtual void updateVariableValues ( const utl::Variables &vars ) override;

      void setc1(double);
      void setc2(double);
      void setc3(double);
      void setRw(double);

      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif
   private:
      double Rw, c1, c2, c3;
      std::string fRwname, fc1name, fc2name, fc3name;

};

/** \brief Implements a cut off polynomial
 *
 * The function is
 *  f(R) = u > 0 ? sum_{i=0}^N ci*u^i : c0 with u = R - Rc
 * where constants ci, Rc, and N are specified.
 */
class CutOffPolynomialRadialProfile : public RadialProfile {
   public:
      /** \brief XML constructor
       *
       * Needs the child attribute <polynomialDegree>N</polynomialDegree>
       *
       * The variable names are: Rc, and c<xx> where <xx> is
       * a number from 0 to N of size 2, zero padded. This limits N to 99.
       */
      CutOffPolynomialRadialProfile ( const utl::Branch &b );

      //! \brief Constructor giving the variable names that should be defined in the pars object.  The degree is specified with the size of cinames.
      CutOffPolynomialRadialProfile ( const utl::Parameters &pars, const std::string& Rcname, const std::vector<std::string> &cinames );

      //! \brief Constructor giving the scaling parameter name.  
      CutOffPolynomialRadialProfile ( double Rc, std::vector<double> cis, const std::string& Rcname, const std::vector<std::string>& cinames );

      virtual double operator () ( double R ) const override;
      virtual void updateVariableValues ( const utl::Variables &vars ) override;

      void setRc(double);
      void setci(std::vector<double>);

      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif
   private:
      double Rcut;
      std::vector<double> ci;
      std::string fRcname;
      std::vector<std::string> ciname;

};

/** \brief Exponential function with hole in the center
 *
 * The function is
 * f(R) = exp( -(R - Rs)/R0 ) * (1 - exp( -pow(R/Rh,hi)))
 * where R0 is the scaling parameter, Rs is the normalization radius, 
 * Rh is the hole radius and hi controls the sharpness of the hole
 */
class ExpHoleRadialProfile : public RadialProfile {
   public:
      /** \brief XML constructor
       *
       * The variable names are R0, Rs, Rh, and hi
       */
      ExpHoleRadialProfile ( const utl::Branch &b );

      /** \brief Constructor given variable names
       *
       * The variable names are given with the four strings.  Those 
       * have to be set in the pars instance.
       */ 
      ExpHoleRadialProfile ( const utl::Parameters &pars, 
			     const std::string &R0name, 
			     const std::string &Rsname, 
			     const std::string &Rhname, 
			     const std::string &hiname);

      /** \brief Constructor given variable names and values
       *
       * The variable names are given with the four strings.
       */ 
      ExpHoleRadialProfile ( double R0, double Rs, double Rh, double hi,
			     const std::string &R0name, const std::string &Rsname, const std::string &Rhname, const std::string &hiname);

      virtual double operator () ( double R ) const override;
      virtual void updateVariableValues ( const utl::Variables &vars ) override;

      void setR0(double r);
      void setRs(double r);
      void setRh(double r);
      void sethi(double i);

      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

   private:
      double oneOver_R0, Rs, oneOver_Rh, hi;
      std::string fR0name, fRsname, fRhname, fhiname;

};

/** \brief Exponential function with gaussian hole in the centre
 *
 * The function is
 * f(R) = f0*exp(-0.5*((R - mu0)/sigma0)^2) for R < Rsmooth
 *           exp(-R/hR)                     for R >= Rsmooth
 *
 * where Rsmooth = sigma0^2/hR + mu0
 *       f0      = exp(-Rsmooth/hR)/(exp(-0.5*((Rsmooth - mu0)/sigma0)^2))
 *
 * and these come from matching the piecewise function at Rsmooth ensuring
 * the continuity. From Robitaille et al. A&A 545, 39 (2012) for their modified
 * dust hole (they use parameters mu0 = 4.5 kpc and sigma0 = 1 kpc for 
 * hR = 3.5 kpc).
 */

class ExpGaussianHoleRadialProfile : public RadialProfile {
   public:
      /** \brief XML constructor
       *
       * The variable names are Sigma0, Mu0, Rh
       */
      ExpGaussianHoleRadialProfile ( const utl::Branch &b );

      /** \brief Constructor given variable names
       *
       * The variable names are given with the four strings.  Those 
       * have to be set in the pars instance.
       */ 
      ExpGaussianHoleRadialProfile ( const utl::Parameters &pars, 
				     const std::string &Sigma0name, 
				     const std::string &Mu0name, 
				     const std::string &Rhname);

      /** \brief Constructor given variable names and values
       *
       * The variable names are given with the four strings.
       */ 
      ExpGaussianHoleRadialProfile ( double Sigma0, double Mu0, double Rh,
				     const std::string &Sigma0name, const std::string &Mu0name, const std::string &Rhname);

      virtual double operator () ( double R ) const override;
      virtual void updateVariableValues ( const utl::Variables &vars ) override;

      void setMu0(double mu);
      void setSigma0(double sigma);
      void setRh(double r);

      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

   private:
      double mu0 = 0., sigma0 = 0., oneOver_Rh = 0., oneOver_sigma0 = 1., rSmooth = 0., f0 = 0.;
      std::string fMu0name, fSigma0name, fRhname;

};

/** \brief Implements the sech2 radial function
 *
 * The function is
 *  f(R) = sech2(R/R0)
 * where R0 is the scaling parameter.
 */
class Sech2RadialProfile : public RadialProfile {
   public:
      /** \brief XML constructor
       *
       * The variable names are: R0
       */
      Sech2RadialProfile ( const utl::Branch &b );

      //! \brief Constructor giving the scaling parameter name.  It has to be defined in the pars instance.
      Sech2RadialProfile ( const utl::Parameters &pars, const std::string &R0name );

      //! \brief Constructor giving the scaling parameter name.
      Sech2RadialProfile ( double R0, const std::string &R0name );

      virtual double operator () ( double R ) const override;
      virtual void updateVariableValues ( const utl::Variables &vars ) override;

      void setR0(double r);

      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

   private:
      double oneOver_R0;
      std::string fR0name;

};

/** \brief Implements the Gaussian shape for a radial profile
 *
 * The function is
 *  f(R) = exp (-(R+Roff)**2/R0**2)
 * where R0 and Roff are parameters.
 */
class GaussianRadialProfile : public RadialProfile {
   public:
      /** \brief XML constructor
       *
       * The variable names are: R0 and Roff
       */
      GaussianRadialProfile ( const utl::Branch &b );

      /** \brief Constructor given variable names
       *
       * The variable names are given with the two strings.  Those 
       * have to be set in the pars instance.
       */ 
      GaussianRadialProfile ( const utl::Parameters &pars, const std::string &R0name, const std::string &Roffname );

      /** \brief Constructor given variable names and values
       *
       * The variable names are given with the two strings.
       */ 
      GaussianRadialProfile ( double R0, double Roff, const std::string &R0name, const std::string &Roffname );

      virtual double operator () ( double R ) const override;
      virtual void updateVariableValues ( const utl::Variables &vars ) override;

      void setR0(double r);
      void setRoff(double r);

      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

   private:
      double oneOver_R02, Roff;
      std::string fR0name, fRoffname;

};

/** \brief Implements multiple Gaussian shape for a radial profile
 *
 * The function is
 *  f(R) = sum_i n_i exp (-(R+Roff_i)**2/R0_i**2)
 * where n_i, R0_i, and Roff_i are parameters.
 */
class MultipleGaussianRadialProfile : public RadialProfile {
   public:
      /** \brief XML constructor
       *
       * Number of Gaussians is given in the element
       * numberOfGaussians
       *
       * The variable names are: n_ii, R0_ii, and Roff_ii one for each Gaussian
       */
      MultipleGaussianRadialProfile ( const utl::Branch &b );

      /** \brief Constructor given variable names
       *
       * The variable names are given in the three vectors and
       * have to be set in the pars instance.
       */ 
      MultipleGaussianRadialProfile ( const utl::Parameters &pars, 
            const std::vector<std::string> &nnames, 
            const std::vector<std::string> &R0names, 
            const std::vector<std::string> &Roffnames );

      /** \brief Constructor given a prefix
       *
       * The variable names are prefix_<x>_<ii> with <x> in [n, R0, Roff] and <ii>
       * a 0 based index for the number of Gaussians given with nGaussians
       */ 
      MultipleGaussianRadialProfile ( const utl::Parameters &pars, const std::string &prefix, const size_t nGaussians );

      virtual double operator () ( double R ) const override;
      virtual void updateVariableValues ( const utl::Variables &vars ) override;

      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

   private:
      std::vector<double> ni, oneOver_R02, Roff;
      std::vector<std::string> fnnames, fR0names, fRoffnames;

};

/** \brief Constant core profile
 *
 * Takes any radial profile and makes it constant core
 *
 * f(R) = 1             if R  < Rcore
 * f(R) = g(R-Rcore)    if R >= Rcore
 *
 * g(R) is the given radial profile
 */
class ConstantCoreRadialProfile : public RadialProfile {
   public:
      /** \brief XML constructor
       *
       * The variable name is Rcore
       * Requires a radial profile.  No need for intent, first is always chosen.
       */
      ConstantCoreRadialProfile ( const utl::Branch &b );

      /** \brief Constructor given variable names
       *
       * The variable names are given with the string.  That
       * has to be set in the pars instance.
       */ 
      ConstantCoreRadialProfile ( std::unique_ptr<RadialProfile> g, const utl::Parameters &pars, const std::string &Rcorename );

      /** \brief Constructor given variable names and values
       *
       * The variable name is given with the string.
       */ 
      ConstantCoreRadialProfile ( std::unique_ptr<RadialProfile> g, double Rcore, const std::string &Rcorename );

      virtual double operator () ( double R ) const override;
      virtual void updateVariableValues ( const utl::Variables &vars ) override;

      void setRcore(double r);
      void setg(std::unique_ptr<RadialProfile> g);

      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

   private:
      std::unique_ptr<RadialProfile> fg;
      double Rcore;
      std::string fRcorename;

};

/** \brief Implements a spline radial profile
 *
 * Uses by default a cubic spline between given points.  Number of points are fixed but everything 
 * else can be a parameter, i.e. both radial location as well as magnitude of points.
 * Needs at least 3 points to be specified.
 *
 * Returns 0 for values that are outside the range of the spline.
 */
class SplineRadialProfile : public RadialProfile {
   public:
      /** \brief XML constructor
       *
       * The number of points and interpolation type are defined in the attributes
       * numberOfPoints and interpolationType, respectively.
       *
       * The variable names are:  radius_ii, and value_ii
       * where ii is the 0 padded 0 based number of the point.  The radii
       * have to be in strictly ascending order.
       *
       * interpolationType is interpreted as a string and can have the values (case insensitive)
       * * LINEAR
       * * POLYNOMIAL
       * * CSPLINE
       * * AKIMA
       */
      SplineRadialProfile ( const utl::Branch &b );

      /** \brief Constructor with a prefix and number of points
       *
       * The variable names are <prefix>_radius_ii and <prefix>_value_ii where
       * ii is an index counting from 0 to numPoints-1.  The radius have to be specified
       * in strictly ascending order with radius_00 being smallest.
       * All variables have to be set in the pars instance.
       */ 
      SplineRadialProfile ( const utl::Parameters &pars, const std::string &prefix, size_t numPoints, const gsl_interp_type *interpType = gsl_interp_cspline );

      /** \brief Constructor with each variable name specified
       *
       * The variable names are given in pairs, radius first, value second.
       * The radius must be given in ascending order, the first index being smallest.
       */
      SplineRadialProfile ( const utl::Parameters &pars, const std::vector< std::pair<std::string,std::string> > &varNames, const gsl_interp_type *interpType = gsl_interp_cspline );

      /** \brief Constructor with each variable name specified and values given
       *
       * The variable names are given in pairs, radius first, value second.
       * The varaible values are also given in such pairs.
       * The radius must be given in ascending order, the first index being smallest.
       */
      SplineRadialProfile ( const std::vector< std::pair<double,double> > &varValues, const std::vector< std::pair<std::string,std::string> > &varNames, const gsl_interp_type *interpType = gsl_interp_cspline );

      virtual ~SplineRadialProfile ();

      virtual double operator () ( double R ) const override;
      virtual void updateVariableValues ( const utl::Variables &vars ) override;

      void setRadii(const std::vector<double> &r);
      void setValues(const std::vector<double> &v);

      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
      virtual size_t getOpenCLNPars() const override;
      virtual std::vector<cl_float> getOpenCLPars() const override;
#endif

   protected:
      double Rmin, Rmax;

      gsl_spline *spline;
      gsl_interp_accel *acc;
      const gsl_interp_type *type;
      
      std::vector< std::pair<std::string,std::string> > fvarNames;

      const std::string XMLtype;

      // So we can re-use it for SplineLog
      void addToDOMinternal( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes,
            const std::string &profileType ) const; 

#ifdef HAVE_OPENCL
      virtual std::string getOpenCLBSearchFun(std::string &name, size_t parIndex) const;
      virtual std::string getOpenCLLinearFun(std::string &name, size_t parIndex) const;
      virtual std::string getOpenCLCSplineFun(std::string &name, size_t parIndex) const;
      virtual std::string getOpenCLAkimaFun(std::string &name, size_t parIndex) const;
#endif
};


/** \brief Implements a spline radial profile, interpolating in log of values
 *
 * Uses by default a cubic spline between logarithm of given points.    
 * Number of points are fixed but everything 
 * else can be a parameter, i.e. both radial location as well as magnitude of points.
 * Of course the values must be positive for logarithm to be valid.
 * Needs at least 3 points to be specified for default interpolation type
 *
 * Returns 0 outside of the interpolation range
 */
class SplineLogRadialProfile : public SplineRadialProfile {
   public:
      /** \brief XML constructor
       *
       * The number of points and interpolation type are defined in the attributes
       * numberOfPoints and interpolationType, respectively.
       *
       * The variable names are:  radius_ii, and value_ii
       * where ii is the 0 padded 0 based number of the point.  The radii
       * have to be in strictly ascending order.
       *
       * interpolationType is interpreted as a string and can have the values (case insensitive)
       * * LINEAR
       * * POLYNOMIAL
       * * CSPLINE
       * * AKIMA
       */
      SplineLogRadialProfile ( const utl::Branch &b );

      /** \brief Constructor with a prefix and number of points
       *
       * The variable names are <prefix>_radius_ii and <prefix>_value_ii where
       * ii is an index counting from 0 to numPoints-1.  The radius have to be specified
       * in strictly ascending order with radius_00 being smallest.  Returns 0 for values that
       * are outside the range of the spline.
       * All variables have to be set in the pars instance.
       */ 
      SplineLogRadialProfile ( const utl::Parameters &pars, const std::string &prefix, size_t numPoints, const gsl_interp_type *interpType = gsl_interp_cspline);

      /** \brief Constructor with each variable name specified
       *
       * The variable names are given in pairs, radius first, value second.
       * The radius must be given in ascending order, the first index being smallest.
       */
      SplineLogRadialProfile ( const utl::Parameters &pars, const std::vector< std::pair<std::string,std::string> > &varNames, const gsl_interp_type *interpType = gsl_interp_cspline );

      /** \brief Constructor with each variable name specified and values given
       *
       * The variable names are given in pairs, radius first, value second.
       * The varaible values are also given in such pairs.
       * The radius must be given in ascending order, the first index being smallest.
       */
      SplineLogRadialProfile ( const std::vector< std::pair<double,double> > &varValues, const std::vector< std::pair<std::string,std::string> > &varNames, const gsl_interp_type *interpType = gsl_interp_cspline );

      virtual void updateVariableValues ( const utl::Variables &vars ) override;
      virtual double operator () ( double R ) const override;

      virtual void addToDOM( xercesc::DOMNode *node, 
            const std::map<std::string,std::string> &attributes ) const override;

#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name, size_t parIndex) const override;
#endif
};



}
#endif
