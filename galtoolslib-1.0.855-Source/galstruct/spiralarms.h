#ifndef SPIRALARMS_H
#define SPIRALARMS_H

#include "radialprofiles.h" //For sech2

#include <cmath>

namespace GalacticStructure {

/** \brief Distribution around spiral arm structure
 *
 * Pure virtual class for a scale function.  It is just
 * a function taking a double and returning a double.
 *
 * In general it should decline with |x| and be 1 at 0.
 *
 * Can also be used for disk scaling, etc.
 */
class ScaleFunction {
   public:
      virtual ~ScaleFunction() {}
      virtual double operator() (double d) const = 0;

      //! The name used in the registry
      virtual std::string name() const = 0;
#ifdef HAVE_OPENCL
      /** \brief Provides an OpenCL code for the function
       *
       * The function prototype is
       * float <name> ( float x )
       * The function name is combined of a class specific prefix and a user
       * specified postfix (see below).
       *
       * \param name is on input the user specified postfix but contains on output
       * the entier function name.
       * \return the code in a string
       */
      virtual std::string getOpenCLFunction(std::string &name) const = 0;
   protected:
      //! Convenience function to create the code with a simple return command and x as the variable
      std::string createFunction(const std::string &name, const std::string &returnCommand) const;
#endif
};

/** \brief Distribution around spiral arm structure
 *
 * exp(-d)
 */
class ExpFunction : public ScaleFunction {
   public:
      virtual double operator() (double d) const override { return exp(-d); }
      virtual std::string name() const override { return "Exponential"; }
#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name) const override;
#endif
};

/** \brief Distribution around spiral arm structure
 *
 * exp(-d**2)
 */
class GaussianFunction : public ScaleFunction {
   public:
      virtual double operator() (double d) const override { return exp(-d*d); }
      virtual std::string name() const override { return "Gaussian"; }
#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name) const override;
#endif
};

/** \brief Distribution around spiral arm structure
 *
 * sech(d)**2
 */
class Sech2Function : public ScaleFunction {
   public:
      virtual double operator() (double d) const override { return sech2(d); }
      virtual std::string name() const override { return "Sech2"; }
#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name) const override;
#endif
};

/** \brief Distribution around spiral arm structure
 *
 * sech(d)
 */
class SechFunction : public ScaleFunction {
   public:
      virtual double operator() (double d) const override { return sech(d); }
      virtual std::string name() const override { return "Sech"; }
#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name) const override;
#endif
};

/** \brief Distribution around spiral arm structure
 *
 * sqrt(sech(d))
 */
class Sech1_2Function : public ScaleFunction {
   public:
      virtual double operator() (double d) const override { return sqrt(sech(d)); }
      virtual std::string name() const override { return "Sech1_2"; }
#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name) const override;
#endif
};

/** \brief Distribution around spiral arm structure
 *
 * H(1-d) 
 */
class StepFunction : public ScaleFunction {
   public:
      virtual double operator() (double d) const override { return d <= 1 ? 1 : 0; }
      virtual std::string name() const override { return "Step"; }
#ifdef HAVE_OPENCL
      virtual std::string getOpenCLFunction(std::string &name) const override;
#endif
};

/** \brief Handle a single spiral arm
 * 
 * Defaults to a logarithmic arm but can be extended to handle more
 * general arms.   
 *
 * ArmScaleFunction gives the behavior as a function of d/w, where
 * d is the distance from the arm and w is the width.  
 * */
class ArmFunction {
  
public:
  ArmFunction( double a,
	       double rMin,
	       double phiMin,
	       double phiExtent,
	       double width,
               std::unique_ptr<ScaleFunction> scale);
  
  virtual double operator () ( double radius, double phi ) const;
  
  double a() const;
  double rMin() const;
  double phiMin() const;
  double phiExtent() const;
  double rMax() const;
  double width() const;
  
  void setA( double a );
  void setRMin( double rMin );
  void setPhiMin( double phiMin );
  void setPhiExtent( double phiExtent );
  void setPhiExtentFromRMax( double rMax );
  void setWidth( double width );
  
  const std::vector< std::pair<double, double> > & RPArmData() const;
  
  //! The theta bins must be linearly distributed.
  void setRPArmData( const std::vector< std::pair<double, double> > &RPArmData );

  //! Test function for comparing armdist methods.  NOT PART OF A STABLE API
  void CompareArmDist( double r, double phi, double &accurate, double &approx ) const;

  //! Return the name of the scale function
  std::string scaleName() { return fscale->name(); }
  
#ifdef HAVE_OPENCL
  /** \brief Provides an OpenCL code for the function
   *
   * WARNING.  This function only implements the logarithmic spiral.
   * The NE2001 model will be imperfect in OpenCL.
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
   * This number should not change dynamically and should be set on class 
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
  double fA, fRMin, fPhiMin, fPhiExtent, fWidth, fOneOverWidth, fCosAlpha;
  std::unique_ptr<ScaleFunction> fscale;

  double (GalacticStructure::ArmFunction::* armDistPointer)(double radius, double phi) const;

  //Brute force method for comparison
  double armDistBF(double radius, double phi) const;
  //Does work and is reasonably fast
  double armDist(double radius, double phi) const;
  void findBestPoint(size_t &imin, double &sqrsmin, size_t ibegin, double rbegin, double rad, double x, double y) const;

  //New method using linear approximation and tangents.
  double armDistLinear(double radius, double phi) const;

  //Method using distance on a circle with same radius
  double armDistRadial(double radius, double phi) const;

  //Old semi-fast method, more reliable
  double armDistOld(double radius, double phi) const;

  double angle(double radius) const;
  double radius(double angle) const;
  
  size_t nStepsPerRing;
  
  void calculateXYdata();
  void calculateRPdata();
  
  std::vector< std::pair<double, double> > fXYArmData, fRPArmData;
  double fphiStep;
  
};

}
#endif
