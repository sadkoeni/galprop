#ifndef Variables_h
#define Variables_h

#include <Parameters.h>
#include <Reader.h>
#include <string>
#include <map>
#include <exception>
#include <iostream>
#include <vector>

namespace utl {

   /** \brief helper function to convert string to number
    */
   template<typename C>
   bool fromString(const std::string &s, C &number) {
      std::istringstream iss(s);
      return !(iss >> number).fail();
   }

/** \brief Handles communication of free variables between models and fitters.
 *
 * Stores the names, value, bounds, and errors of variables.
 */
class Variables {
   public:

      //! The type of the variable
      enum varType {
         FIXED, //!< The variable is fixed, do not include in the fit
         FREE, //!< The variable is free, no boundaries set
         UPPBND, //!< The variable has an upper bound, but not a lower bound
         LOWBND, //!< The variable has a lower bound, but not an upper bound
         BOTHBND //!< The variable has both a lower and an upper bounds
      };

      //! Available priors, NONE is used when adding FIXED variables
      enum PRIOR { NONE, UNIFORM, LOGUNIFORM, GAUSSIAN, CAUCHY, LOGNORMAL };
      
      //! Default constructor to initialize the variables
      Variables() : fmaxLength(13) {} //The length of "Variable Name"
      
      /** \brief Add variables from a utl::Parameters object.
       *
       * \param names is the variable names to look up in the utl::Parameters object
       * \param pars is the utl::Parameters object containing the variables and their value.  
       * The format for the variables in the parameter object has to be
       * \code
       * name = value errorEstimate <lowerBound upperBound>
       * \endcode
       * The limits are optional and to specify only an upperBound make the lowerBound
       * larger than the upperBound
       *
       * For Bayesian analysis it is possible to set priors on the variables with the
       * following syntax in the parameter object
       * \code
       * name_prior = prior par1 par2
       * \endcode
       * prior is the name of the prior and can take the following values
       * * uniform: par1 is the minimum and par2 the maximum
       * * loguniform: par1 is the minimum and par2 the maximum
       * * Cauchy: par1 is the mean and par2 is gamma ( = FWHM/2 )
       * * Gaussian: par1 is the mean and par2 is sigma ( = sqrt(variance) )
       * * lognormal: par1 is the mode and par2 is sigma ( the width parameter )
       *
       * If no prior is given it is assumed from the values and error estimates in 
       * the following way:
       * If no bounds are given we assume a Gaussian prior with value as the mean 
       * and errorEstimate as the sigma.
       * If lower bound is given we assume a LogNormal prior with value as the mode 
       * and errorEstimate as the sigma.  We are basically assuming the lower bound is 0.
       * If both limits are given we assume a uniform prior with lowerBound as minimum
       * and upperBound as maximum.
       * If only upper Bound is given we display an error.
       */
      void add(const std::vector<std::string> & names, const utl::Parameters &pars);

      /** \brief Add variables from an xml branch
       *
       * \param names is the variable ids to look up in the xml branch.  On output it contains the names
       * of the variables added, <name> if it exists otherwise <prefix>_<id>.
       * \param prefix is used as a prefix for the variable names
       * \param b is the xml branch object containing the variables and their value.  
       * The format for the variables in the parameter object has to be
       * \code
       * <variable id="string"> 
       *   <value>double</value> 
       *   <step>double</step>
       *   <min>double</min>
       *   <max>double</max>
       *   <prior>string</prior>
       *   <pval1>double</pval1>
       *   <pval2>double</pval2>
       *   <name>string</name>
       * </variable>
       * \endcode
       * Bounds are only read in if step is given. Single sided bounds are allowed.
       * NOTE: NO OTHER ATTRIBUTE MAY BE IN THE VARIABLE ELEMENT AND NO ATTRIBUTES CAN BE IN THE OTHER ELEMENTS
       *
       * prior is the name of the prior and can take the following values
       * * uniform: pval1 is the minimum and pval2 the maximum
       * * loguniform: pval1 is the minimum and pval2 the maximum
       * * Cauchy: pval1 is the mean and pval2 is gamma ( = FWHM/2 )
       * * Gaussian: pval1 is the mean and pval2 is sigma ( = sqrt(variance) )
       * * lognormal: pval1 is the mode and pval2 is sigma ( the width parameter )
       * For the prior to be set, both pval1 and pval2 must be defined.
       *
       * If no prior is given it is assumed from the values and error estimates in 
       * the following way:
       * If no bounds are given we assume a Gaussian prior with value as the mean 
       * and errorEstimate as the sigma.
       * If lower bound is given we assume a LogNormal prior with value as the mode 
       * and errorEstimate as the sigma.  We are basically assuming the lower bound is 0.
       * If both bounds are given we assume a uniform prior with lowerBound as minimum
       * and upperBound as maximum.
       * If only upper bound is given we display an error.
       *
       * If varName is not given, the variable name is <prefix>_<name>.
       */
      void add( std::vector<std::string> & names, const std::string &prefix, const utl::Branch &b );
      
      /** \brief Add variables from another Variables object 
       *
       * This silently ignores variables already defined in the object by design.
       * This way it is easy to aggregate variables from different models that 
       * may use the same variable.
       */
      void add(const Variables & vars);
      
      /** \brief Add a variable that is fixed
       *
       * \param name is the variable name
       * \param value is its value
       *
       * No prior is used.
       */
      void add(const std::string &name, double value);
      
      /** \brief Add a variable without bounds 
       *
       * \param name is the variable name
       * \param value is its value
       * \param error is its error boundary
       *
       * A Gaussian prior with mean = value and sigma = error is used
       */
      void add(const std::string &name, double value, double error);
      
      /** \brief Add a variable with bounds 
       *
       * \param name is the variable name
       * \param value is its value
       * \param error is its error boundary
       * \param lowerBound is its lower bound
       * \param upperBound is its upper bound
       *
       * A uniform prior is used with min = lowerBound and max = upperBound
       */
      void add(const std::string &name, double value, double error, double lowerBound, double upperBound);

      /** \brief Set the variables from another Variables object
       *
       * Does not add varialbes, only changes the values of current variables.
       * Throws an error if not all of the current variables are defined in the
       * Variables class provided
       */
      void setFromVariables(const Variables &vars);

      /** \brief Return a reference to the variable value */
      double & operator [] (const std::string &name);
      /** \brief Return a constant reference to the variable value */
      const double & operator [] (const std::string &name) const;
      /** \brief Return a reference to the variable value */
      double & value(const std::string &name);
      /** \brief Return a constant reference to the variable value */
      const double & value(const std::string &name) const;

      /** \brief Return the variable error */
      double error(const std::string &name) const;

      /** \brief Return the asymetric errors */
      std::pair<double,double> errors(const std::string &name) const;

      /** \brief Return the upper limit */
      double upperLimit(const std::string &name) const;

      /** \brief Set the error variable.  
       *
       * Note FIXED variables are fixed until the error is set.
       *
       * If the variable was fixed before and the prior is NONE, the prior is set to a GAUSSIAN.
       * Note that the prior is not adjusted in any other case.
       */
      void setError(const std::string &name, double error);

      /** \brief Set asymetric errors
       * No other modification is made
       **/
      void setErrors( const std::string &name, std::pair<double,double> errors);

      /** \brief Set upper limits
       * No other modification is made
       **/
      void setUpperLimit( const std::string &name, double uplim);

      /** \brief Return a vector with all the variable names */
      std::vector<std::string> getNames() const;

      /** \brief Set the upper bound of a variable
       * \param name is the variable name
       * \param bound is the new upper bound.
       *
       * If the upper bound has not been set before and the prior is 
       * still at the default values it is adjusted to uniform if 
       * lower bound is set and none otherwise.
       *
       * Throws an exception if the variable is FIXED
       */
      void setUpperBound(const std::string &name, double bound);
      /** \brief Set the lower bound of a variable
       * \param name is the variable name
       * \param bound is the new lower bound.
       *
       * If the lower bound has not been set before and the prior is 
       * still at the default values it is adjusted to uniform if 
       * lower bound is set and lognormal otherwise.
       *
       * Throws an exception if the variable is FIXED
       */
      void setLowerBound(const std::string &name, double bound);

      /** \brief Remove the upper bound
       *
       * Adjust the priors if they are at default values.  Changes to lognormal
       * if lower bound is set, otherwise it uses a Gaussian.
       */
      void removeUpperBound(const std::string &name);
      /** \brief Remove the lower bound
       *
       * Adjust the priors if they are at default values.  Changes to none
       * if upper bound is set, otherwise it uses a Gaussian.
       */
      void removeLowerBound(const std::string &name);

      /** \brief Set the variables prior
       *
       * Available priors are
       * * UNIFORM: par1 is the minimum and par2 the maximum
       * * LOGUNIFORM: par1 is the minimum and par2 the maximum
       * * CAUCHY: par1 is the mean and par2 is gamma ( = FWHM/2 )
       * * GAUSSIAN: par1 is the mean and par2 is sigma ( = sqrt(variance) )
       * * LOGNORMAL: par1 is the mode and par2 is sigma ( the width parameter )
       */
      void setPrior(const std::string &name, PRIOR prior, double par1, double par2);

      //! Return the prior type
      PRIOR prior(const std::string &name) const;

      //! Return prior parameters
      double priorPar1(const std::string &name) const;
      double priorPar2(const std::string &name) const;

      //! Return the type of the variable
      varType type(const std::string &name) const;

      /** \brief Return a constant reference to the variable upper bound, if it exists */
      double upperBound(const std::string &name) const;
      /** \brief Return a constant reference to the variable lower bound, if it exists */
      double lowerBound(const std::string &name) const;

      //! Fix the variable, set type to FIXED and prior to NONE.  Does not change error or bounds
      void fix(const std::string &name);

      /** \brief Return the values of the variables as a vector.
       *
       * Uses the alphabetical order of the parameters
       */
      std::vector<double> getValueVector() const;

      /** \brief Sets the values of the variables from a vector.
       *
       * Uses the alphabetical order of the parameters and throws an error if the sizes do not match.
       */
      void setValueVector(const std::vector<double> & values);

      //! Clear all set variables
      void clear();

      /** \brief Add the requested variables to the DOMNode
       *
       * This method is designed to be approximately the reverse
       * of the add method using a branch.
       *
       * * \param names is used to identify the variables within the Variables object.
       * * \param ids should be the names used for the id of the variable.
       * 
       * If the variable name is the same as <prefix>_<id> we will not add the 
       * extra <name> element to the branch.  Using the same ids and prefix as
       * in the add command should return the same document structure.
       */
      void addToDOM( xercesc::DOMNode *node, 
            const std::vector<std::string> &names, 
            const std::vector<std::string> &ids, 
            const std::string &prefix ) const;

      /** \brief Default error class */
      class VariableError : public std::exception {
         private:
            std::string eMessage;
         public:
            /** \brief Sets a default error string "Error in Variables class" */
            VariableError() : eMessage("Error in Variables class"){}
            /** \brief Sets a specialized error string */
            VariableError(const std::string & message) : eMessage(message){}
            ~VariableError() throw(){}
            /** \brief Returns the error string */
            virtual const char* what () const throw(){
               return(eMessage.c_str());
            }
      };

      /** \brief Thrown when the variable name does not exist */
      class VariableNotFound : public std::exception {
         private:
            std::string eMessage;
         public:
            /** \brief Sets a default error string "Variable not found" */
            VariableNotFound() : eMessage("Variable not found"){}
            /** \brief Sets a specialized error string */
            VariableNotFound(const std::string & message) : eMessage(message){}
            ~VariableNotFound() throw(){}
            /** \brief Returns the error string */
            virtual const char* what () const throw(){
               return(eMessage.c_str());
            }
      };

      /** \brief Thrown when trying to access a bound that is not there */
      class BoundNotSet : public std::exception {
         private:
            std::string eMessage;
         public:
            /** \brief Sets a default error string "Bound not set" */
            BoundNotSet() : eMessage("Bound not set"){}
            /** \brief Sets a specialized error string */
            BoundNotSet(const std::string & message) : eMessage(message){}
            ~BoundNotSet() throw(){}
            /** \brief Returns the error string */
            virtual const char* what () const throw(){
               return(eMessage.c_str());
            }
      };
      /** \brief Thrown when trying to create a variable with a name that already exists */
      class VariableAlreadyCreated : public std::exception {
         private:
            std::string eMessage;
         public:
            /** \brief Sets a default error string "Variable already created" */
            VariableAlreadyCreated() : eMessage("Variable already created"){}
            /** \brief Sets a specialized error string */
            VariableAlreadyCreated(const std::string & message) : eMessage(message){}
            ~VariableAlreadyCreated() throw(){}
            /** \brief Returns the error string */
            virtual const char* what () const throw(){
               return(eMessage.c_str());
            }
      };
      /** \brief Output operator for the Variables class */
      friend std::ostream & operator << (std::ostream &os, const Variables &vars);

      //! For easy convertion from PRIOR to string
      static std::string priorToString(PRIOR prior);
      //! For easy convertion from string to PRIOR
      static PRIOR stringToPrior(const std::string & priorString);

   private:
      /** \brief A structure to store each variable */
      struct variable {
         double value, error, uBound, lBound;
         double uError, lError, upLim;
         bool uSet, lSet;
         varType type;
         PRIOR prior;
         double priorPar1, priorPar2;
      };

      std::map<std::string, variable> fvars;
      typedef std::map<std::string, variable>::iterator varIt;
      typedef std::map<std::string, variable>::const_iterator cvarIt;

      int fmaxLength; //!< The maximum lenght of a variable name

      inline cvarIt getIterator(const std::string &name) const;
      inline varIt getIterator(const std::string &name);
};

}
#endif
