#ifndef SPECTRAL_DISTRIBUTION_H
#define SPECTRAL_DISTRIBUTION_H

#include "Particle.h"
#include <Parameters.h>

#include <memory>

/** \brief Interface class for the spectral distribution part of the uniform spectrum source class
 *
 * The purpose of this class is to set the injection spectrum of
 * each particle as a function of energy/momentum/rigidity/beta/gamma or
 * some combination thereof.  It provides interfaces for integration
 * over specific ranges of those as well as direct evaluation given 
 * a specific number of one of those.  The base class provides 
 * convenience functions to calculate all from any single one.
 *
 * The class is set up using a parameters object and all public functions 
 * require a particle instance and a double valued vector containing spectral 
 * parameters.
 *
 * The base class by default uses the convenience function to transfer between 
 * the different inputs and calls a protected global value function.  It also
 * performs automatic numerical integration using those functions so only
 * one function needs to be written.  There are, however, possibilities to 
 * override those if analytic integration is possible.  In practice this means you
 * should override the functions most closely related to the variable used to 
 * calculate the spectrum.
 *
 * The classes should register using capabilites from Registry.h in SourceClass_UniformSpectra.cc.
 * There should be lines there you can copy and modify.
 */
class SpectralDistribution 
{
   public:
      //! Construct from a parameters object.  See derived classes for available parameters.
      SpectralDistribution(const utl::Parameters &pars) {}
      virtual ~SpectralDistribution() {}

      /** Return the spectrum at points evaluated at particle.p[ip]
       *
       * This does not need special normalization as the caller is
       * expected to normalize the spectrum.
       *
       * Uses by default the value function
       */
      virtual std::vector<double> evaluate(const std::vector<double> *spPars, const Particle &particle) const;

      /** Validate the spectral parameters
       *
       * Returns false if the parameters are not valid
       */
      virtual bool validateParameters(const std::vector<double> *spPars) const = 0;

      //! Value at a given kinetic energy per nucleon
      virtual double valueEkin(const std::vector<double> *spPars, const Particle &particle, const double Ekin) const;
      //! Value at a given total energy
      virtual double valueEtot(const std::vector<double> *spPars, const Particle &particle, const double Etot) const;
      //! Value at a given total momentum
      virtual double valuepTot(const std::vector<double> *spPars, const Particle &particle, const double pTot) const;
      //! Value at a given momentum per nucleon
      virtual double valuepNuc(const std::vector<double> *spPars, const Particle &particle, const double pNuc) const;
      //! Value at a given rigidity
      virtual double valueRigidity(const std::vector<double> *spPars, const Particle &particle, const double Rigidity) const;
      //! Value at a given Lorentz factor
      virtual double valueGamma(const std::vector<double> *spPars, const Particle &particle, const double Gamma) const;
      //! Value at a given velocity (normalized to c)
      virtual double valueBeta(const std::vector<double> *spPars, const Particle &particle, const double Beta) const;
   
      //! Integral over a range of kinetic energy per nucleon
      virtual double integralEkin(const std::vector<double> *spPars, const Particle &particle, const double EkinMin, const double EkinMax) const;
      //! Integral over a range of total energy
      virtual double integralEtot(const std::vector<double> *spPars, const Particle &particle, const double EtotMin, const double EtotMax) const;
      //! Integral over a range of total momentum
      virtual double integralpTot(const std::vector<double> *spPars, const Particle &particle, const double pTotMin, const double pTotMax) const;
      //! Integral over a range of momentum per nucleon
      virtual double integralpNuc(const std::vector<double> *spPars, const Particle &particle, const double pNucMin, const double pNucMax) const;
      //! Integral over a range of rigidity
      virtual double integralRigidity(const std::vector<double> *spPars, const Particle &particle, const double RigidityMin, const double RigidityMax) const;
      //! Integral over a range of Lorentz factor
      virtual double integralGamma(const std::vector<double> *spPars, const Particle &particle, const double GammaMin, const double GammaMax) const;
      //! Integral over a range of velocity (normalized to c)
      virtual double integralBeta(const std::vector<double> *spPars, const Particle &particle, const double BetaMin, const double BetaMax) const;
   
      //icpc 13.0.1 does not support scoped enums in switch
      enum Type { ETOT, EKIN, PTOT, PNUC, RIGIDITY, GAMMA, BETA};

   protected:
      //Returns a value at a certain point.  All the properties have been precalculated and are inputs
      virtual double value(const std::vector<double> *spPars, const Particle &particle,
            const double Ekin, const double Etot, const double pTot, const double pNuc, const double Rigidity, const double Gamma, const double Beta) const = 0;

      //Convert from one to other
      static void EkinToOther(const Particle &particle, const double Ekin, double & Etot, 
            double & pTot, double & pNuc, double & Rigidity, double & Gamma, double & Beta);
      static void EtotToOther(const Particle &particle, double & Ekin, const double Etot, 
            double & pTot, double & pNuc, double & Rigidity, double & Gamma, double & Beta);
      static void pTotToOther(const Particle &particle, double & Ekin, double & Etot, 
            const double pTot, double & pNuc, double & Rigidity, double & Gamma, double & Beta);
      static void pNucToOther(const Particle &particle, double & Ekin, double & Etot, 
            double & pTot, const double pNuc, double & Rigidity, double & Gamma, double & Beta);
      static void RigidityToOther(const Particle &particle, double & Ekin, double & Etot, 
            double & pTot, double & pNuc, const double Rigidity, double & Gamma, double & Beta);
      static void GammaToOther(const Particle &particle, double & Ekin, double & Etot, 
            double & pTot, double & pNuc, double & Rigidity, const double Gamma, double & Beta);
      static void BetaToOther(const Particle &particle, double & Ekin, double & Etot, 
            double & pTot, double & pNuc, double & Rigidity, double & Gamma, const double Beta);
      
      //Performs numerical integration using the value function from min to max on type.
      double integral( const std::vector<double> *spPars, const Particle &particle, const double min, const double max, Type type) const;

};

#endif
