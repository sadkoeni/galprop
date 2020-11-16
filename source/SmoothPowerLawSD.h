#ifndef SMOOTH_POWER_LAW_ST_H
#define SMOOTH_POWER_LAW_ST_H

#include "PowerLawSD.h"

class SmoothPowerLawSD : public PowerLawSD
{
   public:
      /** Construct a new smoothly joined broken power law spectral distribution
       *
       * The class supports arbitrary number of breaks and the variable 
       * used in the power law can be set with
       * * powerLawType = ETOT | EKIN | PTOT | PNUC | RIGIDITY | GAMMA | BETA
       * with the key defined as
       * * ETOT : total energy
       * * EKIN : kinetic energy per nucleon [ (ETOT-m0)/A ]
       * * PTOT : total momentum
       * * PNUC : momentum per nucleon [ PTOT/A ]
       * * RIGIDITY : rigidity [ PTOT/|Z| ]
       * * GAMMA : Lorentz factor [ ETOT/m0 ]
       * * BETA : speed in units of c
       *
       * Assumes the double vector spectral parameter has the format
       * <index> [<break>, <smoothing>, <index>]*
       * where the first index is mandatory and the breaks must be given
       * in ascending order (smallest first)
       *
       * Smaller smoothing parameter gives a sharper break, a value of 1 means
       * break is smoothed over roufhly a decade.
       *
       * This class follows the convention of negative indices, i.e.
       * F(x) = x**(-index)
       */
      SmoothPowerLawSD(const utl::Parameters &pars);

      /** Validate the spectral parameters
       *
       * Make sure there is even number of pairs and the breaks
       * are ascending.
       */
      virtual bool validateParameters(const std::vector<double> *spPars) const ;

   protected:
      //Need to define the value function
      virtual double value(const std::vector<double> *spPars, const Particle &particle,
            const double Ekin, const double Etot, const double pTot, const double pNuc, 
            const double Rigidity, const double Gamma, const double Beta) const ;

};

#endif
