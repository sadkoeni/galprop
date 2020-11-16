#include "SmoothPowerLawSD.h"

#include <ErrorLogger.h>

#include <algorithm>

SmoothPowerLawSD::SmoothPowerLawSD(const utl::Parameters &pars):
   PowerLawSD(pars)
{
}

bool SmoothPowerLawSD::validateParameters(const std::vector<double> *spPars) const
{
   //The size must be a multiple of 3 plus 1
   if (spPars->size()%3 != 1)
      return false;

   //Smoothing must be larger than zero for all breaks
   for (size_t i(2); i < spPars->size(); i+=3)
      if ((*spPars)[i] <= 0)
         return false;

   //The breaks must be ascending
   for (size_t i(4); i < spPars->size(); i+=3)
      if ((*spPars)[i] <= (*spPars)[i-3])
         return false;

   //Passed all tests
   return true;

}

double SmoothPowerLawSD::value(const std::vector<double> *spPars, const Particle &particle,
      const double Ekin, const double Etot, const double pTot, const double pNuc, 
      const double Rigidity, const double Gamma, const double Beta) const
{

   //Select the correct value to use
   double x(0);
   switch (plType) {
      case ETOT:
         x = Etot;
         break;
      case EKIN:
         x = Ekin;
         break;
      case PTOT:
         x = pTot;
         break;
      case PNUC:
         x = pNuc;
         break;
      case RIGIDITY:
         x = Rigidity;
         break;
      case GAMMA:
         x = Gamma;
         break;
      case BETA:
         x = Beta;
         break;
      default:
         ERROR("Need to fix the code, this should never happen");
         throw(std::logic_error("This should be fixed"));
   }

   //Loop through the breaks and create the function
   //Calculate the log of the function, it is more stable far away from the breaks
   double output =  -(*spPars)[0]*log(x);
   for (size_t i(1); i < spPars->size(); i+=3) {
      const double dind = (*spPars)[i-1]-(*spPars)[i+2];
      const double xlog = log(x/(*spPars)[i]);
      //This should be more than accurate enough for doubles
      if (fabs(dind)/(*spPars)[i+1]*xlog > 100) {
         output += dind*xlog;
      } else {
         output += std::copysign((*spPars)[i+1], dind) * log(1 + exp(fabs(dind)/(*spPars)[i+1]*xlog));
      }
   }

   return exp(output);
}

