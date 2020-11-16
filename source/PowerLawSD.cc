#include "PowerLawSD.h"

#include <ErrorLogger.h>

#include <algorithm>

PowerLawSD::PowerLawSD(const utl::Parameters &pars):
   SpectralDistribution(pars)
{
   //Read the type and set it
   std::string type;
   pars.getParameter("powerLawType", type);
   setType(type);
}

bool PowerLawSD::validateParameters(const std::vector<double> *spPars) const
{
   //The size must be odd
   if (spPars->size()%2 == 0) {
      WARNING("Even numbers of parameters");
      return false;
   }

   //The breaks must be ascending
   for (size_t i(3); i < spPars->size(); i+=2)
      if ((*spPars)[i] <= (*spPars)[i-2]) {
         WARNING("Breaks aren't ascending");
         return false;
      }

   //Passed all tests
   return true;

}

double PowerLawSD::value(const std::vector<double> *spPars, const Particle &particle,
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

   //Loop through the breaks and find the index
   double ind = (*spPars)[0];
   double norm(1);
   for (size_t i(1); i < spPars->size(); i+=2)
      if (x >= (*spPars)[i]) {
         norm *= pow( (*spPars)[i], (*spPars)[i+1]-ind );
         ind = (*spPars)[i+1];
      }

   return norm * pow(x, -ind);
}

void PowerLawSD::setType(std::string name)
{
   std::transform(name.begin(), name.end(), name.begin(), (int(*)(int)) std::toupper);

   if ( name == "EKIN" )
      plType = EKIN;
   else if (name == "ETOT")
      plType = ETOT;
   else if (name == "PTOT")
      plType = PTOT;
   else if (name == "PNUC")
      plType = PNUC;
   else if (name == "RIGIDITY")
      plType = RIGIDITY;
   else if (name == "GAMMA")
      plType = GAMMA;
   else if (name == "BETA")
      plType = BETA;
   else {
      std::ostringstream buf;
      buf<<"No such variable \""<<name<<"\".";
      ERROR(buf.str());
      throw(std::runtime_error("Fix your configuration"));
   }
}
