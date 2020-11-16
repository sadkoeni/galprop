#include "SpectralDistribution.h"

#include <ErrorLogger.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

struct integrationParameters {
   integrationParameters(const std::vector<double> *p, const Particle &pa, const SpectralDistribution::Type t, const SpectralDistribution *s):
      pars(p), particle(pa), type(t), sd(s) {}
   const std::vector<double> *pars;
   const Particle &particle;
   const SpectralDistribution::Type type;
   const SpectralDistribution* sd;
};

static double integrationFunction(double x, void *params)
{
   auto p = static_cast<integrationParameters*>(params);
   switch (p->type) {
      case SpectralDistribution::ETOT:
         return x*x*p->sd->valueEtot(p->pars, p->particle, x);
      case SpectralDistribution::EKIN:
         return x*x*p->sd->valueEkin(p->pars, p->particle, x);
      case SpectralDistribution::PTOT:
         return x*x*p->sd->valuepTot(p->pars, p->particle, x);
      case SpectralDistribution::PNUC:
         return x*x*p->sd->valuepNuc(p->pars, p->particle, x);
      case SpectralDistribution::RIGIDITY:
         return x*x*p->sd->valueRigidity(p->pars, p->particle, x);
      case SpectralDistribution::GAMMA:
         return x*x*p->sd->valueGamma(p->pars, p->particle, x);
      case SpectralDistribution::BETA:
         return p->sd->valueBeta(p->pars, p->particle, x);
      default:
         ERROR("Need to fix the code, this should never happen");
         throw(std::logic_error("This should be fixed"));
   }
}


std::vector<double> SpectralDistribution::evaluate(const std::vector<double> *spPars, const Particle &particle) const
{
   std::vector<double> output(particle.n_pgrid);
   //Use the total energy to convert to others and call the value operator
   double Ekin, pTot, pNuc, Rigidity, Gamma, Beta;

   for (int i(0); i < particle.n_pgrid; ++i) {
      const double Etot = particle.Etot[i];
      EtotToOther(particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
      output[i] = value(spPars, particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
   }

   return output;
}
      
double SpectralDistribution::valueEkin(const std::vector<double> *spPars, const Particle &particle, const double Ekin) const
{
   double Etot, pTot, pNuc, Rigidity, Gamma, Beta;
   EkinToOther(particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
   return value(spPars, particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
}

double SpectralDistribution::valueEtot(const std::vector<double> *spPars, const Particle &particle, const double Etot) const
{
   double Ekin, pTot, pNuc, Rigidity, Gamma, Beta;
   EtotToOther(particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
   return value(spPars, particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
}

double SpectralDistribution::valuepTot(const std::vector<double> *spPars, const Particle &particle, const double pTot) const
{
   double Ekin, Etot, pNuc, Rigidity, Gamma, Beta;
   pTotToOther(particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
   return value(spPars, particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
}

double SpectralDistribution::valuepNuc(const std::vector<double> *spPars, const Particle &particle, const double pNuc) const
{
   double Ekin, Etot, pTot, Rigidity, Gamma, Beta;
   pNucToOther(particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
   return value(spPars, particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
}

double SpectralDistribution::valueRigidity(const std::vector<double> *spPars, const Particle &particle, const double Rigidity) const
{
   double Ekin, Etot, pTot, pNuc, Gamma, Beta;
   RigidityToOther(particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
   return value(spPars, particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
}

double SpectralDistribution::valueGamma(const std::vector<double> *spPars, const Particle &particle, const double Gamma) const
{
   double Ekin, Etot, pTot, pNuc, Rigidity, Beta;
   GammaToOther(particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
   return value(spPars, particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
}

double SpectralDistribution::valueBeta(const std::vector<double> *spPars, const Particle &particle, const double Beta) const
{
   double Ekin, Etot, pTot, pNuc, Rigidity, Gamma;
   BetaToOther(particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
   return value(spPars, particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
}


double SpectralDistribution::integralEkin(const std::vector<double> *spPars, const Particle &particle, const double EkinMin, const double EkinMax) const
{
   return integral(spPars, particle, EkinMin, EkinMax, EKIN);
}

double SpectralDistribution::integralEtot(const std::vector<double> *spPars, const Particle &particle, const double EtotMin, const double EtotMax) const
{
   return integral(spPars, particle, EtotMin, EtotMax, ETOT);
}

double SpectralDistribution::integralpTot(const std::vector<double> *spPars, const Particle &particle, const double pTotMin, const double pTotMax) const
{
   return integral(spPars, particle, pTotMin, pTotMax, PTOT);
}

double SpectralDistribution::integralpNuc(const std::vector<double> *spPars, const Particle &particle, const double pNucMin, const double pNucMax) const
{
   return integral(spPars, particle, pNucMin, pNucMax, PNUC);
}

double SpectralDistribution::integralRigidity(const std::vector<double> *spPars, const Particle &particle, const double RigidityMin, const double RigidityMax) const
{
   return integral(spPars, particle, RigidityMin, RigidityMax, RIGIDITY);
}

double SpectralDistribution::integralGamma(const std::vector<double> *spPars, const Particle &particle, const double GammaMin, const double GammaMax) const
{
   return integral(spPars, particle, GammaMin, GammaMax, GAMMA);
}

double SpectralDistribution::integralBeta(const std::vector<double> *spPars, const Particle &particle, const double BetaMin, const double BetaMax) const
{
   return integral(spPars, particle, BetaMin, BetaMax, BETA);
}

double SpectralDistribution::integral(
      const std::vector<double> *spPars, const Particle &particle,
      const double min, const double max, Type type) const
{
   //Defaults to using integral over total energy
   const size_t limit(10000);
   auto ws = gsl_integration_workspace_alloc(limit);

   double results, error;

   integrationParameters p(spPars, particle, type, this);

   gsl_function F;
   F.function = &integrationFunction;
   F.params = &p;

   const int errorKey = gsl_integration_qag(&F, min, max, 0, 1e-4, limit, GSL_INTEG_GAUSS41, ws, &results, &error);

   if (errorKey != 0) {
      std::string error = gsl_strerror(errorKey);
      ERROR(error);
      throw(std::runtime_error("GSL integration failed"));
   }

   gsl_integration_workspace_free(ws);
   
   return results;
}


void SpectralDistribution::EkinToOther(const Particle &particle, const double Ekin, double & Etot, double & pTot, double & pNuc, double & Rigidity, double & Gamma, double & Beta)
{
   const double A = particle.A == 0 ? 1 : particle.A;
   Etot = Ekin*A + particle.mass;
   double EkinTmp;
   EtotToOther(particle, EkinTmp, Etot, pTot, pNuc, Rigidity, Gamma, Beta);
}

void SpectralDistribution::EtotToOther(const Particle &particle, double &Ekin, const double Etot, double & pTot, double & pNuc, double & Rigidity, double & Gamma, double & Beta)
{
   const double A = particle.A == 0 ? 1 : particle.A;
   Ekin = (Etot-particle.mass)/A;
   pTot = sqrt(Etot*Etot-particle.mass*particle.mass);
   pNuc = pTot/A;
   Gamma = Etot/particle.mass;
   Beta = sqrt((1.-1./Gamma)*(1.+1./Gamma));
   Rigidity = pTot/abs(particle.Z);
}

void SpectralDistribution::pTotToOther(const Particle &particle, double &Ekin, double & Etot, const double pTot, double & pNuc, double & Rigidity, double & Gamma, double & Beta)
{
   Etot = sqrt(pTot*pTot + particle.mass*particle.mass);
   double pTotTmp;
   EtotToOther(particle, Ekin, Etot, pTotTmp, pNuc, Rigidity, Gamma, Beta);
}

void SpectralDistribution::pNucToOther(const Particle &particle, double &Ekin, double & Etot, double & pTot, const double pNuc, double & Rigidity, double & Gamma, double & Beta)
{
   const double A = particle.A == 0 ? 1 : particle.A;
   Etot = sqrt(pNuc*pNuc*A*A + particle.mass*particle.mass);
   double pNucTmp;
   EtotToOther(particle, Ekin, Etot, pTot, pNucTmp, Rigidity, Gamma, Beta);
}

void SpectralDistribution::RigidityToOther(const Particle &particle, double &Ekin, double & Etot, double & pTot, double & pNuc, const double Rigidity, double & Gamma, double & Beta)
{
   Etot = sqrt(Rigidity*Rigidity*particle.Z*particle.Z + particle.mass*particle.mass);
   double RigidityTmp;
   EtotToOther(particle, Ekin, Etot, pTot, pNuc, RigidityTmp, Gamma, Beta);
}

void SpectralDistribution::GammaToOther(const Particle &particle, double &Ekin, double & Etot, double & pTot, double & pNuc, double & Rigidity, const double Gamma, double & Beta)
{
   Etot = Gamma*particle.mass;
   double GammaTmp;
   EtotToOther(particle, Ekin, Etot, pTot, pNuc, Rigidity, GammaTmp, Beta);
}

void SpectralDistribution::BetaToOther(const Particle &particle, double &Ekin, double & Etot, double & pTot, double & pNuc, double & Rigidity, double & Gamma, const double Beta) 
{
   Etot = particle.mass/sqrt( (1-Beta)*(1+Beta) );
   double BetaTmp;
   EtotToOther(particle, Ekin, Etot, pTot, pNuc, Rigidity, Gamma, BetaTmp);
}

