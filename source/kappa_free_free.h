#ifndef KAPPA_FREE_FREE_h
#define KAPPA_FREE_FREE_h

#include "los_integration.h"
#include <valarray>

namespace synchro {

class KappaFreeFree : public SM::LOSfunction < std::valarray<double> > {
  public:
    KappaFreeFree(const std::valarray<double> &nu_synch, const double clumping, const double T) : nu(nu_synch), clumping_factor(clumping), Te(T) {}
    virtual std::valarray<double> operator () ( const double x, const double y, const double z, const vec3 &dir ) const;
  private:
    const std::valarray<double> &nu;
    const double clumping_factor;
    const double Te;
};

class EmissFreeFree : public SM::LOSfunction<std::valarray<double> > {
  public:
    EmissFreeFree(const std::valarray<double> &nu_synch, const double clumping, const double T) : nu(nu_synch), clumping_factor(clumping), Te(T) {}
    virtual std::valarray<double> operator () ( const double x, const double y, const double z, const vec3 &dir ) const;
  private:
    const std::valarray<double> &nu;
    const double clumping_factor;
    const double Te;
};

inline double gff(const double nu, const double Te);
inline double kff(const double nu, const double Ne, const double Te);
inline double eff(const double nu, const double Ne, const double Te);

double kappa_free_free(double nu,double Ne, double Te,  double &emiss_free_free, int options, int debug); //AWS20110704
double kappa_free_free_2D(double nu, double R,           double z,  double Te, double clumping_factor, double &emiss_free_free, int options, int debug); //AWS20110704
double kappa_free_free_3D(double nu, double x, double y, double z,  double Te, double clumping_factor, double &emiss_free_free, int options, int debug); //AWS20110704

}

#endif
