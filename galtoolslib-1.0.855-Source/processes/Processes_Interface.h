#ifndef _processes_interface_h_
#define _processes_interface_h_

#include <string>
#include <valarray>

namespace processes {

  void aprtab_cc(const std::string& path);
  double apspec_cc(double e0, double epbar, int iap, int iat);

  double pp_meson_cc(double Esec, double Pp1, int NA1, int NA2, int key1);
  double sighad_cc(int IS, double PA, double PZ, double TA, double TZ, double E);
  void sigtap_cc(int ISS, const std::string& path);
  double e_loss_compton_cc(double w, double gam);
  double bremss_spec_cc(double Egam, double E0, int IZ1, int Ne1);
  double synchrotron_cc(double gamma, double nu, double B);
  double antiproton_cc(int key, double Pap1, double Pp1, int NZ1, int NA1, int NZ2, int NAZ);
  void nucleon_cs(int option, double Ek, int Zp, int Zt, int At,
		  double *PP_inel, double *PA_inel, 
		  double *aPP_non, double *aPA_non,
		  double *aPP_ann, double *aPA_ann);
  double nucleon_loss(int z, int a, double emev, double nhcm3, 
		      double nhicm3, double he_to_h,
		      double* aion, double* coul);
  double electron_loss(double emev, 
		       double nhcm3, double nhicm3, double he_to_h, 
		       double uevcm3, double bevcm3, 
		       double* aion, double* coul, 
		       double* brem1, double* brem2, 
		       double* sync, double* cmptn);

  std::valarray<double> kamae_gamma_param_nd(double);
  std::valarray<double> kamae_gamma_param_diff(double);
  std::valarray<double> kamae_gamma_param_delta(double);
  std::valarray<double> kamae_gamma_param_res(double);
  std::valarray<double> kamae_elec_param_nd(double);
  std::valarray<double> kamae_elec_param_diff(double);
  void kamae_elec_param_delta(double, double*);
  std::valarray<double> kamae_elec_param_res(double);
  std::valarray<double> kamae_posi_param_nd(double);
  std::valarray<double> kamae_posi_param_diff(double);
  std::valarray<double> kamae_posi_param_delta(double);
  std::valarray<double> kamae_posi_param_res(double);
  // functions for component differential cross section
  double kamae_nd(double, double, std::valarray<double>, int);
  double kamae_diff(double, double, std::valarray<double>);
  double kamae_delta(double, double, std::valarray<double>);
  double kamae_res(double, double, std::valarray<double>);
  // functions for total inelastic differential cross sections
  double 
    kamae(double Esec, double Pp2, int NA11, int NA21, int key, 
	  std::valarray<double> params0, std::valarray<double> params1, 
	  std::valarray<double> params2, std::valarray<double> params3);

  double blattnig_gamma(double Esec, double Pp2, int NA11, int NA21, int key);

}

#endif
