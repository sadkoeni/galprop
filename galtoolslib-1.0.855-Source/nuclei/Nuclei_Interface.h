#ifndef _nuclei_interface_h_
#define _nuclei_interface_h_

#include <string>

namespace nuclei {

  void set_sigma_cc(const std::string& path);
  double wsigma_cc(int Z, int A, int ZF, int AF, double energy);
  double yieldx_cc(int Z, int A, int ZF, int AF, double energy);
  int He_to_H_CS(double E1, int IZI, int IAI, int IZF, int IAF, double* CSratio,double* CStot_ratio);
  double sigma_boron_dec_heinbach_simon(int IZ, int IA, int JZ, int JA, double EJ);
  void Kcapture_cs(double Ek, int Zp, int Zt, double* attach, double *strip);

}

#endif
