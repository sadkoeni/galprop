#include <Nuclei_Interface.h>

using namespace nuclei;
using namespace std;

#include <config.h>

#define SET_SIGMA_F77 F77_FUNC_(set_sigma,SET_SIGMA)
#ifdef __cplusplus
extern "C" 
#endif
void SET_SIGMA_F77(int*, char*, int*);//set_sigma_(int*); // IMOS20020502

#define YIELDX_F77 F77_FUNC(yieldx,YIELDX)
#ifdef __cplusplus
extern "C" 
#endif
void YIELDX_F77(int*, int*, int*, int*, float*, float*); // TS code IMOS20020502

#define WSIGMA_F77 F77_FUNC(wsigma,WSIGMA)
#ifdef __cplusplus
extern "C" 
#endif
double WSIGMA_F77(int*, int*, int*, int*, double*);// Webber's code IMOS20020502

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// Webber's isotopic production cross section  IMOS20020502
double nuclei::wsigma_cc(int IZ, int IA, int IZF, int IAF, double E) {

  return WSIGMA_F77(&IZ,&IA,&IZF,&IAF,&E);

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// Silberberg & Tsao isotopic production cross section  IMOS20020502
double nuclei::yieldx_cc(int IZ, int IA, int IZF, int IAF, double energy) {

  float CSmb, E = float(energy);
  YIELDX_F77(&IZ,&IA,&IZF,&IAF,&E,&CSmb);
  return( 1.*CSmb );

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// initialization of Webber's code
void nuclei::set_sigma_cc(const std::string& path) { 
   
  int cdr = 99;

  const std::string fullFilename = path + "/" + "WNEWTR_082693.CDR.dat";

  char* fn = const_cast<char*>(fullFilename.c_str());
  printf(fn);

  int length = fullFilename.size();

  SET_SIGMA_F77(&cdr, fn, &length); // IMOS20020502

}
