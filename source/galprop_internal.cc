#include <ErrorLogger.h>

#include <galprop_classes.h>
#include <galprop_internal.h>

//#include <config.h>

#include <constants.h>

extern Galprop* gGalprop;

//IMOS20060420 transferred from fort_interface1.cc
// the routine is called by FORTRAN routine emiss(r) (in file cfactor.f)
// hence underscore is appended and extern "C" supplied

#define ISRF_ENERGY_DENSITY_F77 F77_FUNC_(isrf_energy_density,ISRF_ENERGY_DENSITY)
#ifdef __cplusplus
extern "C" 
#endif
void ISRF_ENERGY_DENSITY_F77(float* r, float* z, float* energy_density) {

  *energy_density = gGalprop->isrf_energy_density(*r,*z);
//   cout<<"energy_density = "<<*energy_density<<endl;

}

//#define RHO_DARKSUSY_F77 F77_FUNC_(rho_darksusy,RHO_DARKSUSY)
//#ifdef __cplusplus
//extern "C" 
//#endif
//void RHO_DARKSUSY_F77(double*, double*, double*, double*); //IMOS20060901

#define FJONES_F77 F77_FUNC(fjones,FJONES)
#ifdef __cplusplus
extern "C" 
#endif
double FJONES_F77(double*, double*, double*); // IMOS20060420

#define AIC_F77 F77_FUNC(aic,AIC)
#ifdef __cplusplus
extern "C" 
#endif
void AIC_F77(int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*); // IMOS20060420


//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// IMOS20060901
// a dummy routine to be replaced by DarkSUSY when compiled with GALPROP 

#define RHO_DARKSUSY_F77 F77_FUNC_(rho_darksusy,RHO_DARKSUSY)
#ifdef __cplusplus
extern "C" 
#endif
void RHO_DARKSUSY_F77(double* Xkpc, double* Ykpc, double* Zkpc, double* rho0) {

  *rho0 = -1.;
  
}

double rho_darksusy_cc(double Xkpc, double Ykpc, double Zkpc) {
   double rho;
   RHO_DARKSUSY_F77(&Xkpc, &Ykpc, &Zkpc, &rho);
   return rho;
}

double aic_cc(int key,int kbg,double E0,double E,double gamma,double RG,//IMOS20060420
            double rho,double xi,double z,double RS,double DENS)
{
//   cout<<"aic_cc"<<endl;
  double SPEC;
  if(z==0.) z = 1.e-4;
  AIC_F77(&key,&kbg,&E0,&E,&gamma,&RG,&rho,&xi,&z,&RS, &SPEC,&DENS);
  return(Pi*Rele*Rele/Mele *SPEC);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double fjones_cc(double gam, double E1, double E4)   //IMOS20060420
{
 // AWS20090128 documentation:
  // inverse Compton differential cross-section, cm^2 MeV^-1
  // gam=electron Lorentz factor
  // E1 = energy of    target photon (units of electron mass)
  // E4 = energy of scattered photon (units of electron mass)
 // FJONES_F77 gives the differential cross-section in cm^2 (E/Mele)^-1  / (Pi*Rele^2)
 // hence convert this to cm^2 MeV^-1
 //   cout<<"fjones_cc"<<endl;

  return ( Pi*Rele*Rele/Mele*FJONES_F77(&gam,&E1,&E4) );
}



