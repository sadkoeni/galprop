#include <Processes_Interface.h>

using namespace processes;

#include <config.h>
#include <map>



//Store a cache for the pp_meson routine, very useful for fitting.  Should not be a terrible
//slowdown for single operations.
//We expect key to vary slowly, then Esec, then Pp1, and finally NA1 and NA2
//Merge NA1 and NA2 into a single key by using NA1+1000*NA2
static std::map<int, //key
   std::map<double, //Esec
   std::map<double, //Pp1
   std::map<int, //NAs
   double> > > > pp_meson_cache;

#define NUCLEON_CS_F77 F77_FUNC_(nucleon_cs,NUCLEON_CS)
#ifdef __cplusplus
extern "C" 
#endif
void NUCLEON_CS_F77(int *option, double *Ek1, int *Zp1, int *Zt1, int *At1,
		    double *PP_inel, double *PA_inel, double *aPP_non, 
		    double *aPA_non, double *aPP_ann, double *aPA_ann) {

  double Ek = *Ek1;
  int Zp = *Zp1, Zt = *Zt1, At = *At1,opt = *option;
  nucleon_cs(opt,Ek,Zp,Zt,At,PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann);

}


#define PP_MESON_F77 F77_FUNC_(pp_meson, PP_MESON)
#ifdef __cplusplus
extern "C" 
#endif
double PP_MESON_F77(double*, double*, int*, int*, int*);

double processes::pp_meson_cc(double Esec, double Pp1, int NA1, int NA2, int key1) { 


   auto NAs = NA1 + 1000*NA2;

   double cs = -1;

#pragma omp critical (pp_meson_cache)
   {
      auto it_key = pp_meson_cache.find(key1);
      if (it_key != pp_meson_cache.end()) {
         auto it_Esec = it_key->second.find(Esec);

         if (it_Esec != it_key->second.end()) {
            auto it_Pp1 = it_Esec->second.find(Pp1);

            if (it_Pp1 != it_Esec->second.end()) {
               auto it_NAs = it_Pp1->second.find(NAs);

               if (it_NAs != it_Pp1->second.end()) {
                  cs = it_NAs->second;
               }
            }
         }
      }
   }

   if (cs >= 0)
      return cs;

   cs = PP_MESON_F77(&Esec, &Pp1, &NA1, &NA2, &key1);

#pragma omp critical (pp_meson_cache)
   {
      pp_meson_cache[key1][Esec][Pp1][NAs] = cs;
   }
   return cs;

}

#define SIGHAD_F77 F77_FUNC(sighad,SIGHAD)
#ifdef __cplusplus
extern "C" 
#endif
double SIGHAD_F77(int*, double*, double*, double*, double*, double*); // IMOS20020502

// Barashenkov & Polanski pA total cross section  IMOS20020502
double 
processes::sighad_cc(int IS, double PA, double PZ, double TA, double TZ, double E) { 
  
  return SIGHAD_F77(&IS, &PA, &PZ, &TA, &TZ, &E);

}

#define SIGTAP2_F77 F77_FUNC(sigtap2,SIGTAP2)
#ifdef __cplusplus
extern "C" 
#endif
void SIGTAP2_F77(int*, char*, int*); // IMOS20010511

// initialization of the Barashenkov & Polanski cross section code
void processes::sigtap_cc(int ISS, const std::string& path) {
 
  const std::string fullFilename = path + "/" + "barpol.dat";

  char* fn = const_cast<char*>(fullFilename.c_str());

  int length = fullFilename.size();

  SIGTAP2_F77(&ISS, fn, &length);

}

#define APRTAB_F77 F77_FUNC(aprtab,APRTAB)
#ifdef __cplusplus
extern "C" 
#endif
void APRTAB_F77(char*, int*); 

// initialization of the KACHELRIESS, MOSKALENKI and OSTAPCHENKO interpolation table
void processes::aprtab_cc(const std::string& path) {
 
  const std::string fullFilename = path + "/" + "ap-table.dat";

  char* fn = const_cast<char*>(fullFilename.c_str());

  int length = fullFilename.size();

  APRTAB_F77(fn, &length);

}

#define APSPEC_F77 F77_FUNC(apspec,APSPEC)
#ifdef __cplusplus
extern "C" 
#endif
double APSPEC_F77(double*, double*, int*, int*); 

// KACHELRIESS, MOSKALENKI and OSTAPCHENKO interpolation
double 
processes::apspec_cc(double e0, double epbar, int iap, int iat) { 
  
  return APSPEC_F77(&e0, &epbar, &iap, &iat);

}

#define ANTIPROTON_F77 F77_FUNC(antiproton,ANTIPROTON)
#ifdef __cplusplus
extern "C" 
#endif
double ANTIPROTON_F77(int*, double*, double*, int*, int*, int*, int*); // IMOS20010511

double 
processes::antiproton_cc(int key, double Pap1, double Pp1, int NZ1, int NA1, int NZ2, int NA2) { 

  return ANTIPROTON_F77(&key, &Pap1, &Pp1, &NZ1, &NA1, &NZ2, &NA2);

}

#define SYNCHROTRON_F77 F77_FUNC(synchrotron,SYNCHROTRON)
#ifdef __cplusplus
extern "C" 
#endif
double SYNCHROTRON_F77(double*, double*, double*);

double processes::synchrotron_cc(double gamma, double nu, double B) { 

   return SYNCHROTRON_F77(&gamma, &nu, &B);

}

#define BREMSS_SPEC_F77 F77_FUNC_(bremss_spec,BREMSS_SPEC)
#ifdef __cplusplus
extern "C" 
#endif
void BREMSS_SPEC_F77(double*, double*, int*, int*, double*);

double processes::bremss_spec_cc(double Egam, double E0, int IZ1, int Ne1) { 
  double dSdK;
  BREMSS_SPEC_F77(&Egam, &E0, &IZ1, &Ne1, &dSdK);
  return dSdK;
}

#define E_LOSS_COMPTON_F77 F77_FUNC_(e_loss_compton,E_LOSS_COMPTON)
#ifdef __cplusplus
extern "C" 
#endif
double E_LOSS_COMPTON_F77(double*, double*);

double processes::e_loss_compton_cc(double w, double gam) { 

  return E_LOSS_COMPTON_F77(&w, &gam);

}
