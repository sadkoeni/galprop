
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * D_pp.cc *                                     galprop package * 2/13/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// diffusion in momentum space
// formula from Seo and Ptuskin ApJ 431, 705 (1994); Ptuskin-2003 (Berezinskii et al. book)
// with w=1 (since can be subsumed in v_alfven)
// v_alfven in km s-1
// NB Dpp_constant was defined but not used in galprop v25 and earlier

using namespace std;//AWS20050624
#include<iostream>
#include<cmath>
#include"galprop_classes.h"

double Galprop::D_pp(double p,double a,double v_alfven,double D_xx)
{
   const double Dpp_constant = 4./(3*a*(4.-a)*(4.-a*a));
   //if(galdef.diff_reacc==1) Dpp_constant=2.*2./(3.*(2.-a)*(4.-a)*(2.+a)*a); // Seo & Ptuskin
   //else                     Dpp_constant=1./(a*(4.-a)*(4.-a*a));            // IMOS20030213 Ptuskin-2003
   
   const double D_pp_=Dpp_constant*pow(p,2.0) * pow(v_alfven*1.e5,2) /D_xx; 

// cout<<"D_pp "<<D_pp_<<" "<<a<<" "<<Dpp_constant<<endl;  
   return D_pp_;
}


double Galprop::alfvenVelocity(int ir,int ix, int iy, int iz){
#ifdef DEBUG   
   static int oldir,oldiz;
#endif   
    const int zSolar=int((1e-5-galaxy.z_min)/galaxy.dz);
   
   if (galdef.B_dep_diffusion==0) return galdef.v_Alfven;
   
   double Bfield(0),gasDensity(0),BfieldSolar(0),gasDensitySolar(0);
   if(galdef.n_spatial_dimensions==2){
        const int rSolar=int((8.5-galaxy.r_min)/galaxy.dr);
        Bfield=galaxy.B_field.d2[ir][iz].s[0];
        gasDensity=galaxy.n_HII.d2[ir][iz].s[0];
        BfieldSolar=galaxy.B_field.d2[rSolar][zSolar].s[0];
        gasDensitySolar=galaxy.n_HII.d2[rSolar][zSolar].s[0];
   }
   else if (galdef.n_spatial_dimensions==3){
        const int xSolar=int((8.5-galaxy.x_min)/galaxy.dx);
        const int ySolar=int((1.e-5-galaxy.y_min)/galaxy.dy);
        Bfield=galaxy.B_field.d3[ix][iy][iz].s[0];
        gasDensity=galaxy.n_HII.d3[ix][iy][iz].s[0];
        BfieldSolar=galaxy.B_field.d3[xSolar][ySolar][zSolar].s[0];
        gasDensitySolar=galaxy.n_HII.d3[xSolar][ySolar][zSolar].s[0];
   };

   //Correct the gas density for the filling factor from Gaensler 2008 independent of nHII model selected
   //This is because there is no explicit equation for the filling factor for the other 2 models and
   //it is not possible to determine it from the data given in the papers.
   gasDensity /= 0.04 * exp(abs(galaxy.z[iz])/0.7);
   gasDensitySolar /= 0.04;
   
   gasDensity+=0.0025*exp(-abs(galaxy.z[iz])/4.4);   // thin hot ionized medium with large scale height from Kerp & Kalberla 1998
   gasDensitySolar+=0.0025;
   
   const double vAlfven= Bfield/BfieldSolar / sqrt(gasDensity/gasDensitySolar);
   const double scaled_vAlfven=galdef.v_Alfven*vAlfven;
   
#ifdef DEBUG   
   if(ir!=oldir || iz!=oldiz)
   std::cout<<"alfvenVelocity: r= "<<galaxy.r[ir]<<" z= "<<galaxy.z[iz]<<" Bfield= "<<Bfield<<" gasDensity= "<<gasDensity
            <<" vAlfven= "<<vAlfven<<" scaled_vAlfven= "<<scaled_vAlfven<<std::endl;
   if (ir<oldir) exit(1);
   oldir=ir; oldiz=iz;
#endif

   return scaled_vAlfven;
};
