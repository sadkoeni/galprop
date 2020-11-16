
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * D_xx.cc *                                     galprop package * 02/13/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// Wave damping formalism is described in:
//
// Ptuskin, V.S., et al. 2006, ApJ 642, 902
// Ptuskin, V.S., et al. 2005, Adv. Space Res. 35, 162
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

using namespace std;//AWS20050624
#include<cstdio>
#include<cstdlib>
#include <sstream>

#include <gsl/gsl_sf_bessel.h>

#include"galprop_classes.h"
#include"galprop_internal.h"

static int iprotons,ir,ix,iy,iz,ip;
static int damping_min_ip; // IMOS20060330
static double damping_max_path_L;

//this is to avoid problems of using Galprop class members in static function "fu" IMOS20060322
static Particle *protons;
static double damping_p0;
static int n_spatial_dimensions, diff_reacc;
#pragma omp threadprivate(iprotons,ir,ix,iy,iz,ip,damping_min_ip,protons,damping_p0,n_spatial_dimensions,diff_reacc)

int Galprop::D_xx(Particle &particle,int iprotons_,int ir_,int ix_,int iy_,int iz_,int ip_)
{
   iprotons=iprotons_; ir=ir_; ix=ix_; iy=iy_; iz=iz_; ip=ip_;
   double L_cm, Lp_cm, tmp(0);
// integration parameters
   double a=particle.rigidity[ip],ai(0);

//this is to avoid problems of using Galprop class members in static function "fu" IMOS20060322
   protons=&gcr[iprotons];
   damping_p0=galdef.damping_p0;
   n_spatial_dimensions=galdef.n_spatial_dimensions;
   diff_reacc=galdef.diff_reacc;

   // Fix the diffusion coefficient to match the break rather then reference.
   double D0_xx=galdef.D0_xx; 
   if (galdef.D_rigid_ref != galdef.D_rigid_br) {
      if ( galdef.D_rigid_ref < galdef.D_rigid_br ) {
         D0_xx *= pow(galdef.D_rigid_br/galdef.D_rigid_ref, galdef.D_g_1);
      } else {
         D0_xx *= pow(galdef.D_rigid_br/galdef.D_rigid_ref, galdef.D_g_2);
      }
   }

   if (galdef.B_dep_diffusion>0) {
      D0_xx = Dxx_from_Bfield_model(D0_xx,ir,ix,iy,iz,ip);
   };

   // Apply plane dependencies
   if (galdef.Dxx_plane_scale > 0) 
      D0_xx *= (1.+(galdef.Dxx_plane_scale-1)*exp(-pow(particle.z[iz]/galdef.Dxx_plane_scale_height, 2)));


// STANDARD DIFFUSION COEFFICIENT (galdef.diff_reacc =0, 1, 2, -n==beta^n Dxx)
   if(galdef.diff_reacc < 3)
   {
// test of electron propagation vs analytical calculations IMOS20061030
     if(abs(galdef.DM_int0)==99) 
       particle.Dxx.d2[ir][iz].s[ip]=galdef.D0_xx *pow(particle.Ekin[ip]/galdef.D_rigid_br, galdef.D_g_1); 
       
// end of the test area
     else //IMOS20070110
       {

          double Dxx = pow(particle.beta[ip], galdef.D_eta) * D0_xx;
          if (galdef.diff_reacc<0) Dxx = pow(particle.beta[ip], galdef.diff_reacc) * D0_xx;
          if(particle.rigidity[ip]< galdef.D_rigid_br)
             Dxx *= pow(particle.rigidity[ip]/galdef.D_rigid_br, galdef.D_g_1);
          else
             Dxx *= pow(particle.rigidity[ip]/galdef.D_rigid_br, galdef.D_g_2);


         //Wave damping approximation using equation (15) of Ptuskin 2005 and b=3.  Assumes Kraichnan
         if (galdef.diff_reacc==2) {

            //Reset damping_min_ip when new grid point detected.  Really depends on the loop in fill_transport_array loop to be correct
            if(ip==particle.n_pgrid-1) {
               damping_min_ip=0; 
            }

            if ( ip < damping_min_ip ) {

               Dxx = particle.beta[ip] * C*galdef.damping_max_path_L/3.;

            } else {

               //Use g Psi_0(p) = k p**(3/2(1-b))
               //Most of damping done by protons, use their momentum
               const double omega = 2./3. * sqrt( galdef.damping_const_K * pow(protons->p[ip], -3./2.) );
               const double bess = gsl_sf_bessel_J1(2*omega);

               if ( bess > 0 ) {

                  Dxx *= omega / bess;

               } else {

                  Dxx = particle.beta[ip] * C*galdef.damping_max_path_L/3.;

                  damping_min_ip = ip;
               }
            }
         }

	 if(n_spatial_dimensions==2)
	   {
	     particle.Dxx.d2[ir][iz].s[ip] = Dxx;
	   }
	 if(n_spatial_dimensions==3)
	   {
	     particle.Dxx.d3[ix][iy][iz].s[ip] = Dxx;
	   }

       }
     return 0;
   }

// WAVE DAMPING (see Ptuskin et al. astro-ph/0301420)

   // Find out if we are at the solar location
   bool solar_location = false;
   if ( fabs(particle.z[iz]) < particle.dz/2. ) {
      if (n_spatial_dimensions==2)
         if ( fabs(particle.r[ir] - 8.49) < particle.dr/2. )
            solar_location = true;
      if (n_spatial_dimensions==3) 
         if ( fabs(particle.y[iy]) < particle.dy/2. )
            if ( fabs(particle.x[ix] - 8.49) < particle.dx/2. )
               solar_location = true;
   }

   //Reset damping_min_ip when new grid point detected.  Really depends on the loop in fill_transport_array loop to be correct
   if(ip==particle.n_pgrid-1) {
      damping_min_ip=0; 
      //This should not make a difference, this will be reset below before it is used
      damping_max_path_L=1e20;
   }

   L_cm  = damping_max_path_L;                               // max free path

   //TODO: This does not depend on anything but galdef.damping_p0 and the grid so should be moved out of the loop.
   //Possible destination is create_gcr
   int idp0(0);
   for(idp0=1;idp0<gcr[iprotons].n_pgrid;idp0++) 
      if(gcr[iprotons].rigidity[idp0]> galdef.damping_p0)
      {
	 damping_p0 = gcr[iprotons].rigidity[idp0];  // re-definition of galdef.damping_p0
         break;
      }

   //Free path according to the initial unaffected turbulence
   Lp_cm = 3./C *D0_xx*pow(particle.rigidity[ip]/galdef.D_rigid_br,galdef.D_g_1);

   if ( ip >= damping_min_ip ) {

      /*
      //Old method using numerical integrator with only half analytical.  Kept for reference only
      if(n_spatial_dimensions==2)
        if(iz>0 && iz<particle.n_zgrid-1 && ir<particle.n_rgrid-1)
          if(a<damping_p0) ai=sim(damping_p0,a,a/100.,0.01,1.e-10,&Galprop::fu);  //IMOS20060330
      if(n_spatial_dimensions==3)
        if(iz>0 && iz<particle.n_zgrid-1 && ix>0 && ix<particle.n_xgrid-1 && iy>0 && iy<particle.n_ygrid-1)
          if(a<damping_p0) ai=sim(damping_p0,a,a/100.,0.01,1.e-10,&Galprop::fu);  //IMOS20060330
      //   if(a<damping_p0) ai=sim(damping_p0,a,h,reps,aeps,&fu);  //IMOS20060330
      //   cout<<" Dxx------>"<<damping_p0<<" "<<a<<" "<<ai<<endl;
      ai = -ai; 
      */


      //New method using analytical integration assuming a piecewise broken power-law.  Must check that inner integral should go 
      //to damping_p0 only.  Only matters close to damping_p0, keep it that way
   
      //Nothing to do unless a < damping_p0
      if ( a < damping_p0 ) {
         //Lower index determined by a
         int ilow(0);
         for (ilow = 0; ilow < idp0; ++ilow)
            if (gcr[iprotons].p[ilow] > a)
               break;
         if (ilow > 0) --ilow;

         //Precalculate the power-law index and store the inner integral
         std::vector<double> plindex(idp0-ilow, 0.);
         double innerInt(0);

         // alpha + 1, never need alpha
         const double alphap1 = galdef.diff_reacc==11 ? 5./3. : 1.5;  //Kolmogorov vs. Kraichnan.

         if ( n_spatial_dimensions == 2 ) {

            for (int i(ilow); i < idp0; ++i)
               if ( protons->cr_density.d2[ir]    [iz].s[i] > 0 && protons->cr_density.d2[ir]    [iz].s[i+1] > 0 )
                  plindex[i-ilow] =log(gcr[iprotons].cr_density.d2[ir]    [iz].s[i]/gcr[iprotons].cr_density.d2[ir]    [iz].s[i+1])
                     /log(gcr[iprotons].                          p[i]/gcr[iprotons].                          p[i+1]);

            //Do the integration from top to bottom to be able to recursively calculate the inner integral
            for (int i(idp0-1); i >= ilow; --i) {
               if ( protons->cr_density.d2[ir]    [iz].s[i] > 0 && protons->cr_density.d2[ir]    [iz].s[i+1] > 0 ) {
                  innerInt += gcr[iprotons].cr_density.d2[ir]    [iz].s[i] / plindex[i-ilow] * pow(gcr[iprotons].p[i+1]/gcr[iprotons].p[i], plindex[i-ilow]);
                  //First one is different, second term in innerInt is non-existent
                  if ( i < idp0-1 )
                     innerInt -= gcr[iprotons].cr_density.d2[ir]    [iz].s[i+1] / plindex[i-ilow+1];

                  //First term in outer integral
                  if (i > ilow) {
                     ai += innerInt/alphap1 * ( pow(gcr[iprotons].p[i+1], alphap1) - pow(gcr[iprotons].p[i], alphap1) );
                     //second term
                     ai -= gcr[iprotons].cr_density.d2[ir]    [iz].s[i]/(plindex[i-ilow]*pow(gcr[iprotons].p[i],plindex[i-ilow])) / (alphap1 + plindex[i-ilow]) *
                        (pow(gcr[iprotons].p[i+1],plindex[i-ilow]+alphap1) - pow(gcr[iprotons].p[i],plindex[i-ilow]+alphap1));
                  } else {
                     ai += innerInt/alphap1 * ( pow(gcr[iprotons].p[i+1], alphap1) - pow(a, alphap1) );
                     //second term
                     ai -= gcr[iprotons].cr_density.d2[ir]    [iz].s[i]/(plindex[i-ilow]*pow(gcr[iprotons].p[i],plindex[i-ilow])) / (alphap1 + plindex[i-ilow]) *
                        (pow(gcr[iprotons].p[i+1],plindex[i-ilow]+alphap1) - pow(a,plindex[i-ilow]+alphap1));
                  }
               }

               //No need to go further 
               if (galdef.damping_const_G*ai > 1)
                  break;
            }

         }

         if ( n_spatial_dimensions == 3 ) {

            for (int i(ilow); i < idp0; ++i)
               if ( protons->cr_density.d3[ix][iy][iz].s[i] > 0 && protons->cr_density.d3[ix][iy][iz].s[i+1] > 0 )
                  plindex[i-ilow] =log(gcr[iprotons].cr_density.d3[ix][iy][iz].s[i]/gcr[iprotons].cr_density.d3[ix][iy][iz].s[i+1])
                     /log(gcr[iprotons].                          p[i]/gcr[iprotons].                          p[i+1]);

            //Do the integration from top to bottom to be able to recursively calculate the inner integral
            for (int i(idp0-1); i >= ilow; --i) {
               if ( protons->cr_density.d3[ix][iy][iz].s[i] > 0 && protons->cr_density.d3[ix][iy][iz].s[i+1] > 0 ) {
                  innerInt += gcr[iprotons].cr_density.d3[ix][iy][iz].s[i] / plindex[i-ilow] * pow(gcr[iprotons].p[i+1]/gcr[iprotons].p[i], plindex[i-ilow]);
                  //First one is different, second term in innerInt is non-existent
                  if ( i < idp0-1 )
                     innerInt -= gcr[iprotons].cr_density.d3[ix][iy][iz].s[i+1] / plindex[i-ilow+1];

                  if ( i > ilow) {
                     //First term in outer integral
                     ai += innerInt/alphap1 * ( pow(gcr[iprotons].p[i+1], alphap1) - pow(gcr[iprotons].p[i], alphap1) );
                     //second term
                     ai -= gcr[iprotons].cr_density.d3[ix][iy][iz].s[i]/(plindex[i-ilow]*pow(gcr[iprotons].p[i],plindex[i-ilow])) / (alphap1 + plindex[i-ilow]) *
                        (pow(gcr[iprotons].p[i+1],plindex[i-ilow]+alphap1) - pow(gcr[iprotons].p[i],plindex[i-ilow]+alphap1));
                  } else {
                     //First term in outer integral
                     ai += innerInt/alphap1 * ( pow(gcr[iprotons].p[i+1], alphap1) - pow(a, alphap1) );
                     //second term
                     ai -= gcr[iprotons].cr_density.d3[ix][iy][iz].s[i]/(plindex[i-ilow]*pow(gcr[iprotons].p[i],plindex[i-ilow])) / (alphap1 + plindex[i-ilow]) *
                        (pow(gcr[iprotons].p[i+1],plindex[i-ilow]+alphap1) - pow(a,plindex[i-ilow]+alphap1));
                  }
               }

               //No need to go further 
               if (galdef.damping_const_G*ai > 1)
                  break;
            }

         }

         if(solar_location) {
            std::ostringstream oss;
            oss << "Index between "<<ilow<<" and "<<ilow+1<<" is: "<<plindex[0];
            DEBUGLOG(oss.str());
         }
      }


      tmp = 1;

      // Kolmogorov diffusion with wave damping ## 
      if(galdef.diff_reacc==11)
      {
         tmp = 1. -galdef.damping_const_G*ai; //    /pow(particle.Z, 5./3.)
         if(tmp>0.) L_cm = Lp_cm/pow(tmp, 2);
      }

      // Kraichnan diffusion with wave damping ## 
      if(galdef.diff_reacc==12)
      {
         tmp = 1. -galdef.damping_const_G*ai; //    /pow(particle.Z, 3./2.)
         if(tmp>0.) L_cm = Lp_cm/tmp;
      }


      //Fix the free path to a maximum value as soon as wave damping becomes too strong
      if (a < damping_p0 && L_cm >= 0.5*Lp_cm/galdef.damping_const_G) 
      { 

         std::ostringstream oss;
         oss<<" New maximum free path: "<<ir<<" "<<iz<<" "<<ip<<" "<<damping_max_path_L<<" "<<Lp_cm/galdef.damping_const_G<<" "<<tmp<<" "<<L_cm;
         DEBUGLOG(oss.str());

         damping_max_path_L = 0.5*Lp_cm/galdef.damping_const_G;
         L_cm = damping_max_path_L;

         damping_min_ip=ip;
      }

   }

   const double newDxx = particle.beta[ip] *C*L_cm/3.;

   if(solar_location) {
      std::ostringstream oss;
      if ( n_spatial_dimensions == 2 ) {
         oss<<" D_xx>>>> "<<particle.rigidity[ip]<<" "<<particle.Ekin[ip]
            <<" "<<newDxx<<" "<<ai<<" "<<1. - particle.Dxx.d2[ir]    [iz].s[ip]/newDxx;
      }
      if ( n_spatial_dimensions == 3 ) {
         oss<<" D_xx>>>> "<<particle.rigidity[ip]<<" "<<particle.Ekin[ip]
            <<" "<<newDxx<<" "<<ai<<" "<<1. - particle.Dxx.d3[ix][iy][iz].s[ip]/newDxx;
      }
      DEBUGLOG(oss.str());
   }

   // Use damping during convergenge.  
   // Have to make sure we only apply it during the convergenge iterations, this is handled in propagate particles by nulling the Dxx distribution.
   // logarithmic fraction of D_xx change to apply each iteration to damp the convergence and make it more stable and faster
   const double frac = 0.2;

   if(n_spatial_dimensions==2) {
      if ( particle.Dxx.d2[ir]    [iz].s[ip] > 0 ) {
         if (fabs(log(newDxx/particle.Dxx.d2[ir]    [iz].s[ip])) > 0.001*frac) {
            particle.Dxx.d2[ir]    [iz].s[ip] = pow(newDxx, frac)*pow(particle.Dxx.d2[ir]    [iz].s[ip], 1-frac);
         } else {
            particle.Dxx.d2[ir]    [iz].s[ip] = newDxx;
         }

      } else {
         particle.Dxx.d2[ir]    [iz].s[ip] = newDxx;
      }
   }
   if(n_spatial_dimensions==3) {
      if ( particle.Dxx.d3[ix][iy][iz].s[ip] > 0 ) {
         if (fabs(log(newDxx/particle.Dxx.d3[ix][iy][iz].s[ip])) > 0.001*frac) {
            particle.Dxx.d3[ix][iy][iz].s[ip] = pow(newDxx, frac)*pow(particle.Dxx.d3[ix][iy][iz].s[ip], 1-frac);
         } else {
            particle.Dxx.d3[ix][iy][iz].s[ip] = newDxx;
         }

      } else {
         particle.Dxx.d3[ix][iy][iz].s[ip] = newDxx;
      }
   }
   return 0;
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

double Galprop::fu(double x)
{
   int i,n,m;
   double int_psi=0., y;

//     cout<<" ir,iz,ip,iprotons,x= "<<ir<<" "<<iz<<" "<<ip<<" "<<iprotons<<" "<<x<<endl;
//     for(int i=0; i<protons->n_pgrid; i++) cout<<" "<<protons->cr_density.d2[ir][iz-1].s[i];
//     cout<<endl;

// search in the grid
   for(n=1;n<protons->n_pgrid;n++) if(protons->p[n]> x) break;
   for(m=1;m<protons->n_pgrid;m++) if(protons->p[m]> damping_p0) break;
//   cout<<">>>> x,n= "<<x<<" "<<n<<" gcr[iprotons].p[n]= "<<gcr[iprotons].p[n]<<endl;

// integration over rigidity (= momentum for protons)

   if(n_spatial_dimensions==2)
   { 
// fit each interval with power-law and integrate analytically
      for(int_psi=0., i=n-1; i<m; i++)
      {
	 if ( protons->cr_density.d2[ir]    [iz].s[i] > 0 && protons->cr_density.d2[ir]    [iz].s[i+1] > 0 )
	 {
// derive y (= power-law index) 
         y =log(protons->cr_density.d2[ir]    [iz].s[i]/protons->cr_density.d2[ir]    [iz].s[i+1])
           /log(protons->                          p[i]/protons->                          p[i+1]);
// integrate (psi/p) dp
         if(i>n-1) int_psi +=protons->cr_density.d2[ir][iz].s[i  ]       // fit norm.
		        *pow(protons-> p[i  ],-y)/y
                       *(pow(protons-> p[i+1], y)
                        -pow(protons-> p[i  ], y));
         else      int_psi +=protons->cr_density.d2[ir][iz].s[i  ]       // fit norm.
                        *pow(protons-> p[i  ],-y)/y
                       *(pow(protons-> p[i+1], y)
                        -pow(x,               y));
	 }
      }
   }

   if(n_spatial_dimensions==3)
   { 
// fit each interval with power-law and integrate analytically
      for(int_psi=0., i=n-1; i<m; i++)
      {
// derive y (= power-law index) 
	 if ( protons->cr_density.d3[ix][iy][iz].s[i] > 0 && protons->cr_density.d3[ix][iy][iz].s[i+1] > 0 )
	 {
         y =log(protons->cr_density.d3[ix][iy][iz].s[i]/protons->cr_density.d3[ix][iy][iz].s[i+1])
           /log(protons->                          p[i]/protons->                          p[i+1]);
// integrate (psi/p) dp
         if(i>n-1) int_psi +=protons->cr_density.d3[ix][iy][iz].s[i  ]       // fit norm.
                        *pow(protons-> p[i  ],-y)/y
                       *(pow(protons-> p[i+1], y)
                        -pow(protons-> p[i  ], y));
         else      int_psi +=protons->cr_density.d3[ix][iy][iz].s[i  ]       // fit norm.
                        *pow(protons-> p[i  ],-y)/y
                       *(pow(protons-> p[i+1], y)
                        -pow(x,               y));
	 }
      }
   }

//if(ir==0 && iz==80 && ip==40) cout<<" integral = "<<n<<" "<<y<<" "<<int_psi<<" "<<protons->p[n]<<" "<<x<<" "<<damping_p0 <<endl; //exit(1);
   if(diff_reacc==11) return ( pow(x,2./3.)*int_psi );   // Kolmogorov
   if(diff_reacc==12) return (sqrt(x)      *int_psi );   // Kraichnan
   return 0; // in case of an error return 0.
}


double Galprop::Dxx_from_Bfield_model(double d0, int ir, int ix, int iy, int iz, int ip){
    
    double bRandom(0),bRegular(0),bRandomSolar(0),bRegularSolar(0);
#ifdef DEBUG    
    static int oldir,oldiz;
#endif
    const int zSolar=int((1e-5-galaxy.z_min)/galaxy.dz);
    
    if (galdef.n_spatial_dimensions==2){
        const int rSolar=int((8.5-galaxy.r_min)/galaxy.dr);
        bRandom=galaxy.B_field_random.d2[ir][iz].s[0];
        bRegular=galaxy.B_field_regular.d2[ir][iz].s[0];
        bRandomSolar=galaxy.B_field_random.d2[rSolar][zSolar].s[0];
        bRegularSolar=galaxy.B_field_regular.d2[rSolar][zSolar].s[0];
    }
    else if (galdef.n_spatial_dimensions==3){
        const int xSolar=int((8.5-galaxy.x_min)/galaxy.dx);
        const int ySolar=int((1.e-5-galaxy.y_min)/galaxy.dy);
        bRandom=galaxy.B_field_random.d3[ix][iy][iz].s[0];
        bRegular=galaxy.B_field_regular.d3[ix][iy][iz].s[0];
        bRandomSolar=galaxy.B_field_random.d3[xSolar][ySolar][zSolar].s[0];
        bRegularSolar=galaxy.B_field_regular.d3[xSolar][ySolar][zSolar].s[0];
    };

    if (bRandom == 0) {
       WARNING("Random field 0, cannot scale Dxx");
       return d0;
    }
    if (bRandomSolar == 0) {
       WARNING("Random field at solar location is 0, cannot scale Dxx");
       return d0;
    }
    if (bRegular == 0) {
       WARNING("Regular field is 0, cannot scale Dxx");
       return d0;
    }
    if (bRegularSolar == 0) {
       WARNING("Regular field at solar location is 0, cannot scale Dxx");
       return d0;
    }
    
    double D_xx_max = galdef.D_xx_max;
    double D_xx_min = galdef.D_xx_min;
    // Transform it to be at reference rigidity
    if (galdef.D_rigid_br < 1e4) {
       D_xx_max *= pow(galdef.D_rigid_br/1e4, galdef.D_g_2);
       D_xx_min *= pow(galdef.D_rigid_br/1e4, galdef.D_g_2);
    } else {
       D_xx_max *= pow(galdef.D_rigid_br/1e4, galdef.D_g_1);
       D_xx_min *= pow(galdef.D_rigid_br/1e4, galdef.D_g_1);
    }

    double dxx=d0;
    switch (galdef.B_dep_diffusion) {
       // model: diffusion coefficient scales like (B/dB)^2 r_G^gamma = (B/dB)^2 (R/B)^gamma
       // use only the D_g_1 index for the magnetic field scaling, since a break would cause a jump in the diffusion coefficient at the break energy in this model--> unphysical  
       // use D_g_1 to be consistent with wave damping
       // limit diffusion coefficient @ reference rigidity to < D_xx_max at 10 GV.  Defaults to 10^30, corresponding to 30 pc mean free path.
       // magnetic field dependence is normalized so that Dxx(r_Solar,z_Solar)=galdef.D_xx
       case 1:  
          dxx=d0*pow(bRegular/bRegularSolar,2.-galdef.D_g_1)*pow(bRandomSolar/bRandom,2);
          dxx=std::min(dxx,D_xx_max);
          dxx=std::max(dxx,D_xx_min);
          break;
       default:
          break;
    };
#ifdef DEBUG
    if(ir!=oldir || iz!=oldiz) std::cout<<"D0_xx: r= "<<galaxy.r[ir]<<" z= "<<galaxy.z[iz]<<" D0_xx= "<<d0
                         <<" bRandom= "<<(1e10*bRandom)<<" uG bRegular= "<<(1e10*bRegular)
                         <<" uG D0_new= "<<dxx <<std::endl;
    oldir=ir;oldiz=iz;
#endif    
    
    return dxx;
};

