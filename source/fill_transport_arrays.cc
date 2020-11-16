
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * create_transport_arrays.cc *                  galprop package * 10/12/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"
#include <cstring>

#include <ErrorLogger.h>
#include <Processes_Interface.h>
#include <Nuclei_Interface.h>

#include "constants.h"

#include <string>
#include <sstream>

static int key=-1;

int Galprop::fill_transport_arrays(Particle& particle) {
   
  INFO("Entry");
  int stat=0, A1,Z2,A2,K_electron; // IMOS20010816
  int galdef_network_par=0;  // imos network, use everywhere IMOS20010816
  double t_half[2]; // IMOS20010816
  double fragment_p,fragment_He,PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann; // IMOS20000607
   
  // ASSIGNING PRIMARY SOURCE FUNCTION
  
  INFO("Assigning primary source function");
  particle.primary_source_function = 0;
  particle.secondary_source_function = 0;
  std::for_each(galdef.source_classes.begin(), galdef.source_classes.end(), [&](const std::unique_ptr<SourceClass> &sc) { sc->addSource(particle); });

  //Multiply with the normalization
  particle.primary_source_function *= particle.normalization_factor;
  
  //Loop over the source distribution and apply cuts on rigidity and local bubble.
  const bool local_bubble(galdef.local_bubble_radius > 0 && galdef.local_bubble_source_fraction != 1.0) ;
  if (  local_bubble ||
        particle.rigidity[0] <= galdef.rigid_min ||
        particle.rigidity[particle.rigidity.size()-1] > galdef.rigid_max ||
        particle.Ekin[0] <= galdef.inj_Ekin_min ||
        particle.Ekin[particle.Ekin.size()-1] > galdef.inj_Ekin_max )
  {
     size_t ipmin(0);
     size_t ipmax(particle.rigidity.size());

     while ( ( particle.rigidity[ipmin] <= galdef.rigid_min || particle.Ekin[ipmin] <= galdef.inj_Ekin_min ) && ipmin < particle.rigidity.size()-1 )
        ++ipmin;
     while ( ( particle.rigidity[ipmax-1] > galdef.rigid_max || particle.Ekin[ipmax-1] <= galdef.inj_Ekin_max ) && ipmax > 1 )
        --ipmax;

     if ( particle.n_spatial_dimensions == 2 ) {

        for ( int ir(0); ir < particle.n_rgrid; ++ir ) {
           for ( int iz(0); iz < particle.n_zgrid; ++iz ) {

              //0 everything above and below the given energy range
              for ( size_t ip(0); ip < ipmin; ++ip ) 
                 particle.primary_source_function.d2[ir][iz].s[ip] = 0;
              for ( size_t ip(ipmax); ip < particle.rigidity.size(); ++ip ) 
                 particle.primary_source_function.d2[ir][iz].s[ip] = 0;

              //multiply within local bubble
              if ( local_bubble ) {
                 const double dSun = sqrt( (particle.r[ir]-Rsun)*(particle.r[ir]-Rsun) + particle.z[iz]*particle.z[iz] );

                 if (dSun < galdef.local_bubble_radius) 
                    for ( size_t ip(ipmin); ip < ipmax; ++ip ) 
                       particle.primary_source_function.d2[ir][iz].s[ip] *= galdef.local_bubble_source_fraction;
              }
           }
        }
     }

     if ( particle.n_spatial_dimensions == 3 ) {

#pragma omp parallel for schedule(static) default(shared)
        for ( int ix=0; ix < particle.n_xgrid; ++ix ) {
           for ( int iy(0); iy < particle.n_ygrid; ++iy ) {
              for ( int iz(0); iz < particle.n_zgrid; ++iz ) {

                 //0 everything above and below the given energy range
                 for ( size_t ip(0); ip < ipmin; ++ip ) 
                    particle.primary_source_function.d3[ix][iy][iz].s[ip] = 0;
                 for ( size_t ip(ipmax); ip < particle.rigidity.size(); ++ip ) 
                    particle.primary_source_function.d3[ix][iy][iz].s[ip] = 0;

                 //multiply within local bubble
                 if ( local_bubble ) {
                    const double dSun = sqrt( (particle.x[ix]-Rsun)*(particle.x[ix]-Rsun) + particle.y[iy]*particle.y[iy] + particle.z[iz]*particle.z[iz] );

                    if (dSun < galdef.local_bubble_radius) 
                       for ( size_t ip(ipmin); ip < ipmax; ++ip ) 
                          particle.primary_source_function.d3[ix][iy][iz].s[ip] *= galdef.local_bubble_source_fraction;
                 }
              }
           }
        }
     }
  }

  if (galdef.verbose == -123) {
     std::ostringstream os;
     os<<particle.name<<" Primary source function:";
     INFO(os.str());
     particle.primary_source_function.print();
  }

  // ASSIGNING DIFFUSION COEFFICIENT
  
  INFO("      assigning diffusion coefficient");
  
  ////////////////////////V//IMOS20030214 all region
  int iprotons=-1;
  
  std::ostringstream ost;
  if(galdef.diff_reacc > 1)
    {
      // identify CR protons
      for(int i=n_species-1; i>=0; --i)  
	if(101==100*gcr[i].Z+gcr[i].A)
	  {
            iprotons=i;
	    ost.str("");
	    ost<<"  CR protons found as species #"<<iprotons;
	    INFO(ost.str());
	    break;
	  }
      if(iprotons==-1) { WARNING("CR protons not found!"); return 1; }
      ost.str("");
      ost<<"create_transport_arrays>> "<<particle.Z*100+particle.A<<" "<<particle.p[0];
      INFO(ost.str());
      
      // Zero approximation proton spectrum for calculation of damping
      /*
      if(gcr[iprotons].cr_density.max() == 0)
	{
	  if(galdef.n_spatial_dimensions==2)
            for(int ir=0; ir<particle.n_rgrid; ir++)
	      for(int iz=0; iz<particle.n_zgrid; iz++)
		for(int ip=0; ip<particle.n_pgrid; ip++)
		  {
		    gcr[iprotons].cr_density.d2[ir]    [iz].s[ip] = 
		      (1.-(galdef.z_min+iz*galdef.dz)/galdef.z_max)
		      *galdef.proton_norm_flux *pow(gcr[iprotons].Etot[ip]/galdef.proton_norm_Ekin, -2.75);
		    //cout<<"create_transport_arrays>> "<<gcr[iprotons].cr_density.d2[ir]    [iz].s[ip]<<endl;
                  }
	  
	  if(galdef.n_spatial_dimensions==3)
            for(int ix=0; ix<particle.n_xgrid; ix++)
	      for(int iy=0; iy<particle.n_ygrid; iy++)
		for(int iz=0; iz<particle.n_zgrid; iz++)
		  for(int ip=0; ip<particle.n_pgrid; ip++)
		    gcr[iprotons].cr_density.d3[ix][iy][iz].s[ip] = 
		      (1.-(galdef.z_min+iz*galdef.dz)/galdef.z_max)
		      *galdef.proton_norm_flux *pow(gcr[iprotons].Etot[ip]/galdef.proton_norm_Ekin, -2.75);
	}
        */
      key=1;
    }
  //   cout<<"create_transport_arrays>> "<<iprotons<<" "<<gcr[iprotons].cr_density.d2[9][9].s[10]<< " "<<protons.d2[9][9].s[10]<<endl;
  //   for(int ip=0; ip<particle.n_pgrid; ip++) cout<<" "<<protons.d2[9][9].s[ip];
  //   cout<<endl;
  ////////////////////////^//IMOS20030214

  
  if(galdef.n_spatial_dimensions==2)
    {
      ost.str("");
      ost<<" Calculating diffusion coefficients (ir, iz, ip): ("<<particle.n_rgrid<<", "<<particle.n_zgrid<<", "<<particle.n_pgrid<<")";
      INFO(ost.str());

      //Store Dxx distribution for debugging info in wavedamping
      Distribution D_xx_old(particle.Dxx);

      //TODO: possibly add some precalculations for the wave damping sum, like the index, inner integral, etc.
//#pragma omp parallel for schedule(dynamic) default(shared)
      for(int ir=0; ir<particle.n_rgrid-1; ++ir)
	for(int iz=1; iz<particle.n_zgrid-1; ++iz)
	  for(int ip=particle.n_pgrid-1; ip>=0; --ip) // IMOS20060330 changed to reverse order (Wave-particle interactions)
	    {
	      D_xx(particle,iprotons,ir, 0, 0,iz,ip); //array assigned in D_xx IMOS20030129
	      particle.Dpp.d2[ir][iz].s[ip] =D_pp(particle.p[ip],galdef.D_g_1,alfvenVelocity(ir,0,0,iz),particle.Dxx.d2[ir][iz].s[ip]);
	    }

      //Assign the boundary of D_xx
      for(int iz=0; iz<particle.n_zgrid; ++iz)
         for(int ip=particle.n_pgrid-1; ip>=0; --ip) // IMOS20060330 changed to reverse order (Wave-particle interactions)
         {
            particle.Dxx.d2[particle.n_rgrid-1][iz].s[ip] = particle.Dxx.d2[particle.n_rgrid-2][iz].s[ip];
            particle.Dpp.d2[particle.n_rgrid-1][iz].s[ip] = D_pp(particle.p[ip],galdef.D_g_1,alfvenVelocity(particle.n_rgrid-1,0,0,iz),particle.Dxx.d2[particle.n_rgrid-1][iz].s[ip]);
         }
      for(int ir=0; ir<particle.n_rgrid; ++ir)
         for(int ip=particle.n_pgrid-1; ip>=0; --ip) // IMOS20060330 changed to reverse order (Wave-particle interactions)
         {
            particle.Dxx.d2[ir][0].s[ip] = particle.Dxx.d2[ir][1].s[ip];
            particle.Dpp.d2[ir][0].s[ip] = D_pp(particle.p[ip],galdef.D_g_1,alfvenVelocity(ir,0,0,0),particle.Dxx.d2[ir][0].s[ip]);
            particle.Dxx.d2[ir][particle.n_zgrid-1].s[ip] = particle.Dxx.d2[ir][particle.n_zgrid-2].s[ip];
            particle.Dpp.d2[ir][particle.n_zgrid-1].s[ip] = D_pp(particle.p[ip],galdef.D_g_1,alfvenVelocity(ir,0,0,particle.n_zgrid-1),particle.Dxx.d2[ir][particle.n_zgrid-1].s[ip]);
         }

      
      //Calculate the maximum relative error for information purposes
      //TODO: Use something along this lines as criteria for stopping the wave damping iterations
      double maxRelErr(0);
      for(int ir=0; ir<particle.n_rgrid; ++ir)
	for(int iz=0; iz<particle.n_zgrid; ++iz)
	  for(int ip=particle.n_pgrid-1; ip>=0; --ip) 
          {
             const double relErr = (particle.Dxx.d2[ir][iz].s[ip] - D_xx_old.d2[ir][iz].s[ip])/particle.Dxx.d2[ir][iz].s[ip];
             if (fabs(relErr) > fabs(maxRelErr))
                maxRelErr = relErr;
          }
      ost.str("");
      ost<<"Maximum relative error in D_xx: "<<maxRelErr;
      INFO(ost.str());
     

    }
  if(galdef.n_spatial_dimensions==3)
    {
      ost.str("");
      ost<<" Calculating diffusion coefficients (ix, iy, iz, ip): ("<<particle.n_xgrid<<", "<<particle.n_ygrid<<", "<<particle.n_zgrid<<", "<<particle.n_pgrid<<")";
      INFO(ost.str());
#pragma omp parallel for schedule(static) default(shared)
      for(int ix=0; ix<particle.n_xgrid; ++ix)
	for(int iy=0; iy<particle.n_ygrid; ++iy)
	  for(int iz=0; iz<particle.n_zgrid; ++iz)
	    for(int ip=particle.n_pgrid-1; ip>=0; --ip) // IMOS20060330 changed to reverse order
	      {
		D_xx(particle,iprotons, 0,ix,iy,iz,ip); //array assigned in D_xx IMOS20030129
		particle.Dpp.d3[ix][iy][iz].s[ip] =D_pp(particle.p[ip], galdef.D_g_1, alfvenVelocity(0,ix,iy,iz), particle.Dxx.d3[ix][iy][iz].s[ip]);
	      }  // p
    } 
  if(galdef.verbose>=2)
    {
      ost.str("");
      ost<< "spatial   diffusion coefficient Dxx  for species "<<particle.name;
      INFO(ost.str());
      particle.Dxx.print();
      ost.str("");
      ost<< "momentum diffusion coefficient Dpp  for species "<<particle.name;
      INFO(ost.str());
      particle.Dpp.print();
    }
  
  // ASSIGNING FRAGMENTATION RATE

   INFO("======== assigning fragmentation rate ======== ");
   int ZH=1, ZHe=2;                                                    //  IMOS20010816
   double attach_H=0., attach_He=0., strip_H=0., strip_He=0.;          //  IMOS20010816

   particle.fragment=0.0;
   if(particle.A!=0)
   {
     double CSratio,CStot_ratio;
     
     for(int ip=0; ip<particle.n_pgrid; ip++)
       {                                                             // IMOS20000607 whole segment
	 A1 = 1;                                                    // nucleus
	 Z2 = particle.Z;  A2 = particle.A;                         // - " -
	 if(101==100*Z2+A2) { Z2 = 2;  A2 = 4; }                    // protons
	 if(-99==100*Z2+A2) { A1 =-1;  Z2 = 2;  A2 = 4; }           // antiprotons
	 processes::nucleon_cs(galdef.total_cross_section,particle.Ekin[ip]*1.e-3,A1,Z2,A2, // AWS20010620
		    &PP_inel,&PA_inel,&aPP_non,&aPA_non,&aPP_ann,&aPA_ann);
	 nuclei::He_to_H_CS(particle.Ekin[ip]*1.e-3,particle.Z,particle.A,999,999,&CSratio,&CStot_ratio);
	 
	 fragment_p = PA_inel;                                      // nuclei
	 fragment_He= PA_inel*CStot_ratio;                          // -"- 

// ELECTRON ATTACHMENT/STRIPPING CROSS SECTION                                  IMOS20010816
	 
	 if(galdef.K_capture)
	   {
	     for(K_electron=0;K_electron<=galdef.K_capture;K_electron++) 
	       nucdata(galdef_network_par,particle.Z,particle.A,K_electron,particle.Z,particle.A,&Z2,&A2,&t_half[K_electron]);
	   
	     if(t_half[0] != t_half[1])
	       {
		 nuclei::Kcapture_cs(particle.Ekin[ip],particle.Z,ZH, &attach_H ,&strip_H );
		 nuclei::Kcapture_cs(particle.Ekin[ip],particle.Z,ZHe,&attach_He,&strip_He);
		 if(particle.K_electron)
		   {
                     fragment_p += strip_H ;
                     fragment_He+= strip_He;
		   } else
		     {
		       fragment_p += attach_H ;
		       fragment_He+= attach_He;
		     }
		 
		 if(galdef.verbose==-502)// selectable debug
		 {
		   ost.str("");
		   ost<<"create_transport_arrays: Ekin Z,A,K_electron,attach_H,strip_H: "
		       <<particle.Ekin[ip]<<" "<<particle.Z<<" "<<particle.A<<" "
		       <<particle.K_electron<<" "<<strip_H<<" "<<attach_H;
		   INFO(ost.str());
		 }
               }
	   }
	 if(101==100*particle.Z+particle.A)                         // protons
	   { 
	     fragment_p = PP_inel;
	     fragment_He= PA_inel;
	   }
	 if(-99==100*particle.Z+particle.A)                         // antiprotons
	   { 
	     fragment_p = aPP_non+aPP_ann;
	     fragment_He= aPA_non+aPA_ann;
	   }

	 if(galdef.n_spatial_dimensions==2)
	   for(int ir=0; ir<particle.n_rgrid; ir++)                   // IMOS20010816
	     for(int iz=0; iz<particle.n_zgrid; iz++)
	       particle.fragment.d2[ir] [iz].s[ip]= particle.beta[ip]*C
		 *(galaxy.n_HI.d2[ir] [iz].s[0] +2*galaxy.n_H2.d2[ir] [iz].s[0]+galaxy.n_HII.d2[ir] [iz].s[0]) 
		 *(fragment_p+galdef.He_H_ratio*fragment_He) *1.0e-27;
     	 
	 if(galdef.n_spatial_dimensions==3)
#pragma omp parallel for schedule(static) default(shared)
	   for(int ix=0; ix<particle.n_xgrid; ix++)                // IMOS20010816
	     for(int iy=0; iy<particle.n_ygrid; iy++)
	       for(int iz=0; iz<particle.n_zgrid; iz++)
		 particle.fragment.d3[ix][iy][iz].s[ip]= particle.beta[ip]*C
		   *(galaxy.n_HI.d3[ix][iy][iz].s[0] +2*galaxy.n_H2.d3[ix][iy][iz].s[0]+galaxy.n_HII.d3[ix][iy][iz].s[0])
		   *(fragment_p+galdef.He_H_ratio*fragment_He) *1.0e-27;
       }  //  p
   }  //  A!=0
     
   if(galdef.verbose>=2)
     {
       ost.str("");
       ost<< "fragmentation for species "<<particle.name;
       INFO(ost.str());
       particle.fragment.print();
     }

// ASSIGNING MOMENTUM LOSS RATE

   INFO("======== assigning momentum loss rate ======== ");

   if(galdef.n_spatial_dimensions==2)
   {
      for(int ir=0; ir<particle.n_rgrid; ir++)
      {
         for(int iz=0; iz<particle.n_zgrid; iz++)
         {
            for(int ip=0; ip<particle.n_pgrid; ip++)
            {
	      double aion(0.),coul(0.);                                      // NUCLEONS

               if(particle.A!=0) { 
                  //Enable turning off Coulomb and ionization losses individually
                  nucleon_loss(particle.Z,particle.mass,particle.Etot[ip],
                        galaxy.n_HI .d2[ir] [iz].s[0] +2*galaxy.n_H2.d2[ir] [iz].s[0],
                        galaxy.n_HII.d2[ir] [iz].s[0], galdef.He_H_ratio, 
                        &aion,  &coul); // energy eV s-1 -> momentum MeV s-1

                   particle.dpdt.d2[ir] [iz].s[ip] = 0;

                   if (galdef.ionization_losses) 
                      particle.dpdt.d2[ir] [iz].s[ip] += aion;

                   if (galdef.coulomb_losses) 
                      particle.dpdt.d2[ir] [iz].s[ip] += coul;

                   particle.dpdt.d2[ir] [iz].s[ip] /= particle.beta[ip]*1.0e6;
               }
      
               if(particle.A==0) {
		 
		 auto uevcm3(0.), bevcm3(0.), brem1(0.),brem2(0.),sync(0.),cmptn(0.);   // ELECTRONS
		 // IMOS200008016
		 // test of electron propagation vs analytical calculations (only IC losses) IMOS20061030
		 if(abs(galdef.DM_int0)==99)
		   {
		     particle.dpdt.  d2[ir] [iz].s[ip]= 
		       electron_loss( particle.Etot[ip], 0., 0., galdef.He_H_ratio, galdef.DM_double7, 0.,
				      &aion, &coul,&brem1,&brem2,&sync,&cmptn) *1.0e-6; // energy eV s-1 -> momentum MeV s-1
		     //32./9.*Pi*pow(Rele/Mele*1.e3,2)*C*galdef.DM_double7*pow(particle.Ekin[ip],2)*1.e-12;
		     
		     continue;
		   }
		 // end of the test area
		 
		 bevcm3=pow(galaxy.B_field.d2[ir][iz].s[0]*10000.,2)/8./Pi *erg_to_eV;// mag. energy density eV cm-3

                 //Enable to turn off various losses individually.
                 particle.  dpdt.d2[ir][iz].s[ip]= 0;

                 electron_loss(particle.Etot[ip],
                       galaxy.n_HI .d2[ir][iz].s[0] +2*galaxy.n_H2.d2[ir][iz].s[0],
                       galaxy.n_HII.d2[ir][iz].s[0],galdef.He_H_ratio, uevcm3, bevcm3,
                       &aion,&coul,&brem1,&brem2,&sync,&cmptn);

                 if (galdef.ionization_losses) 
                    particle.dpdt.d2[ir][iz].s[ip] += aion;

                 if (galdef.coulomb_losses) 
                    particle.dpdt.d2[ir][iz].s[ip] += coul;

                 if (galdef.bremss_losses)
                    particle.dpdt.d2[ir][iz].s[ip] += brem1 + brem2;

                 // This is added later using the e_KN_loss routine
                 //if (galdef.IC_losses)
                 //   particle.dpdt.d2[ir][iz].s[ip] += cmptn;

                 if (galdef.sync_losses)
                    particle.dpdt.d2[ir][iz].s[ip] += sync;

                 particle.  dpdt.d2[ir][iz].s[ip] /= particle.beta[ip]*1.0e6; // energy eV s-1 -> momentum MeV s-1
               }  //  A==0
	       // cout<<" dpdt="<<particle.dpdt.d2[ix]    [iz].s[ip]<<" aion="<<aion<<endl;
            }  //  p
         }  //  z
      }  //  r
   }
   
   if(galdef.n_spatial_dimensions==3)
   {
#pragma omp parallel for schedule(static) default(shared)
      for(int ix=0; ix<particle.n_xgrid; ix++)
      {
         for(int iy=0; iy<particle.n_ygrid; iy++)
         {
            for(int iz=0; iz<particle.n_zgrid; iz++)
            {
               for(int ip=0; ip<particle.n_pgrid; ip++)
               {
		 auto aion(0.),coul(0.);                                        // NUCLEONS

                  if(particle.A!=0) { 
                     //Enable turning off Coulomb and ionization losses individually
                     nucleon_loss(particle.Z,particle.mass,particle.Etot[ip],
                           galaxy.n_HI .d3[ix][iy][iz].s[0] +2*galaxy.n_H2.d3[ix][iy][iz].s[0],
                           galaxy.n_HII.d3[ix][iy][iz].s[0],galdef.He_H_ratio, 
                           &aion,&coul);

                     particle.dpdt.d3[ix][iy][iz].s[ip] = 0;

                     if (galdef.ionization_losses) 
                        particle.dpdt.d3[ix][iy][iz].s[ip] += aion;

                     if (galdef.coulomb_losses) 
                        particle.dpdt.d3[ix][iy][iz].s[ip] += coul;

                     particle.dpdt.d3[ix][iy][iz].s[ip] /= particle.beta[ip]*1.0e6; // energy eV s-1 -> momentum MeV s-1
                  }

                  if(particle.A==0)
                  {
		    auto uevcm3(0.), bevcm3(0.),brem1(0.),brem2(0.),sync(0.),cmptn(0.);   // ELECTRONS
		    // IMOS200008016
		    bevcm3=pow(galaxy.B_field.d3[ix][iy][iz].s[0]*10000.,2)/8./Pi *erg_to_eV;// mag. energy density eV cm-3

                    //Enable to turn off various losses individually.
		    particle.  dpdt.d3[ix][iy][iz].s[ip]= 0;
                    electron_loss(particle.Etot[ip],
                          galaxy.n_HI .d3[ix][iy][iz].s[0] +2*galaxy.n_H2.d3[ix][iy][iz].s[0],
                          galaxy.n_HII.d3[ix][iy][iz].s[0],galdef.He_H_ratio, uevcm3, bevcm3,
                          &aion,&coul,&brem1,&brem2,&sync,&cmptn);

                     if (galdef.ionization_losses) 
                        particle.dpdt.d3[ix][iy][iz].s[ip] += aion;

                     if (galdef.coulomb_losses) 
                        particle.dpdt.d3[ix][iy][iz].s[ip] += coul;

                    if (galdef.bremss_losses)
                        particle.dpdt.d3[ix][iy][iz].s[ip] += brem1 + brem2;

                    // This is added later using the e_KN_loss routine
                    //if (galdef.IC_losses)
                    //    particle.dpdt.d3[ix][iy][iz].s[ip] += cmptn;

                    if (galdef.sync_losses)
                        particle.dpdt.d3[ix][iy][iz].s[ip] += sync;

		    particle.  dpdt.d3[ix][iy][iz].s[ip] /= particle.beta[ip]*1.0e6; // energy eV s-1 -> momentum MeV s-1
        	  }  //  A==0
		  //  cout<<" dpdt="<<particle.dpdt.d3[ix][iy][iz].s[ip]<<" p="<<particle.p[ip] <<endl;
               }  //  p
            }  //  z
         }  //  y
      }  //  x
   }
   
   // IF ELECTRON ADD KLEIN_NISHINA LOSSES
   
   if(abs(galdef.DM_int0)!=99) 
      if(particle.A==0) 
         if (galdef.IC_losses)
            e_KN_loss(particle);  // MeV s-1 IMOS20061030
   
   if(galdef.verbose>=2)
     {
       ost.str("");
       ost<< "dpdt for species "<<particle.name;
       INFO(ost.str());
       particle.dpdt.print();
     }
   
   // ASSIGNING DECAY RATE
   
   particle.decay = 0.0;
   if(particle.t_half!=0.0)
     {
       INFO("======== assigning decay rate ======== ");

      if(galdef.n_spatial_dimensions==2)
      {
         for(int ir=0; ir<particle.n_rgrid; ir++)
         {
            for(int iz=0; iz<particle.n_zgrid; iz++)
            {
               for(int ip=0; ip<particle.n_pgrid; ip++)
                  particle.decay.d2[ir][iz].s[ip]=1.0/(particle.gamma[ip]*particle.t_half*year2sec/log(2.0));
            }  // z
         }  //  r
      }

      if(galdef.n_spatial_dimensions==3)
      {
#pragma omp parallel for schedule(static) default(shared)
         for(int ix=0; ix<particle.n_xgrid; ix++)
         {
            for(int iy=0; iy<particle.n_ygrid; iy++)
            {
               for(int iz=0; iz<particle.n_zgrid; iz++)
               {
                  for(int ip=0; ip<particle.n_pgrid; ip++)
                     particle.decay.d3[ix][iy][iz].s[ip]=1.0/(particle.gamma[ip]*particle.t_half*year2sec/log(2.0));
               }  //  z
            }  //  y
         }  //  x
      }

      if(galdef.verbose>=1)
      {
         ost.str("");
         ost<< "decay for species "<<particle.name;
	 INFO(ost.str());
         particle.decay.print();
      }
   }  //  t_half!=0.0

   if(galdef.verbose>=1)
   {
      ost.str("");
      ost<< "primary source function for species "<<particle.name;
      INFO(ost.str());
      particle.primary_source_function.print();
   }

//particle.print();
   ost.str("");
   ost<<"============== completed creation of transport arrays for "<<particle.name;
   INFO(ost.str());
   INFO("Exit");
   return stat;
}




