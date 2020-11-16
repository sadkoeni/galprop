
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_pi0_decay_emiss.cc *                      galprop package * 2/16/2005
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

/*
CR density  gcr.cr_density is in c/4pi * n(E) [ cm s^-1 sr^-1 cm^-3 MeV^-1]

emissivity (cm^-3 s^-1 sr^-1 MeV^-1)=
(c/4pi)*integral[sigma{Egamma,p_beam }  ) n(E)E dlog(E)]

pp_meson has Egamma in GeV, beam momentum in GeV
The particle spectra are assumed to be on equal kinetic energy per nucleon grids
which is the standard for galprop.
BUT UNITS OF density/momentum = flux/(KE/nucleon)..... CHECK CAREFULLY, ALSO factor A
cross section from pp_meson in barns GeV^-1
factor= 1.0e-24* *1.0e-3 log(Ekin_factor)
*/

#include <cassert>
#include <string>
#include <cstring>
#include <sstream>
#include <valarray>
#include <stdexcept>

using namespace std;//AWS20050624

#include "galprop_classes.h"
#include "galprop_internal.h"

//#include <fort_interface.h>

#include <Processes_Interface.h>

#include <ErrorLogger.h>

int Galprop::gen_pi0_decay_emiss() {   // generate pi0-decay emissivity

  INFO("Entry");

  ostringstream dBuf;
  dBuf << "Generating pi0-decay emissivity for n_spatial_dimensions="<<gcr[0].n_spatial_dimensions;
  INFO(dBuf.str());

  double factor,cs[2][2], y, P=-1.,Ekin=1.,Etot,beta,gamma,rig, yp,fp; // IMOS20061114 removed cs_p_HI etc.
  int  NA[2]={1,4}, iprotons, iHelium, i, key,test=0;                   // IMOS20061114 NA-array
  Distribution protons;

  // identify CR protons
  if (2 == galdef.n_spatial_dimensions) 
    protons.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
  
  if(3 == galdef.n_spatial_dimensions) 
    protons.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
  
  protons=0.;
  
  for(i=0, iprotons=-1; i<n_species; i++)  
    if(101==100*gcr[i].Z+gcr[i].A)
      {
	iprotons=i;
	protons+=gcr[iprotons].cr_density;
        std::ostringstream oss;
	oss<<"  CR protons found as species #"<<iprotons;
        INFO(oss.str());
      }
  if(iprotons==-1) { WARNING("  CR protons not found!"); protons.delete_array(); return 1; }
  
  // identify CR Helium
  for(i=0, iHelium =-1; i<n_species; i++) if(204==100*gcr[i].Z+gcr[i].A) iHelium =i;
  if(iHelium ==-1) WARNING("  CR Helium  not found!");
  else {
     std::ostringstream oss;
     oss<<"  CR Helium  found as species #"<<iHelium;
     INFO(oss.str());
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////

  if( galdef.pi0_decay==1 ||galdef.pi0_decay==2||galdef.pi0_decay==3) // AWS20120416
  {


  int key1=0; // pi0-decay

  valarray<double> p_gamma_matrix;
  
  valarray<double> param_a(0., 9);
  valarray<double> param_b(0., 8);
  valarray<double> param_c(0., 5);
  valarray<double> param_d(0., 5);

  if(galdef.pi0_decay==3)// fill in the p-gamma matrix: Kamae et al. 2006, ApJ 647, 692 // IMOS20061114
    {
      p_gamma_matrix.resize(gcr[iprotons].n_pgrid*galaxy.n_E_gammagrid);
      p_gamma_matrix = 0;
			    // = new float[gcr[iprotons].n_pgrid*galaxy.n_E_gammagrid];
      for(int ip=0; ip<gcr[iprotons].n_pgrid; ip++)
	{

	  param_a = 0;
	  param_b = 0;
	  param_c = 0;
	  param_d = 0;

	  param_a = processes::kamae_gamma_param_nd(gcr[iprotons].Ekin[ip]*1.e-3);//, param_a);
	  param_b = processes::kamae_gamma_param_diff(gcr[iprotons].Ekin[ip]*1.e-3);//, param_b);
	  param_c = processes::kamae_gamma_param_delta(gcr[iprotons].Ekin[ip]*1.e-3);//, param_c);
	  param_d = processes::kamae_gamma_param_res(gcr[iprotons].Ekin[ip]*1.e-3);//, param_d);
	  
	  for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++) {

	    const int index = ip + iEgamma*gcr[iprotons].n_pgrid;

	    p_gamma_matrix[index] = processes::kamae(galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip]*1.e-3, NA[0], NA[0], key1, param_a, param_b, param_c, param_d);

	    //*(p_gamma_matrix+ip+iEgamma*gcr[iprotons].n_pgrid) =
	    //  kamae(galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip]*1.e-3, NA[0], NA[0], key1, param_a, param_b, param_c, param_d);

	    //cout << ip << " " << gcr[iprotons].Ekin[ip] << " " << iEgamma << " " << galaxy.E_gamma[iEgamma] << " " << index << " " << p_gamma_matrix[index] << endl;

	  }

	}
    }
  
  //exit(0);

  galaxy.pi0_decay_emiss=0.;
  
  //Prepare the loop for more than the two elements by looping over the beam.

  for (int ibeam = 0; ibeam < 2; ++ibeam) {

     //Set up pointers to the right values
     Distribution *crDensity;
     Particle *particle;

     //If helium isn't found
     if (ibeam == 1 && iHelium == -1)
        continue;

     switch (ibeam) {
        case 0:
           crDensity = &protons;
           particle = &gcr[iprotons];
           break;
        case 1:
           crDensity = &gcr[iHelium].cr_density;
           particle = &gcr[iHelium];
           break;
     } 

#pragma omp parallel for schedule(dynamic) default(shared) private(cs,i,key,Ekin,P,yp,fp,y,Etot,beta,gamma,rig)
     for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
     {

        //Reset the cs values
        for(i=0;i<2;i++)  cs[0][i] = cs[1][i] = 0.;

        //Loop over momentum first because the cs are independent of position.
        for(int ip=0; ip<particle->n_pgrid; ip++)
        {

           //Reset the interpolation key
           key = 0;

           //Index for cross section array to calculate
           size_t ics = 0;

           //Need to adjust things if new integration mode
           if (galdef.integration_mode != 1) {

              ics = 1;

              //Copy upper boundary to lower
              for(i=0;i<2;i++) cs[0][i] = cs[1][i];
           }
           
           switch(galdef.pi0_decay) 
           { 
              default:
              case 1: // standard formalism (Moskalenko & Strong 1998, ApJ 493,694)
                 for (i=0;i<2;i++) // proton beam
                    cs[ics][i]  = processes::pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, particle->p[ip]*1.e-3, particle->A, NA[i], key1);
                 break;

              case 2: // Blattnig et al. 2000, PRD 62, #094030
                 for (i=0;i<2;i++) // proton beam
                    cs[ics][i]  = blattnig_gamma(galaxy.E_gamma[iEgamma]*1.e-3, particle->p[ip]*1.e-3, particle->A, NA[i], key1);
                 break;

              case 3: // Kamae et al. 2006, ApJ 647, 692
                 //Only for p-p 
                 if (ibeam == 0) {
                    // proton beam
                    const int index = ip+iEgamma*gcr[iprotons].n_pgrid;
                    cs[ics][0]  = p_gamma_matrix[index];//*(p_gamma_matrix+ip+iEgamma*gcr[iprotons].n_pgrid); // Kamae et al. 2006, ApJ 647, 692
                 } else {
                    cs[ics][0] = processes::pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, particle->p[ip]*1.e-3, particle->A, NA[0], key1);
                 }
                 cs[ics][1] = processes::pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, particle->p[ip]*1.e-3, particle->A, NA[1], key1);
                 break;
           }

           //Simply do nothing if both cross sections are 0
           if (cs[ics][0] == 0 && cs[ics][1] == 0)
              continue;

           //Skip the first step in new integration mode
           if (galdef.integration_mode != 1 && ip == 0)
              continue;

           //If new integration mode and the old cross section is 0 we find the minimum energy
           if (galdef.integration_mode != 1 && cs[0][0] == 0 && cs[0][1] == 0) {
              //Use brute force bisection method, 5 steps should be accurate enough
              //but continue until we have non-zero

              key = 1;
              size_t istep = 0;
              double emin = particle->Ekin[ip-1];
              double emax = particle->Ekin[ip];
              while (istep < 5 || ( cs[0][0] == 0 && cs[0][1] == 0 ) )
              {
                 //Calculate new point in momentum
                 Ekin = 0.5*(emin+emax);
                 P =-1.;
                 kinematic(particle->Z,particle->A,particle->mass,P,Ekin,Etot,beta,gamma,rig,test);

                 switch(galdef.pi0_decay) 
                 { 
                    default:
                    case 1: // standard formalism (Moskalenko & Strong 1998, ApJ 493,694)
                       for (i=0;i<2;i++) // proton beam
                          cs[0][i]  = processes::pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, P*1.e-3, particle->A, NA[i], key1);
                       break;

                    case 2: // Blattnig et al. 2000, PRD 62, #094030
                       for (i=0;i<2;i++) // proton beam
                          cs[0][i]  = blattnig_gamma(galaxy.E_gamma[iEgamma]*1.e-3, P*1.e-3, particle->A, NA[i], key1);
                       break;

                    case 3: // Kamae et al. 2006, ApJ 647, 692
                       //Only for p-p 
                       if (ibeam == 0) {
                          // proton beam
#pragma omp critical 
                          {

                             param_a = 0;
                             param_b = 0;
                             param_c = 0;
                             param_d = 0;

                             param_a = processes::kamae_gamma_param_nd(Ekin*1.e-3);//, param_a);
                             param_b = processes::kamae_gamma_param_diff(Ekin*1.e-3);//, param_b);
                             param_c = processes::kamae_gamma_param_delta(Ekin*1.e-3);//, param_c);
                             param_d = processes::kamae_gamma_param_res(Ekin*1.e-3);//, param_d);

                             cs[0][0] = processes::kamae(galaxy.E_gamma[iEgamma]*1.e-3, P*1.e-3, NA[0], NA[0], key1, param_a, param_b, param_c, param_d);
                          }
                       } else {
                          cs[0][0] = processes::pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, P*1.e-3, particle->A, NA[0], key1);
                       }
                       cs[0][1] = processes::pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, P*1.e-3, particle->A, NA[1], key1);
                       break;
                 }
                 
                 ++istep;

                 if (cs[0][0] == 0 && cs[0][1] == 0) {
                    emin = Ekin;
                 } else {
                    emax = Ekin;
                 }

                 //Emergency break
                 if (istep > 100) {
                    ERROR("Cannot find lower energy point");
                    std::ostringstream oss;
                    oss << "Cross sections: "<<cs[0][0]<<", "<<cs[0][1]<<", "<<cs[1][0]<<", "<<cs[1][1];
                    ERROR(oss.str());
                    oss.str("");
                    oss << "Gamma energy: "<<galaxy.E_gamma[iEgamma]<<" and beam energy "<<Ekin<<", "<<particle->Ekin[ip-1]<<", "<<particle->Ekin[ip];
                    ERROR(oss.str());
                    break;
                 }
              }
              if (istep > 100)
                 continue;
           }

           //Make assertion regarding the validity of the css
           for (size_t ii=0; ii < 2; ++ii) {
              bool stop = false;
              if (cs[ics][ii] < 0) {
                 ERROR("Cross section < 0");
                 stop = true;
              }
              if (std::isnan(cs[ics][ii]) || std::isinf(cs[ics][ii])) {
                 ERROR("Cross section is NaN");
                 stop = true;
              }
              if (stop) {
                 std::ostringstream oss;
                 oss<<"Error for CS at indices "<<ics<<" and "<<ii;
                 ERROR(oss.str());
                 oss.str("");
                 oss<<"Gamma ray energy is "<<galaxy.E_gamma[iEgamma]<<" and proton momentum "<<gcr[iprotons].Ekin[ip]<<" - "<<gcr[iprotons].Ekin[ip+1];
                 ERROR(oss.str());
                 oss.str("");
                 throw(std::runtime_error("Nans in routine"));
              }
           }

           //Perform the actual integration
           if(galaxy.n_spatial_dimensions==2)
           {
              for(int ir=0; ir<gcr[iprotons].n_rgrid-1; ir++)
              {
                 for(int iz=1; iz<gcr[iprotons].n_zgrid-1; iz++)
                 { 
                    //************************************** TEST
                    if(galdef.verbose==-217)
                    {
                       if(galdef.integration_mode==1 && ip>0) Ekin=0.5*(particle->Ekin[ip-1]+particle->Ekin[ip  ]);
                       for(i=0;i<2;i++)  cs[0][i] = cs[1][i] = 1.;
                       crDensity->d2[ir][iz].s[ip-1] = pow(particle->Ekin[ip-1],-3);
                       crDensity->d2[ir][iz].s[ip  ] = pow(particle->Ekin[ip  ],-3);
                       galdef.He_H_ratio=0.;
                    }
                    //***************************************/
		  
                    if( galdef.integration_mode==1 ) // ### old integration ###
                    {
                       galaxy.pi0_decay_emiss.d2[ir][iz].s[iEgamma]+=       
                        (cs[0][0] +cs[0][1] *galdef.He_H_ratio) *crDensity->d2[ir][iz].s[ip] *particle->Ekin[ip]*particle->A;
		      
                       if(galdef.verbose==-701)//selectable debug AWS20040303
                       {
                          if(iEgamma==10)
                          {
                             cout<<"ir iz E_gamma pi0_decay_emiss "
                                <<ir<<" "<<iz<<" "<<" "<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.pi0_decay_emiss.d2[ir][iz].s[iEgamma];
                             cout<<" gcr[iprotons].Ekin[ip]="<<gcr[iprotons].Ekin[ip];
                             cout<<" cs_p_HI cs_p_He cs_He_HI+cs_He_He "<<cs[0][0]<<"  "<<cs[0][1]<<"  "<<cs[1][0]<<"  "<<cs[1][1];
                             cout<<" protons= "<<protons.d2[ir][iz].s[ip]<<" Helium= "<<gcr[iHelium ].cr_density.d2[ir][iz].s[ip]<<endl;
                          }
                       }
                       continue;
                    } 
                    else 
                    {
                       // ### new integration - analytical ###
                       if(key==1)  // case: lower integration limit falls between the grid points
                       { 
                          // interpolate beam spectrum (power-law) or linear
                          if (crDensity->d2[ir][iz].s[ip-1] > 0 && crDensity->d2[ir][iz].s[ip] > 0) 
                          {
                             yp = log( crDensity->d2[ir][iz].s[ip-1] / crDensity->d2[ir][iz].s[ip] )     // yp= power-law index
                                / log( particle->Ekin[i-1] / particle->Ekin[ip] );
                             fp = crDensity->d2[ir][iz].s[ip-1]*pow(particle->Ekin[ip-1],-yp);// normalization
                             fp*= pow(Ekin,yp);                                               // proton flux @ Ekin
                          }
                          else
                          {
                             fp = (Ekin-particle->Ekin[ip-1])*
                                (crDensity->d2[ir][iz].s[ip]-crDensity->d2[ir][iz].s[ip-1]) / 
                                (particle->Ekin[ip] - particle->Ekin[i-1]) + crDensity->d2[ir][iz].s[ip-1];
                          }
			  
                       } 
                       else
                       {  // case: integration between the grid points
                          Ekin=particle->Ekin[ip-1]; // Ekin -kin.energy/nucleon, same for p & He
                          fp = crDensity->d2[ir][iz].s[ip-1];
                       }
		      
                       // protons: derive y (= power-law index of the total expression)
                       const double fmin = (cs[0][0] +cs[0][1]*galdef.He_H_ratio) *fp;
                       const double fmax = (cs[1][0] +cs[1][1]*galdef.He_H_ratio) *crDensity->d2[ir][iz].s[ip];

                        if (fmin <= 0 || fmax <= 0) {
                           WARNING("Values <= 0 in analytic integration");
                           std::ostringstream oss;
                           oss << "Cross sections: "<<cs[0][0]<<", "<<cs[0][1]<<", "<<cs[1][0]<<", "<<cs[1][1];
                           WARNING(oss.str());
                           oss.str("");
                           oss << "Beam is "<<ibeam;
                           WARNING(oss.str());
                        } else {
                           y = log(fmin/fmax) / log(Ekin/particle->Ekin[ip]);
                      
                           //Check for nans
                           if (std::isnan(y) || std::isinf(y)) {
                              ERROR("NaNs in analytic integration");
                              std::ostringstream oss;
                              oss << "Momentum index "<<ip<<" with energy "<<Ekin<<" - "<<particle->Ekin[ip];
                              ERROR(oss.str());
                              oss.str("");
                              oss << "beam density "<<fp<<" - "<<crDensity->d2[ir][iz].s[ip];
                              ERROR(oss.str());
                              oss.str("");
                              oss << "Beam is "<<ibeam;
                              ERROR(oss.str());
                              throw(std::runtime_error("Nans in routine"));
                           }

                           // integrate analytically
                           galaxy.pi0_decay_emiss.d2[ir][iz].s[iEgamma] +=(cs[1][0] +cs[1][1]*galdef.He_H_ratio) 
                              *crDensity->d2[ir][iz].s[ip]*particle->A
                              *particle->Ekin[ip]/(y+1.)
                              *(1.-pow(Ekin/ particle->Ekin[ip],y+1.));

                           if (std::isnan(galaxy.pi0_decay_emiss.d2[ir][iz].s[iEgamma]) || std::isinf(galaxy.pi0_decay_emiss.d2[ir][iz].s[iEgamma])) {
                              ERROR("NaNs in analytic integration");
                              std::ostringstream oss;
                              oss << "Momentum index "<<ip<<" with energy "<<Ekin<<" - "<<particle->Ekin[ip];
                              ERROR(oss.str());
                              oss.str("");
                              oss << "power law index "<<y;
                              ERROR(oss.str());
                              oss.str("");
                              oss << "Beam is "<<ibeam;
                              ERROR(oss.str());
                              throw(std::runtime_error("Nans in routine"));
                           }
                        }
		      
		      
                    }
		  
                 }//iz 
              }//ir 
           }//particle.n_spatial_dimensions==2
	  
	  
           //************************************** TEST
           if(galdef.verbose==-217)
              if(iEgamma==galaxy.n_E_gammagrid/2)
                 cout<<" Ekin= "<<Ekin<<"  "<<gcr[iprotons].Ekin[ip]<<" integral= "<<galaxy.pi0_decay_emiss.d2[gcr[iprotons].n_rgrid/2][gcr[iprotons].n_zgrid/2].s[iEgamma]<<"  >>>"<<pow(Ekin,-2)/2.<<"  " <<pow(Ekin,-2)*5./2.<<endl;
           //***************************************/
	  
           if(galaxy.n_spatial_dimensions==3)
           {
              for(int ix=1; ix<gcr[iprotons].n_xgrid-1; ix++)
              {
                 for(int iy=1; iy<gcr[iprotons].n_ygrid-1; iy++)
                 {
                    for(int iz=1; iz<gcr[iprotons].n_zgrid-1; iz++)
                    {
                       if( galdef.integration_mode==1 ) // ### old integration ###
                       {
                          galaxy.pi0_decay_emiss.d3[ix][iy][iz].s[iEgamma]+=       
                             (cs[0][0] +cs[0][1] *galdef.He_H_ratio) *crDensity->d3[ix][iy][iz].s[ip] *particle->Ekin[ip] *particle->A;

                          if(galdef.verbose==-701)//selectable debug AWS20040303
                          {
                             if(iEgamma==10)
                             {
                                cout<<"ix iy iz E_gamma pi0_decay_emiss "
                                   <<ix<<" "<<iy<<" "<<iz<<" "<<" "<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.pi0_decay_emiss.d3[ix][iy][iz].s[iEgamma];
                                cout<<" gcr[iprotons].Ekin[ip]="<<gcr[iprotons].Ekin[ip];
                                cout<<" cs_p_HI cs_p_He cs_He_HI+cs_He_He "<<cs[0][0]<<"  "<<cs[0][1]<<"  "<<cs[1][0]<<"  "<<cs[1][1];
                                cout<<" protons= "<<protons.d3[ix][iy][iz].s[ip]<<" Helium= "<<gcr[iHelium ].cr_density.d3[ix][iy][iz].s[ip]<<endl;
                             }
                          }
                          continue;
                       } 
                       else 
                       {
                          // ### new integration - analytical ###
                          if(key==1)  // case: lower integration limit falls between the grid points
                          { 
                             // interpolate beam spectrum (power-law) or linear
                             if (crDensity->d3[ix][iy][iz].s[ip-1] > 0 && crDensity->d3[ix][iy][iz].s[ip] > 0) 
                             {
                                yp = log( crDensity->d3[ix][iy][iz].s[ip-1] / crDensity->d3[ix][iy][iz].s[ip] )     // yp= power-law index
                                   / log( particle->Ekin[i-1] / particle->Ekin[ip] );
                                fp = crDensity->d3[ix][iy][iz].s[ip-1]*pow(particle->Ekin[ip-1],-yp);// normalization
                                fp*= pow(Ekin,yp);                                               // proton flux @ Ekin
                             }
                             else
                             {
                                fp = (Ekin-particle->Ekin[ip-1])*
                                   (crDensity->d3[ix][iy][iz].s[ip]-crDensity->d3[ix][iy][iz].s[ip-1]) / 
                                   (particle->Ekin[ip] - particle->Ekin[i-1]) + crDensity->d3[ix][iy][iz].s[ip-1];
                             }

                          } 
                          else
                          {  // case: integration between the grid points
                             Ekin=particle->Ekin[ip-1]; // Ekin -kin.energy/nucleon, same for p & He
                             fp = crDensity->d3[ix][iy][iz].s[ip-1];
                          }

                          // protons: derive y (= power-law index of the total expression)
                          const double fmin = (cs[0][0] +cs[0][1]*galdef.He_H_ratio) *fp;
                          const double fmax = (cs[1][0] +cs[1][1]*galdef.He_H_ratio) *crDensity->d3[ix][iy][iz].s[ip];

                          if (fmin <= 0 || fmax <= 0) {
                             WARNING("Values <= 0 in analytic integration");
                             std::ostringstream oss;
                             oss << "Cross sections: "<<cs[0][0]<<", "<<cs[0][1]<<", "<<cs[1][0]<<", "<<cs[1][1];
                             WARNING(oss.str());
                             oss.str("");
                             oss << "Beam is "<<ibeam;
                             WARNING(oss.str());
                          } else {
                             y = log(fmin/fmax) / log(Ekin/particle->Ekin[ip]);

                             //Check for nans
                             if (std::isnan(y) || std::isinf(y)) {
                                ERROR("NaNs in analytic integration");
                                std::ostringstream oss;
                                oss << "Momentum index "<<ip<<" with energy "<<Ekin<<" - "<<particle->Ekin[ip];
                                ERROR(oss.str());
                                oss.str("");
                                oss << "beam density "<<fp<<" - "<<crDensity->d3[ix][iy][iz].s[ip];
                                ERROR(oss.str());
                                oss.str("");
                                oss << "Beam is "<<ibeam;
                                ERROR(oss.str());
                                throw(std::runtime_error("Nans in routine"));
                             }

                             // integrate analytically
                             galaxy.pi0_decay_emiss.d3[ix][iy][iz].s[iEgamma] +=(cs[1][0] +cs[1][1]*galdef.He_H_ratio) 
                                *crDensity->d3[ix][iy][iz].s[ip] * particle->A
                                *particle->Ekin[ip]/(y+1.)
                                *(1.-pow(Ekin/ particle->Ekin[ip],y+1.));

                             if (std::isnan(galaxy.pi0_decay_emiss.d3[ix][iy][iz].s[iEgamma]) || std::isinf(galaxy.pi0_decay_emiss.d3[ix][iy][iz].s[iEgamma])) {
                                ERROR("NaNs in analytic integration");
                                std::ostringstream oss;
                                oss << "Momentum index "<<ip<<" with energy "<<Ekin<<" - "<<particle->Ekin[ip];
                                ERROR(oss.str());
                                oss.str("");
                                oss << "power law index "<<y;
                                ERROR(oss.str());
                                oss.str("");
                                oss << "Beam is "<<ibeam;
                                ERROR(oss.str());
                                throw(std::runtime_error("Nans in routine"));
                             }
                          }

                       }

                    }//iz
                 }//iy
              }//ix
           } //3D
	  
           //cout<<"gen_pi0_decay_emiss: iEgamma ip "<<iEgamma<<" "<<ip<<endl;
        } //ip

        //************************************** TEST
        if(galdef.verbose==-217)
           if(iEgamma==galaxy.n_E_gammagrid/2) 
           {
              factor= log(galdef.Ekin_factor);
              galaxy.pi0_decay_emiss        *= factor;
              if(galdef.integration_mode==1) cout<<"* integral= "<<galaxy.pi0_decay_emiss.d2[gcr[iprotons].n_rgrid/2][gcr[iprotons].n_zgrid/2].s[iEgamma]<<endl;
              //exit(1);
           }
        //***************************************/
     } //iEgamma
  }
     

  
  factor= (galdef.integration_mode==1) ? 1.0e-24*1.0e-3* log(galdef.Ekin_factor) : 1.0e-24*1.0e-3;
  galaxy.pi0_decay_emiss *= factor;
  //if(galdef.pi0_decay==3) delete[] p_gamma_matrix; // IMOS20061114 //Gulli20070810
  protons.delete_array();                        // IMOS20030217
  

  }//  galdef.pi0_decay==1 ||galdef.pi0_decay==2||galdef.pi0_decay==3


  if(galdef.verbose==-702) // selectable debug
 { cout<<"   pi0-decay emissivity "<<endl; galaxy.pi0_decay_emiss.print(); }
  //cout<<" <<<< gen_pi0_decay_emiss"<<endl;

  if(galdef.verbose==-478) exit(0); 


  INFO("Exit");

  return 0;

}
