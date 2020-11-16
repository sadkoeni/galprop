
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_secondary_antiproton_source.cc *          galprop package * 2001/05/11
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// The routine to calculate the antiproton source function. 
//
// CR density gcr.cr_density is in c/4pi * n(E) [cm s^-1 sr^-1 cm^-3 MeV^-1]
//
// The routine ANTIPROTON written in FORTRAN-77 is designed to calculate
// the antiproton (+antineutron) production spectrum vs. momentum (barn/GeV). 
// Antiproton momentum and nucleus momentum (GeV) per nucleon are used as input
// parameters as well as beam and target nuclei atomic numbers.
//
// The antiproton source function [cm^-2 s^-2 sr^-1 MeV^-1] as used in galprop is 
// defined as following (c/4pi * q)  [q = cm^-3 s^-1 MeV^-1]:
//                ___      ___  
//         c      \        \    /                   c   d sigma_ij(p,p')
// q(p) * --- = c /__  n_i /__  \ dp' beta n_j(p') ---  --------------- ,
//        4pi    i=H,He     j   /                  4pi        dp  
// 
// where n_i is the gas density, d sigma_ij(p,p')/dp is
// the production cross section, n_j(p') is the CR species density, 
// and p' is the total momentum of a nucleus.
// Substitution of dp' with d(log Ekin) gives:
//                ___                          ___  
//       c        \        /                   \               c  d sigma_ij(p,Ekin)
// q(p)*--- = c A /__ n_i  \ d(log Ekin)  Ekin /__  n_j(Ekin) --- -----------------
//      4pi      i=H,He    /                    j             4pi       dp  
//                         ___     ___      ___
//                         \       \        \              c   d sigma_ij(p,Ekin)
//      = c A /\(log Ekin) /__ n_i /__ Ekin /__ n_j(Ekin) ---  -----------------,
//                        i=H,He   Ekin      j            4pi        dp             
// 
// where /\=Delta, and we used dp'=1/beta A Ekin d(log Ekin).
// 
// To transfer to units cm^2/MeV we need a factor= 1.0e-24 *1.0e-3.
// Ref.: Moskalenko I.V. et al. 2002, ApJ 565, 280
//=="====!===="====!===="====!===="====!===="====!===="====!===="====!===="====!

using namespace std;//AWS20050624

#include"galprop_classes.h"
#include"galprop_internal.h"

//#include <fort_interface.h>

#include <Processes_Interface.h>

#include <cstring>
#include <sstream>

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galprop::gen_secondary_antiproton_source(Particle &particle)
{
   INFO("ENTRY");

   ostringstream buf;
   buf<<"Generating "<<particle.name<<" source function for n_spatial_dimensions="
       <<gcr[0].n_spatial_dimensions;
   INFO(buf.str());

   if ("secondary_antiprotons" != particle.name) 
   {  
      buf.str("");
      buf<<"invalid particle "<<particle.name; 
      WARNING(buf.str());
      return 2; 
   }

   if (galdef.secondary_antiprotons < 1 || galdef.secondary_antiprotons > 3) 
   {
      // Cross section not known
      ERROR("option for secondary_antiprotons is not known");
      ERROR("Available options: 1, 2, and 3");
      return 3;
   }

   int stat=0, iprotons=-1, iHelium =-1,  Z1, A1, Z2, A2;

   //Store the cross sections and indices in vectors
   //An index of -1 means use protons distribution
   //The first value in all vectors should always exist and the corresponding index should be to protons.
   std::vector<double> cs_HI, cs_He;
   std::vector<int> ind_HI, ind_He;
   cs_HI.reserve(n_species);
   cs_He.reserve(n_species);
   ind_HI.reserve(n_species);
   ind_He.reserve(n_species);

   Distribution protons;                 // IMOS20000606.6

// identify CR protons                   // IMOS20000606.7
   if(galdef.n_spatial_dimensions==2) protons.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
   if(galdef.n_spatial_dimensions==3) protons.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
   protons=0.;
   for(int i=0; i<n_species; i++)  
      if(101==100*gcr[i].Z+gcr[i].A)
      {
         iprotons=i;
 	 protons+=gcr[iprotons].cr_density;
         buf.str("");
         buf<<"  CR protons found as species #"<<iprotons;
         INFO(buf.str());
      }
   if(iprotons==-1) { ERROR("CR protons not found!"); return 1; }
 
// identify CR Helium
   for(int i=0; i<n_species; i++) if(204 == 100*gcr[i].Z+gcr[i].A) iHelium =i;
   if(iHelium ==-1) { ERROR("CR Helium  not found!"); return 1; }
   else if(galdef.verbose>=1) {
      buf.str("");
      buf<<"  CR Helium  found as species #"<<iHelium;
      INFO(buf.str());
   }

//Gulli20070821 
#pragma omp parallel for schedule(dynamic) default(shared) firstprivate(Z1,Z2,A1,A2,cs_HI,cs_He,ind_HI,ind_He)
   for(int ip_sec=0; ip_sec<particle.n_pgrid; ip_sec++)
   {
      for(int ip=0; ip<gcr[iprotons].n_pgrid; ip++)
      {
         cs_HI.clear();
         cs_He.clear();
         ind_HI.clear();
         ind_He.clear();
         if (galdef.secondary_antiprotons == 1 || galdef.secondary_antiprotons == 2) {

            Z1=gcr[iprotons].Z;  A1=gcr[iprotons].A;  Z2=1;  A2=1;    // beam+target: p+HI
            cs_HI.push_back(processes::antiproton_cc(galdef.total_cross_section,particle.p[ip_sec]*1.e-3, gcr[iprotons].p[ip]*1.e-3, Z1,A1,Z2,A2));
            ind_HI.push_back(-1); //Use the protons

            // secondary_antiprotons =1 uses scaling to calc.; =2 uses factors by Simon et al. 1998
            if(galdef.secondary_antiprotons == 2)                                                     // IMOS20000802.2
            {
               cs_HI.back() *= 0.12/pow(particle.Ekin[ip_sec]/1000,1.67)+1.78;
               cs_He.push_back(0);
               ind_He.push_back(-1);
            }
            else
            {
               Z1=gcr[iprotons].Z;  A1=gcr[iprotons].A;  Z2=2;  A2=4;    // beam+target: p+He
               cs_He.push_back(processes::antiproton_cc(galdef.total_cross_section,particle.p[ip_sec]*1.e-3, gcr[iprotons].p[ip]*1.e-3, Z1,A1,Z2,A2)); 
               ind_He.push_back(-1);

               Z1=gcr[iHelium ].Z;  A1=gcr[iHelium ].A;  Z2=1;  A2=1;    // beam+target: He+HI
               cs_HI.push_back(processes::antiproton_cc(galdef.total_cross_section,particle.p[ip_sec]*1.e-3, gcr[iHelium ].p[ip]*1.e-3, Z1,A1,Z2,A2));
               ind_HI.push_back(iHelium);

               Z1=gcr[iHelium ].Z;  A1=gcr[iHelium ].A;  Z2=2;  A2=4;    // beam+target: He+He
               cs_He.push_back(processes::antiproton_cc(galdef.total_cross_section,particle.p[ip_sec]*1.e-3, gcr[iHelium ].p[ip]*1.e-3, Z1,A1,Z2,A2)); 
               ind_He.push_back(iHelium);
            }
         } 
         else if (galdef.secondary_antiprotons == 3) 
         {
            // Kachelriess, Moskalenko & Ostapchenko, arXiv:1502.04158
            // Convert from Etpn d sigma/dEtpn to d sigma/dp and mbarn to barn.  Etpn is Etot/A for the anti-proton.  
            const double conv = particle.p[ip_sec]/particle.Etot[ip_sec]/particle.Etot[ip_sec];
            
            //Start by adding protons
            cs_HI.push_back(processes::apspec_cc(gcr[iprotons].Etot[ip]*1e-3/gcr[iprotons].A, particle.Etot[ip_sec]*1e-3, gcr[iprotons].A, 1)*conv);
            cs_He.push_back(processes::apspec_cc(gcr[iprotons].Etot[ip]*1e-3/gcr[iprotons].A, particle.Etot[ip_sec]*1e-3, gcr[iprotons].A, 4)*conv);
            ind_HI.push_back(-1);
            ind_He.push_back(-1);

            // Use all available nuclei up to iron
            for ( int j = 0; j < n_species; ++j ) {
               if ( (gcr[j].A > 1 && gcr[j].Z > 1 && gcr[j].A < 61) ) {
                  cs_HI.push_back(processes::apspec_cc(gcr[j].Etot[ip]*1e-3/gcr[j].A, particle.Etot[ip_sec]*1e-3, gcr[j].A, 1)*conv);
                  cs_He.push_back(processes::apspec_cc(gcr[j].Etot[ip]*1e-3/gcr[j].A, particle.Etot[ip_sec]*1e-3, gcr[j].A, 4)*conv);
                  ind_HI.push_back(j);
                  ind_He.push_back(j);
               }
            }
         }

         if(galaxy.n_spatial_dimensions==2)
         {
            for(int ir=0; ir<gcr[iprotons].n_rgrid; ir++)
            {
               for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
               {
                  const double gas = galaxy.n_HI.d2[ir][iz].s[0]+2.0*galaxy.n_H2.d2[ir][iz].s[0]+galaxy.n_HII.d2[ir][iz].s[0];

                  //First is always protons and always exists
                  particle.secondary_source_function.d2[ir][iz].s[ip_sec ] +=  
                     gas * (cs_HI[0] + galdef.He_H_ratio * cs_He[0]) * 
                     protons.d2[ir][iz].s[ip] * gcr[iprotons].Ekin[ip];

                  //Loop over the rest
                  for (size_t ics=1; ics < cs_HI.size(); ++ics) {
                     particle.secondary_source_function.d2[ir][iz].s[ip_sec ]+=  
                        gas * cs_HI[ics] * gcr[ind_HI[ics]].cr_density.d2[ir][iz].s[ip] * 
                        gcr[ind_HI[ics]].Ekin[ip] * gcr[ind_HI[ics]].A;
                  }
                  for (size_t ics=1; ics < cs_He.size(); ++ics) {
                     particle.secondary_source_function.d2[ir][iz].s[ip_sec ]+=  
                        gas * galdef.He_H_ratio * cs_He[ics] * gcr[ind_He[ics]].cr_density.d2[ir][iz].s[ip] * 
                        gcr[ind_He[ics]].Ekin[ip] * gcr[ind_HI[ics]].A;
                  }
                        
               }  //  iz
            }  //  ir
         }  //  particle.n_spatial_dimensions==2

         if(galaxy.n_spatial_dimensions==3)
         {
            for(int ix=0; ix<gcr[iprotons].n_xgrid; ix++)
            {
               for(int iy=0; iy<gcr[iprotons].n_ygrid; iy++)
               {
                  for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
                  {
                     const double gas = galaxy.n_HI.d3[ix][iy][iz].s[0]+2.0*galaxy.n_H2.d3[ix][iy][iz].s[0]+galaxy.n_HII.d3[ix][iy][iz].s[0];

                     //First is always protons and always exists
                     particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec ] +=  
                        gas * (cs_HI[0] + galdef.He_H_ratio * cs_He[0]) * 
                        protons.d3[ix][iy][iz].s[ip] * gcr[iprotons].Ekin[ip];

                     //Loop over the rest
                     for (size_t ics=1; ics < cs_HI.size(); ++ics) {
                        particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec ]+=  
                           gas * cs_HI[ics] * gcr[ind_HI[ics]].cr_density.d3[ix][iy][iz].s[ip] * 
                           gcr[ind_HI[ics]].Ekin[ip] * gcr[ind_HI[ics]].A;
                     }
                     for (size_t ics=1; ics < cs_He.size(); ++ics) {
                        particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec ]+=  
                           gas * galdef.He_H_ratio * cs_He[ics] * gcr[ind_He[ics]].cr_density.d3[ix][iy][iz].s[ip] * 
                           gcr[ind_He[ics]].Ekin[ip] * gcr[ind_HI[ics]].A;
                     }
                  }  //  iz
               }  //  iy
            }  //  ix
         }  //  particle.n_spatial_dimensions==3
      }  //  ip
   }  //  ip_sec


   const double factor=1.e-24 *1.e-3 *C *log(galdef.Ekin_factor); // transformation to cm2/MeV and constant factors
   particle.secondary_source_function *= factor;

   protons.delete_array();                  // IMOS20000606.10

   if(galdef.verbose>=2)
   {
      cout<<"   particle.secondary_source_function for "<<particle.name<<endl;
      particle.secondary_source_function.print();
    }
   INFO("Exit");
   return stat;
}
