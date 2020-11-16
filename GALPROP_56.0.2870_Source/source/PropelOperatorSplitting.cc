#include "PropelOperatorSplitting.h"
#include <ErrorLogger.h>
#include <sstream>
#include <constants.h>

//This has 512 bit alignement by default
#include <Malloc.h>

//The operator initializes the memory, assume that does not change
template <typename T>
PropelOperatorSplitting<T>::PropelOperatorSplitting( const Particle &particle, const Galdef &_galdef) :
   PropelBase(particle, _galdef),
   alpha1_p(nullptr), alpha1_z(nullptr), alpha1_r(nullptr), alpha1_x(nullptr), alpha1_y(nullptr),
   alpha2_p(nullptr), alpha2_z(nullptr), alpha2_r(nullptr), alpha2_x(nullptr), alpha2_y(nullptr),
   alpha3_p(nullptr), alpha3_z(nullptr), alpha3_r(nullptr), alpha3_x(nullptr), alpha3_y(nullptr)
{ 
   //Allocate the memory for the arrays
   if (2 == particle.n_spatial_dimensions) { // ==== 2D ====

      alpha1_p = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha1_r = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha1_z = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha2_p = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha2_r = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha2_z = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha3_p = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha3_r = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha3_z = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));

   } 
   else if (3 == particle.n_spatial_dimensions) { // ==== 3D ====

      alpha1_p = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha1_x = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha1_y = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha1_z = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha2_p = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha2_x = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha2_y = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha2_z = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha3_p = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha3_x = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha3_y = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      alpha3_z = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));

   }

}

template <typename T>
PropelOperatorSplitting<T>::~PropelOperatorSplitting()
{
   mlc::FreeAligned(alpha1_p);
   mlc::FreeAligned(alpha1_z);
   mlc::FreeAligned(alpha1_r);
   mlc::FreeAligned(alpha1_x);
   mlc::FreeAligned(alpha1_y);
   mlc::FreeAligned(alpha2_p);
   mlc::FreeAligned(alpha2_z);
   mlc::FreeAligned(alpha2_r);
   mlc::FreeAligned(alpha2_x);
   mlc::FreeAligned(alpha2_y);
   mlc::FreeAligned(alpha3_p);
   mlc::FreeAligned(alpha3_z);
   mlc::FreeAligned(alpha3_r);
   mlc::FreeAligned(alpha3_x);
   mlc::FreeAligned(alpha3_y);
}

template <typename T>
void PropelOperatorSplitting<T>::tridiag(
      const T* a,
      const T* b,
      const T* c,
      const T* r,
      T* u,
      T* gam,
      const size_t n)
{
   //Assume b[0] != 0, it should not be and the test will slow us down
   T bet = b[0];
   u[0] = r[0]/bet;
   for (size_t j = 1; j < n; ++j) {
      gam[j] = c[j-1]/bet;
      bet = b[j] - a[j]*gam[j];
      //Again remove the test for speed concern.  The result will be junk if there is no solution.
      u[j] = (r[j] - a[j]*u[j-1])/bet;
   }

   //This has to be int
   for (int j = n-2; j >=0; --j)
      u[j] = u[j] - gam[j+1]*u[j+1];

}

template <typename T>
void PropelOperatorSplitting<T>::vectorTridiag(
      const T* a,
      const T* b,
      const T* c,
      const T* r,
      T* u,
      T* gam,
      T* bet,
      const size_t n,
      const size_t stride)
{

#if _OPENMP >= 201307
#pragma omp simd aligned(b,bet,u,r:64)
#endif
   for (size_t k = 0; k < stride; ++k) {
      bet[k] = 1./b[k];
      u[k] = r[k]*bet[k];
   }

   for (size_t j = 1; j < n; ++j) {
      T* pgam = gam+j*stride;
      const T* pc = c+(j-1)*stride;
      const T* pa = a+j*stride;
      const T* pb = b+j*stride;
      const T* pr = r+j*stride;
      T* pu = u+j*stride;
      T* pu1 = u+(j-1)*stride;
#if _OPENMP >= 201307
#pragma omp simd aligned(pa,pb,pc,pr,pu,pu1,bet,pgam:64)
#endif
      for (size_t k = 0; k < stride; ++k) {
         pgam[k] = pc[k] * bet[k];
         bet[k] = 1./(pb[k] - pa[k] * pgam[k]);
         pu[k] = (pr[k] - pa[k] * pu1[k]) * bet[k];
      }
   }

   for (size_t j = n-1; j >0; --j) {
      T* pgam = gam+j*stride;
      T* pu = u+j*stride;
      T* pu1 = u+(j-1)*stride;
#if _OPENMP >= 201307
#pragma omp simd aligned(pu,pu1,pgam:64)
#endif
      for (size_t k = 0; k < stride; ++k) {
         pu1[k] -= pgam[k] * pu[k];
      }
   }

}

//Set up the alphas
template <typename T>
bool PropelOperatorSplitting<T>::initializeAlpha(Particle &particle)
{
   //Code mostly copied from propel
   if (0 == particle.primary_source_function.max() && 
         0 == particle.secondary_source_function.max()) {

      std::ostringstream buf;
      buf << particle.name<<": Zero primary and secondary source functions. No propagation will be done";
      WARNING(buf.str());
      return false;

   }

   const T factor = 1./(kpc2cm*kpc2cm);//pow(kpc2cm, -2.);  

   if (2 == particle.n_spatial_dimensions) { // ==== 2D ====

      //Create an indexer
      const IndexConverter2D index(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);

      // DIFFUSION term: method according to ApJ 509, 212
#pragma omp parallel for default(shared) schedule(static)
      for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
         for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
            T* a1z = alpha1_z+index(ir,iz,0);
            T* a2z = alpha2_z+index(ir,iz,0);
            T* a3z = alpha3_z+index(ir,iz,0);
            T* a1r = alpha1_r+index(ir,iz,0);
            T* a2r = alpha2_r+index(ir,iz,0);
            T* a3r = alpha3_r+index(ir,iz,0);
            T* a1p = alpha1_p+index(ir,iz,0);
            T* a2p = alpha2_p+index(ir,iz,0);
            T* a3p = alpha3_p+index(ir,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a1z, a2z, a3z, a2r, a1r, a3r, a1p, a2p, a3p:64)
#endif
            for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
               const T zvalue = particle.Dxx.d2[ir][iz].s[ip]*pow(particle.dz,-2.)*factor;
               a1z[ip] = zvalue;
               a2z[ip] = zvalue*2.;
               a3z[ip] = zvalue;

               const T rvalue = particle.Dxx.d2[ir][iz].s[ip]*pow(particle.dr,-2.)*factor;
               a2r[ip] = rvalue*2.;
               //Moved the 0 == ir condition outside the loop

               // use: Dxx/(Ri dR)*{R[i+1/2](U[i+1]-Ui)-R[i-1/2](Ui-U[i-1])}
               a1r[ip] = rvalue*(1. - T(particle.dr/2./particle.r[ir]));
               a3r[ip] = rvalue*(1. + T(particle.dr/2./particle.r[ir]));
               //Initialize alpha_p
               a1p[ip] = 0;
               a2p[ip] = 0;
               a3p[ip] = 0;

            }   //ip
         }   //iz
      }   //ir

      //Special case for ir == 0
      // use: Dxx/(R[i-1/2] dR)*{Ri(U[i+1]-Ui)-R[i-1](Ui-U[i-1])}; R[i-1]=0== 0
      for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
         T* a1r = alpha1_r+index(0,iz,0);
         T* a2r = alpha2_r+index(0,iz,0);
         T* a3r = alpha3_r+index(0,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a2r, a1r, a3r:64)
#endif
         for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
            a1r[ip] = 0;
            a3r[ip] = a2r[ip];
         }
      }

      // modification in propagation equation if diffusion coefficient depends on z,r
#pragma omp parallel for default(shared) schedule(static)
      for (size_t ir = 0; ir < size_t(particle.n_rgrid)-1; ++ir) {
         //We add first element for NUMA awareness and memory access
         if (ir == 0)
            continue;

         for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
            T* a1z = alpha1_z+index(ir,iz,0);
            T* a3z = alpha3_z+index(ir,iz,0);
            T* a1r = alpha1_r+index(ir,iz,0);
            T* a3r = alpha3_r+index(ir,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a1z, a3z, a1r, a3r:64)
#endif
            for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
               const T dDxxdr=(particle.Dxx.d2[ir+1][iz].s[ip]-particle.Dxx.d2[ir-1][iz].s[ip]) *pow(2.*particle.dr,-2.)*factor;
               const T dDxxdz=(particle.Dxx.d2[ir][iz+1].s[ip]-particle.Dxx.d2[ir][iz-1].s[ip]) *pow(2.*particle.dz,-2.)*factor;

               a1r[ip] -= dDxxdr*(1.-T(particle.dr/particle.r[ir]));
               a3r[ip] += dDxxdr*(1.+T(particle.dr/particle.r[ir]));
               a1z[ip] -= dDxxdz;
               a3z[ip] += dDxxdz;
            }//ip
         }//iz
      }//ir  


      // CONVECTION term: method according to ApJ 509, 212 with corrections  IMOS20010307
      // convection introduced in v38 using imos code with modifications AWS20010325
      if (galdef.convection) {
         //Reorder the loop to be more memory access efficient
#pragma omp parallel for default(shared) schedule(static)
         for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
            for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
               const T az1 = (particle.z[iz] > 0. && fabs(particle.z[iz]) >= 0.5*particle.dz ) ? 
                  T((galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz-1]))*1.e5/particle.dz/kpc2cm) : 0.;
               const T az2 = (fabs(particle.z[iz]) >= 0.5*particle.dz) ?
                  T((galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz]  ))*1.e5/particle.dz/kpc2cm) : 0;
               const T az3 = (particle.z[iz] < 0. && fabs(particle.z[iz]) >= 0.5*particle.dz ) ? 
                  T((galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz+1]))*1.e5/particle.dz/kpc2cm) : 0.;

               T* a1z = alpha1_z+index(ir,iz,0);
               T* a2z = alpha2_z+index(ir,iz,0);
               T* a3z = alpha3_z+index(ir,iz,0);
               T* a2p = alpha2_p+index(ir,iz,0);
               T* a3p = alpha3_p+index(ir,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a1z, a2z, a3z, a2p, a3p:64)
#endif
               for (size_t ip = 0; ip < size_t(particle.n_pgrid)-1; ++ip) {
                  const T ap2 = T(particle.p[ip]  /3.*galdef.dvdz_conv/(particle.p[ip+1]-particle.p[ip])*1.e5/kpc2cm);
	  
                  const T ap3 = T(particle.p[ip+1]/3.*galdef.dvdz_conv/(particle.p[ip+1]-particle.p[ip])*1.e5/kpc2cm);

                  a1z[ip] += az1;
                  a2z[ip] += az2;
                  a3z[ip] += az3;

                  a2p[ip] += ap2;
                  a3p[ip] += ap3;
               }
               //Add the last momentum element
               const size_t ip = particle.n_pgrid-1;
               const T ap2 = T(particle.p[ip]  /3.*galdef.dvdz_conv/(particle.p[ip]-particle.p[ip-1])*1.e5/kpc2cm);

               alpha1_z[index(ir,iz,ip)] += az1;
               alpha2_z[index(ir,iz,ip)] += az2;
               alpha3_z[index(ir,iz,ip)] += az3;

               alpha2_p[index(ir,iz,ip)] += ap2;
            }
         }
      }

      // DIFFUSIVE REACCELERATION
      if (galdef.diff_reacc > 0) { 
#pragma omp parallel for default(shared) schedule(static)
         for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
            for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
               T* a1p = alpha1_p+index(ir,iz,0);
               T* a2p = alpha2_p+index(ir,iz,0);
               T* a3p = alpha3_p+index(ir,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a1p, a2p, a3p:64)
#endif
               for (size_t ip = 1; ip < size_t(particle.n_pgrid)-1; ++ip) {

                  // alternative scheme #2, most detailed
                  a1p[ip] += 
                     T((-(particle.Dpp.d2[ir][iz].s[ip] - particle.Dpp.d2[ir][iz].s[ip-1])
                      /(particle.p[ip] - particle.p[ip-1])
                      + 2.*particle.Dpp.d2[ir][iz].s[ip]
                      /(particle.p[ip+1] - particle.p[ip-1])
                      + 2.*particle.Dpp.d2[ir][iz].s[ip-1]/particle.p[ip-1]) 
                     /(particle.p[ip]-particle.p[ip-1]));

                  a2p[ip] +=
                     T(-(particle.Dpp.d2[ir][iz].s[ip]-particle.Dpp.d2[ir][iz].s[ip-1])
                     /pow(particle.p[ip] - particle.p[ip-1], 2.)
                     + 2.*particle.Dpp.d2[ir][iz].s[ip]
                     /(particle.p[ip+1] - particle.p[ip-1])
                     *(1./(particle.p[ip+1] - particle.p[ip])
                           + 1./(particle.p[ip] - particle.p[ip-1]))
                     + 2.*particle.Dpp.d2[ir][iz].s[ip]
                     /(particle.p[ip] - particle.p[ip-1])
                     /particle.p[ip]);

                  a3p[ip] +=
                     T(2.*particle.Dpp.d2[ir][iz].s[ip]
                     /(particle.p[ip+1] - particle.p[ip-1])
                     /(particle.p[ip+1] - particle.p[ip]));

               }   //ip

               //Set boundaries for p
               alpha1_p[index(ir,iz,0)] += alpha1_p[index(ir,iz,1)];
               alpha2_p[index(ir,iz,0)] += alpha2_p[index(ir,iz,1)];
               alpha3_p[index(ir,iz,0)] += alpha3_p[index(ir,iz,1)];
               alpha1_p[index(ir,iz,particle.n_pgrid-1)] += alpha1_p[index(ir,iz,particle.n_pgrid-2)];
               alpha2_p[index(ir,iz,particle.n_pgrid-1)] += alpha2_p[index(ir,iz,particle.n_pgrid-2)];
               alpha3_p[index(ir,iz,particle.n_pgrid-1)] += alpha3_p[index(ir,iz,particle.n_pgrid-2)];

            }   //iz
         }   //ir
      }   // diffusive reacceleration

      // MOMENTUM LOSSES  IMOS20010307 minor change to make as in ApJ 509, 212
      // AWS20010622 correction to alpha2, consistent with ApJ 509, 212 Table 3 
      // AWS20010622 but note signs of alpha2, alpha3 incorrect in ApJ 509, 212 Table 3
      // AWS20010622 code is correct since dpdt = momentum loss rate is code
      if (galdef.momentum_losses) {
#pragma omp parallel for default(shared) schedule(static)
         for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
            for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
               T* a2p = alpha2_p+index(ir,iz,0);
               T* a3p = alpha3_p+index(ir,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a2p, a3p:64)
#endif
               for (size_t ip = 0; ip < size_t(particle.n_pgrid)-1; ++ip) {
                  a2p[ip] += T(particle.dpdt.d2[ir][iz].s[ip]/(particle.p[ip+1] - particle.p[ip]));
                  a3p[ip] += T(particle.dpdt.d2[ir][iz].s[ip+1]/(particle.p[ip+1] - particle.p[ip]));

               }   //ip

               const size_t ip = particle.n_pgrid-1;
               alpha2_p[index(ir,iz,ip)] += T(particle.dpdt.d2[ir][iz].s[ip]/(particle.p[ip] - particle.p[ip-1]));
               //Alpha 3 does not need to be set, phi vanishes implicitly beyond the boundary.
            }   //iz
         }   //ir 
      }   //momentum losses

      const T f_use = 1./(galdef.prop_r + galdef.prop_z + galdef.prop_p); // used only by method = 1

      //Add fragmentation
      if (galdef.fragmentation) {
#pragma omp parallel for default(shared) schedule(static)
         for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
            for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
               T* a2z = alpha2_z+index(ir,iz,0);
               T* a2r = alpha2_r+index(ir,iz,0);
               T* a2p = alpha2_p+index(ir,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a2z, a2r, a2p:64)
#endif
               for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                  const T value = T(particle.fragment.d2[ir][iz].s[ip]*f_use);
                  a2r[ip] += value;
                  a2z[ip] += value;
                  a2p[ip] += value;
               }
            }
         }
      }

      //Add decay
      if (particle.t_half != 0 && galdef.radioactive_decay) {
         std::ostringstream rBuf;
         rBuf << "Radioactive decay of " << particle.name << " with half life " << particle.t_half;
         INFO(rBuf.str());
#pragma omp parallel for default(shared) schedule(static)
         for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
            for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
               T* a2z = alpha2_z+index(ir,iz,0);
               T* a2r = alpha2_r+index(ir,iz,0);
               T* a2p = alpha2_p+index(ir,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a2z, a2r, a2p:64)
#endif
               for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                  const T value = T(particle.decay.d2[ir][iz].s[ip]*f_use);
                  a2r[ip] += value;
                  a2z[ip] += value;
                  a2p[ip] += value;
               }
            }
         }
      }

      return true;

   } //End 2D

   else if (3 == particle.n_spatial_dimensions) {

      //Create an indexer
      const IndexConverter3D index(particle.n_xgrid, particle.n_ygrid, particle.n_zgrid, particle.n_pgrid);

      // DIFFUSION term: method according to ApJ 509, 212
#pragma omp parallel for default(shared) schedule(static)
      for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
         for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
            for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
               T* a1z = alpha1_z+index(ix,iy,iz,0);
               T* a2z = alpha2_z+index(ix,iy,iz,0);
               T* a3z = alpha3_z+index(ix,iy,iz,0);
               T* a1x = alpha1_x+index(ix,iy,iz,0);
               T* a2x = alpha2_x+index(ix,iy,iz,0);
               T* a3x = alpha3_x+index(ix,iy,iz,0);
               T* a1y = alpha1_y+index(ix,iy,iz,0);
               T* a2y = alpha2_y+index(ix,iy,iz,0);
               T* a3y = alpha3_y+index(ix,iy,iz,0);
               T* a1p = alpha1_p+index(ix,iy,iz,0);
               T* a2p = alpha2_p+index(ix,iy,iz,0);
               T* a3p = alpha3_p+index(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a1z, a2z, a3z, a1x, a2x, a3x, a1y, a2y, a3y, a1p, a2p, a3p:64)
#endif
               for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                  const T xvalue = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dx,-2.)*factor;
                  a1x[ip] = xvalue;
                  a2x[ip] = xvalue*2.;
                  a3x[ip] = xvalue;

                  const T yvalue = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dy,-2.)*factor;
                  a1y[ip] = yvalue;
                  a2y[ip] = yvalue*2.;
                  a3y[ip] = yvalue;

                  const T zvalue = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dz,-2.)*factor;
                  a1z[ip] = zvalue;
                  a2z[ip] = zvalue*2.;
                  a3z[ip] = zvalue;

                  //Initialize alpha_p
                  a1p[ip] = 0;
                  a2p[ip] = 0;
                  a3p[ip] = 0;

               }   //ip
            }   //iz
         }   //iy
      }   //ix

#pragma omp parallel for default(shared) schedule(static)
      for (size_t ix = 1; ix < size_t(particle.n_xgrid-1); ++ix) {
         for (size_t iy = 1; iy < size_t(particle.n_ygrid-1); ++iy) {
            for (size_t iz = 1; iz < size_t(particle.n_zgrid-1); ++iz) {
               T* a1z = alpha1_z+index(ix,iy,iz,0);
               T* a3z = alpha3_z+index(ix,iy,iz,0);
               T* a1x = alpha1_x+index(ix,iy,iz,0);
               T* a3x = alpha3_x+index(ix,iy,iz,0);
               T* a1y = alpha1_y+index(ix,iy,iz,0);
               T* a3y = alpha3_y+index(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a1z, a3z, a1x, a3x, a1y, a3y:64)
#endif
               for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                  const T dDxxdx=(particle.Dxx.d3[ix+1][iy][iz].s[ip]-particle.Dxx.d3[ix-1][iy][iz].s[ip]) *pow(2.*particle.dx,-2.)*factor;
                  const T dDxxdy=(particle.Dxx.d3[ix][iy+1][iz].s[ip]-particle.Dxx.d3[ix][iy-1][iz].s[ip]) *pow(2.*particle.dy,-2.)*factor;
                  const T dDxxdz=(particle.Dxx.d3[ix][iy][iz+1].s[ip]-particle.Dxx.d3[ix][iy][iz-1].s[ip]) *pow(2.*particle.dz,-2.)*factor;

                  a1x[ip] -= dDxxdx;
                  a3x[ip] += dDxxdx;

                  a1y[ip] -= dDxxdy;
                  a3y[ip] += dDxxdy;

                  a1z[ip] -= dDxxdz;
                  a3z[ip] += dDxxdz;

               }   //ip
            }   //iz
         }   //iy
      }   //ix

      // CONVECTION term: method according to ApJ 509, 212 with corrections  IMOS20010307
      // convection introduced in v38 using imos code with modifications AWS20010325
      if (galdef.convection) {
         //Reorder the loop to be more memory access efficient
#pragma omp parallel for default(shared) schedule(static)
         for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
            for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
               for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
                  const T az1 = (particle.z[iz] > 0. && fabs(particle.z[iz]) >= 0.5*particle.dz ) ? 
                     T((galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz-1]))*1.e5/particle.dz/kpc2cm) : 0.;
                  const T az2 = (fabs(particle.z[iz]) >= 0.5*particle.dz) ?
                     T((galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz]  ))*1.e5/particle.dz/kpc2cm) : 0;
                  const T az3 = (particle.z[iz] < 0. && fabs(particle.z[iz]) >= 0.5*particle.dz ) ? 
                     T((galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz+1]))*1.e5/particle.dz/kpc2cm) : 0.;

                  T* a1z = alpha1_z+index(ix,iy,iz,0);
                  T* a2z = alpha2_z+index(ix,iy,iz,0);
                  T* a3z = alpha3_z+index(ix,iy,iz,0);
                  T* a2p = alpha2_p+index(ix,iy,iz,0);
                  T* a3p = alpha3_p+index(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a1z, a2z, a3z, a2p, a3p:64)
#endif
                  for (size_t ip = 0; ip < size_t(particle.n_pgrid)-1; ++ip) {
                     const T ap2 = T(particle.p[ip]  /3.*galdef.dvdz_conv/(particle.p[ip+1]-particle.p[ip])*1.e5/kpc2cm);

                     const T ap3 = T(particle.p[ip+1]/3.*galdef.dvdz_conv/(particle.p[ip+1]-particle.p[ip])*1.e5/kpc2cm);

                     a1z[ip] += az1;
                     a2z[ip] += az2;
                     a3z[ip] += az3;

                     a2p[ip] += ap2;
                     a3p[ip] += ap3;
                  }
                  //Add the last momentum element
                  const size_t ip = particle.n_pgrid-1;
                  const T ap2 = T(particle.p[ip]  /3.*galdef.dvdz_conv/(particle.p[ip]-particle.p[ip-1])*1.e5/kpc2cm);

                  alpha1_z[index(ix,iy,iz,ip)] += az1;
                  alpha2_z[index(ix,iy,iz,ip)] += az2;
                  alpha3_z[index(ix,iy,iz,ip)] += az3;

                  alpha2_p[index(ix,iy,iz,ip)] += ap2;
               }
            }
         }
      }

      // DIFFUSIVE REACCELERATION
      if (galdef.diff_reacc > 0) { 
#pragma omp parallel for default(shared) schedule(static)
         for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
            for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
               for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
                  T* a1p = alpha1_p+index(ix,iy,iz,0);
                  T* a2p = alpha2_p+index(ix,iy,iz,0);
                  T* a3p = alpha3_p+index(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a1p, a2p, a3p:64)
#endif
                  for (size_t ip = 1; ip < size_t(particle.n_pgrid)-1; ++ip) {

                     // alternative scheme #2, most detailed
                     a1p[ip] += 
                        T((-(particle.Dpp.d3[ix][iy][iz].s[ip] - particle.Dpp.d3[ix][iy][iz].s[ip-1])
                                 /(particle.p[ip] - particle.p[ip-1])
                                 + 2.*particle.Dpp.d3[ix][iy][iz].s[ip]
                                 /(particle.p[ip+1] - particle.p[ip-1])
                                 + 2.*particle.Dpp.d3[ix][iy][iz].s[ip-1]/particle.p[ip-1]) 
                              /(particle.p[ip]-particle.p[ip-1]));

                     a2p[ip] +=
                        T(-(particle.Dpp.d3[ix][iy][iz].s[ip]-particle.Dpp.d3[ix][iy][iz].s[ip-1])
                              /pow(particle.p[ip] - particle.p[ip-1], 2.)
                              + 2.*particle.Dpp.d3[ix][iy][iz].s[ip]
                              /(particle.p[ip+1] - particle.p[ip-1])
                              *(1./(particle.p[ip+1] - particle.p[ip])
                                 + 1./(particle.p[ip] - particle.p[ip-1]))
                              + 2.*particle.Dpp.d3[ix][iy][iz].s[ip]
                              /(particle.p[ip] - particle.p[ip-1])
                              /particle.p[ip]);

                     a3p[ip] +=
                        T(2.*particle.Dpp.d3[ix][iy][iz].s[ip]
                              /(particle.p[ip+1] - particle.p[ip-1])
                              /(particle.p[ip+1] - particle.p[ip]));

                  }   //ip

                  //Set boundaries for p
                  alpha1_p[index(ix,iy,iz,0)] += alpha1_p[index(ix,iy,iz,1)];
                  alpha2_p[index(ix,iy,iz,0)] += alpha2_p[index(ix,iy,iz,1)];
                  alpha3_p[index(ix,iy,iz,0)] += alpha3_p[index(ix,iy,iz,1)];
                  alpha1_p[index(ix,iy,iz,particle.n_pgrid-1)] += alpha1_p[index(ix,iy,iz,particle.n_pgrid-2)];
                  alpha2_p[index(ix,iy,iz,particle.n_pgrid-1)] += alpha2_p[index(ix,iy,iz,particle.n_pgrid-2)];
                  alpha3_p[index(ix,iy,iz,particle.n_pgrid-1)] += alpha3_p[index(ix,iy,iz,particle.n_pgrid-2)];

               }   //iz
            }   //iy
         }// ix
      }   // diffusive reacceleration

      // MOMENTUM LOSSES  IMOS20010307 minor change to make as in ApJ 509, 212
      // AWS20010622 correction to alpha2, consistent with ApJ 509, 212 Table 3 
      // AWS20010622 but note signs of alpha2, alpha3 incorrect in ApJ 509, 212 Table 3
      // AWS20010622 code is correct since dpdt = momentum loss rate is code
      if (galdef.momentum_losses) {
#pragma omp parallel for default(shared) schedule(static)
         for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
            for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
               for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
                  T* a2p = alpha2_p+index(ix,iy,iz,0);
                  T* a3p = alpha3_p+index(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a2p, a3p:64)
#endif
                  for (size_t ip = 0; ip < size_t(particle.n_pgrid)-1; ++ip) {
                     a2p[ip] += T(particle.dpdt.d3[ix][iy][iz].s[ip]/(particle.p[ip+1] - particle.p[ip]));
                     a3p[ip] += T(particle.dpdt.d3[ix][iy][iz].s[ip+1]/(particle.p[ip+1] - particle.p[ip]));

                  }   //ip

                  const size_t ip = particle.n_pgrid-1;
                  alpha2_p[index(ix,iy,iz,ip)] += T(particle.dpdt.d3[ix][iy][iz].s[ip]/(particle.p[ip] - particle.p[ip-1]));
                  //Alpha 3 does not need to be set, phi vanishes implicitly beyond the boundary.
               }   //iz
            }   //iy 
         }   //ix 
      }   //momentum losses

      const T f_use = 1./(galdef.prop_x + galdef.prop_y + galdef.prop_z + galdef.prop_p); // used only by method = 1

      //Add fragmentation
      if (galdef.fragmentation) {
#pragma omp parallel for default(shared) schedule(static)
         for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
            for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
               for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
                  T* a2z = alpha2_z+index(ix,iy,iz,0);
                  T* a2x = alpha2_x+index(ix,iy,iz,0);
                  T* a2y = alpha2_y+index(ix,iy,iz,0);
                  T* a2p = alpha2_p+index(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a2z, a2x, a2y, a2p:64)
#endif
                  for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                     const T value = T(particle.fragment.d3[ix][iy][iz].s[ip]*f_use);
                     a2x[ip] += value;
                     a2y[ip] += value;
                     a2z[ip] += value;
                     a2p[ip] += value;
                  }
               }
            }
         }
      }

      //Add decay
      if (particle.t_half != 0 && galdef.radioactive_decay) {
         std::ostringstream rBuf;
         rBuf << "Radioactive decay of " << particle.name << " with half life " << particle.t_half;
         INFO(rBuf.str());
#pragma omp parallel for default(shared) schedule(static)
         for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
            for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
               for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
                  T* a2z = alpha2_z+index(ix,iy,iz,0);
                  T* a2x = alpha2_x+index(ix,iy,iz,0);
                  T* a2y = alpha2_y+index(ix,iy,iz,0);
                  T* a2p = alpha2_p+index(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a2z, a2x, a2y, a2p:64)
#endif
                  for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                     const T value = T(particle.decay.d3[ix][iy][iz].s[ip]*f_use);
                     a2x[ip] += value;
                     a2y[ip] += value;
                     a2z[ip] += value;
                     a2p[ip] += value;
                  }
               }
            }
         }
      }

      return true;
   } //End 3D

   return false;

}

template <typename T>
PropelCrankNicolson<T>::PropelCrankNicolson(const Particle &particle, const Galdef &_galdef) :
   PropelOperatorSplitting<T>(particle, _galdef),
   total_source_function(nullptr),
   cr_density(nullptr),
   malpha1_z(nullptr), malpha1_r(nullptr), malpha1_x(nullptr), malpha1_y(nullptr),
   malpha2_z(nullptr), malpha2_r(nullptr), malpha2_x(nullptr), malpha2_y(nullptr),
   malpha3_z(nullptr), malpha3_r(nullptr), malpha3_x(nullptr), malpha3_y(nullptr),
   mz_tsf(nullptr), mr_tsf(nullptr), mx_tsf(nullptr), my_tsf(nullptr)
{
   //Allocate the memory
   if (2 == particle.n_spatial_dimensions) { // ==== 2D ====

      total_source_function = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      cr_density = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha1_z = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha2_z = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha3_z = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha1_r = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha2_r = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha3_r = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      mz_tsf = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      mr_tsf = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));

   }
   else if (3 == particle.n_spatial_dimensions) { // ==== 3D ====

      total_source_function = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      cr_density = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha1_z = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha2_z = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha3_z = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha1_x = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha2_x = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha3_x = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha1_y = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha2_y = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha3_y = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      mz_tsf = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      mx_tsf = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      my_tsf = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));

   }
}

template <typename T>
PropelCrankNicolson<T>::~PropelCrankNicolson()
{
   mlc::FreeAligned(total_source_function);
   mlc::FreeAligned(cr_density);
   mlc::FreeAligned(malpha1_z);
   mlc::FreeAligned(malpha2_z);
   mlc::FreeAligned(malpha3_z);
   mlc::FreeAligned(malpha1_r);
   mlc::FreeAligned(malpha2_r);
   mlc::FreeAligned(malpha3_r);
   mlc::FreeAligned(malpha1_x);
   mlc::FreeAligned(malpha2_x);
   mlc::FreeAligned(malpha3_x);
   mlc::FreeAligned(malpha1_y);
   mlc::FreeAligned(malpha2_y);
   mlc::FreeAligned(malpha3_y);
   mlc::FreeAligned(mz_tsf);
   mlc::FreeAligned(mr_tsf);
   mlc::FreeAligned(mx_tsf);
   mlc::FreeAligned(my_tsf);
}

template <typename T>
void PropelCrankNicolson<T>::operator() (Particle &particle) 
{
   //Initialize the alpha and abort if necessary
   if ( ! PropelOperatorSplitting<T>::initializeAlpha(particle) )
      return;

   //Initialize the timestep
   const T timestepFactor = T(PropelBase::galdef.timestep_factor);
   T dt = PropelBase::galdef.start_timestep*year2sec;
   const T end_timestep_sec = PropelBase::galdef.end_timestep*year2sec;
   
   if (2 == particle.n_spatial_dimensions) { // ==== 2D ====

      //Needed by the method
      const T f_use = 1./(PropelBase::galdef.prop_r + PropelBase::galdef.prop_z + PropelBase::galdef.prop_p);

      //Create the index functions
      const typename PropelOperatorSplitting<T>::IndexConverter2D indp(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
      const IndexConverter2Dz indz(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
      const IndexConverter2Dr indr(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);

      //populate it
#pragma omp parallel for default(shared) schedule(static)
      for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
         for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
            for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
               total_source_function[indp(ir,iz,ip)] = T(particle.primary_source_function.d2[ir][iz].s[ip]+particle.secondary_source_function.d2[ir][iz].s[ip]);
               cr_density[indp(ir,iz,ip)] = T(particle.cr_density.d2[ir][iz].s[ip]);
            }
         }
      }

      //Need to look into these timestep modes Andy introduced
      //Ignore them for now
      
      //Storage and temporary reordering of alpha arrays and source function
      //Only need to store for z and r propagation 
      //This could be avoided by overloading initilizeAlpha, but then we'd
      //need some mechanism for not copying the code because that is difficult to maintain properly
      //Do the reordering using NUMA aware looping
#pragma omp parallel for default(shared) schedule(static) 
      for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
         for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
            for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
               malpha1_z[indz(ir,iz,ip)] = alpha1_z[indp(ir,iz,ip)];
               malpha2_z[indz(ir,iz,ip)] = alpha2_z[indp(ir,iz,ip)];
               malpha3_z[indz(ir,iz,ip)] = alpha3_z[indp(ir,iz,ip)];
               mz_tsf[indz(ir,iz,ip)] = total_source_function[indp(ir,iz,ip)];
            }
         }
      }
      //Do the reordering using NUMA aware looping
#pragma omp parallel for default(shared) schedule(static)
      for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
         for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
            for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
               malpha1_r[indr(ir,iz,ip)] = alpha1_r[indp(ir,iz,ip)];
               malpha2_r[indr(ir,iz,ip)] = alpha2_r[indp(ir,iz,ip)];
               malpha3_r[indr(ir,iz,ip)] = alpha3_r[indp(ir,iz,ip)];
               mr_tsf[indr(ir,iz,ip)] = total_source_function[indp(ir,iz,ip)];
            }
         }
      }

      std::ostringstream tBuf;
      tBuf << "galdef_ID " << PropelBase::galdef.galdef_ID << " start timestep dt = " << dt/year2sec << " yr " << particle.name << endl;
      INFO(tBuf.str());

      while ( dt > end_timestep_sec*0.9999 ) {

         for (size_t irept = 1; irept <= size_t(PropelBase::galdef.timestep_repeat); ++irept) {
            // z propagation
            if (1 == PropelBase::galdef.prop_z) {
#pragma omp parallel default(shared) 
               {
                  T* Nz0 = static_cast<T*>(mlc::MallocAligned(particle.n_zgrid*sizeof(T)));
                  T* Nz1 = static_cast<T*>(mlc::MallocAligned(particle.n_zgrid*sizeof(T)));
                  T* Nz2 = static_cast<T*>(mlc::MallocAligned(particle.n_zgrid*sizeof(T)));
                  T* Nz3 = static_cast<T*>(mlc::MallocAligned(particle.n_zgrid*sizeof(T)));
                  T* Rz = static_cast<T*>(mlc::MallocAligned(particle.n_zgrid*sizeof(T)));
                  T* gamz = static_cast<T*>(mlc::MallocAligned(particle.n_zgrid*sizeof(T)));
#pragma omp for schedule(static)
                  for (size_t ir = 0; ir < size_t(particle.n_rgrid)-1; ++ir) {
                     //There is no need to do the last slice, because the boundary condition forces it to 0
                     for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                        T* a1 = malpha1_z+indz(ir,0,ip);
                        T* a2 = malpha2_z+indz(ir,0,ip);
                        T* a3 = malpha3_z+indz(ir,0,ip);
                        T* tsf = mz_tsf+indz(ir,0,ip);
                        //Cannot guarantee alignment of the a's and tsf because it may be a whole multiple of n_zgrid which is not always divisible with 8
#if _OPENMP >= 201307
#pragma omp simd aligned(Nz1, Nz2, Nz3, Nz0:64)
#endif
                        for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
                           Nz1[iz] =    - 0.5*a1[iz]*dt;
                           Nz2[iz] = 1. + 0.5*a2[iz]*dt;
                           Nz3[iz] =    - 0.5*a3[iz]*dt;

                           Nz0[iz] = tsf[iz]*f_use*dt;
                        }
                        for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
                           Nz0[iz] += (1. - 0.5*malpha2_z[indz(ir,iz,ip)]*dt)*cr_density[indp(ir,iz,ip)];
                        }
                        for (size_t iz = 1; iz < size_t(particle.n_zgrid); ++iz) {
                           Nz0[iz] += 0.5*malpha1_z[indz(ir,iz,ip)]*dt*cr_density[indp(ir,iz-1,ip)];
                        }
                        for (size_t iz = 0; iz < size_t(particle.n_zgrid)-1; ++iz) {
                           Nz0[iz] += 0.5*malpha3_z[indz(ir,iz,ip)]*dt*cr_density[indp(ir,iz+1,ip)];
                        }

                        PropelOperatorSplitting<T>::tridiag(Nz1, Nz2, Nz3, Nz0, Rz, gamz, particle.n_zgrid);

                        for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
                           cr_density[indp(ir,iz,ip)] = (Rz[iz] < 0.) ? 0. : Rz[iz];
                        }

                        //Boundary condition
                        cr_density[indp(ir,0,ip)] = 0.;
                        cr_density[indp(ir,particle.n_zgrid-1,ip)] = 0.;

                     }
                  }
                  mlc::FreeAligned(Nz0);
                  mlc::FreeAligned(Nz1);
                  mlc::FreeAligned(Nz2);
                  mlc::FreeAligned(Nz3);
                  mlc::FreeAligned(Rz);
                  mlc::FreeAligned(gamz);
               }
            } // z propagation

            // momentum propagation
            if (1 == PropelBase::galdef.prop_p) {
#pragma omp parallel default(shared) 
               {
                  T* Np0 = static_cast<T*>(mlc::MallocAligned(particle.n_pgrid*sizeof(T)));
                  T* Np1 = static_cast<T*>(mlc::MallocAligned(particle.n_pgrid*sizeof(T)));
                  T* Np2 = static_cast<T*>(mlc::MallocAligned(particle.n_pgrid*sizeof(T)));
                  T* Np3 = static_cast<T*>(mlc::MallocAligned(particle.n_pgrid*sizeof(T)));
                  T* Rp = static_cast<T*>(mlc::MallocAligned(particle.n_pgrid*sizeof(T)));
                  T* gamp = static_cast<T*>(mlc::MallocAligned(particle.n_pgrid*sizeof(T)));
#pragma omp for schedule(static)
                  for (size_t ir = 0; ir < size_t(particle.n_rgrid)-1; ++ir) {
                     //There is no need to do the last slice, because the boundary condition forces it to 0
                     for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
                        //There is no need to do the first and last slices, because the boundary condition forces them to 0
                        T* a1 = alpha1_p+indp(ir,iz,0);
                        T* a2 = alpha2_p+indp(ir,iz,0);
                        T* a3 = alpha3_p+indp(ir,iz,0);
                        T* tsf = total_source_function+indp(ir,iz,0);
                        T* crd = cr_density+indp(ir,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a1, a2, a3, Np1, Np2, Np3, Np0, crd, tsf:64)
#endif
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           Np1[ip] =    - 0.5*a1[ip]*dt;
                           Np2[ip] = 1. + 0.5*a2[ip]*dt;
                           Np3[ip] =    - 0.5*a3[ip]*dt;

                           Np0[ip] = tsf[ip]*f_use*dt + (1. - 0.5*a2[ip]*dt)*crd[ip];
                        }
#if _OPENMP >= 201307
#pragma omp simd aligned(a1, Np0, crd:64)
#endif
                        for (size_t ip = 1; ip < size_t(particle.n_pgrid); ++ip) {
                           Np0[ip] += 0.5*a1[ip]*dt*crd[ip-1];
                        }
#if _OPENMP >= 201307
#pragma omp simd aligned(a3, Np0, crd:64)
#endif
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid)-1; ++ip) {
                           Np0[ip] += 0.5*a3[ip]*dt*crd[ip+1];
                        }

                        PropelOperatorSplitting<T>::tridiag(Np1, Np2, Np3, Np0, Rp, gamp, particle.n_pgrid);

#if _OPENMP >= 201307
#pragma omp simd aligned(Rp, crd:64)
#endif
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           crd[ip] = (Rp[ip] < 0.) ? 0. : Rp[ip];
                        }

                        //No boundary conditions here
                     }
                  }
                  mlc::FreeAligned(Np0);
                  mlc::FreeAligned(Np1);
                  mlc::FreeAligned(Np2);
                  mlc::FreeAligned(Np3);
                  mlc::FreeAligned(Rp);
                  mlc::FreeAligned(gamp);
               }
            } // momentum propagation

            // r propagation
            if (1 == PropelBase::galdef.prop_r) {
#pragma omp parallel default(shared) 
               {
                  T* Nr0 = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*sizeof(T)));
                  T* Nr1 = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*sizeof(T)));
                  T* Nr2 = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*sizeof(T)));
                  T* Nr3 = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*sizeof(T)));
                  T* Rr = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*sizeof(T)));
                  T* gamr = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*sizeof(T)));
#pragma omp for schedule(static)
                  for (size_t iz = 0; iz < size_t(particle.n_zgrid)-1; ++iz) {
                     //There is no need to do the last slice, because the boundary condition forces it to 0
                     if (iz == 0)
                        continue; //We do this with an if statment for NUMA purposes
                     for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                        T* a1 = malpha1_r+indr(0,iz,ip);
                        T* a2 = malpha2_r+indr(0,iz,ip);
                        T* a3 = malpha3_r+indr(0,iz,ip);
                        T* tsf = mr_tsf+indr(0,iz,ip);
#if _OPENMP >= 201307
#pragma omp simd aligned(a1, a2, a3, Nr1, Nr2, Nr3, Nr0, tsf:64)
#endif
                        for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
                           Nr1[ir] =    - 0.5*a1[ir]*dt;
                           Nr2[ir] = 1. + 0.5*a2[ir]*dt;
                           Nr3[ir] =    - 0.5*a3[ir]*dt;

                           Nr0[ir] = tsf[ir]*f_use*dt;
                        }
                        for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
                           Nr0[ir] += (1. - 0.5*malpha2_r[indr(ir,iz,ip)]*dt)*cr_density[indp(ir,iz,ip)];
                        }
                        for (size_t ir = 1; ir < size_t(particle.n_rgrid); ++ir) {
                           Nr0[ir] += 0.5*malpha1_r[indr(ir,iz,ip)]*dt*cr_density[indp(ir-1,iz,ip)];
                        }
                        for (size_t ir = 0; ir < size_t(particle.n_rgrid)-1; ++ir) {
                           Nr0[ir] += 0.5*malpha3_r[indr(ir,iz,ip)]*dt*cr_density[indp(ir+1,iz,ip)];
                        }

                        PropelOperatorSplitting<T>::tridiag(Nr1, Nr2, Nr3, Nr0, Rr, gamr, particle.n_rgrid);

                        for (size_t ir = 0; ir < size_t(particle.n_rgrid)-1; ++ir) {
                           cr_density[indp(ir,iz,ip)] = (Rr[ir] < 0.) ? 0. : Rr[ir];
                        }

                        //Boundary condition
                        cr_density[indp(particle.n_rgrid-1,iz,ip)] = 0.;

                     }
                  }
                  mlc::FreeAligned(Nr0);
                  mlc::FreeAligned(Nr1);
                  mlc::FreeAligned(Nr2);
                  mlc::FreeAligned(Nr3);
                  mlc::FreeAligned(Rr);
                  mlc::FreeAligned(gamr);
               }
            } // r propagation

            //The boundary conditions are properly set above, no need to re-apply them

         } //repeat

         //Fix dt for the next timestep 

         dt *= timestepFactor;
         tBuf.str("");
         tBuf<<" galdef_ID "<<PropelBase::galdef.galdef_ID<<" new timestep dt="<<dt/year2sec<<" yr  "
            <<particle.name;
         INFO(tBuf.str());

      }//Propagation

      //Assign the cr density
#pragma omp parallel for default(shared) schedule(static)
      for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
         for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
#if _OPENMP >= 201307
#pragma omp simd 
#endif
            for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
               particle.cr_density.d2[ir][iz].s[ip] = cr_density[indp(ir,iz,ip)];
            }
         }
      }
   } //2D

   else if (3 == particle.n_spatial_dimensions) { // ==== 3D ====

      //Needed by the method
      const T f_use = 1./(PropelBase::galdef.prop_x + PropelBase::galdef.prop_y + PropelBase::galdef.prop_z + PropelBase::galdef.prop_p);

      //Create the index functions
      const typename PropelOperatorSplitting<T>::IndexConverter3D indp(particle.n_xgrid, particle.n_ygrid, particle.n_zgrid, particle.n_pgrid);
      const IndexConverter3Dz indz(particle.n_xgrid, particle.n_ygrid, particle.n_zgrid, particle.n_pgrid);
      const IndexConverter3Dy indy(particle.n_xgrid, particle.n_ygrid, particle.n_zgrid, particle.n_pgrid);
      const IndexConverter3Dx indx(particle.n_xgrid, particle.n_ygrid, particle.n_zgrid, particle.n_pgrid);

      //populate it
#pragma omp parallel for default(shared) schedule(static)
      for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
         for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
            for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
#if _OPENMP >= 201307
#pragma omp simd 
#endif
               for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                  total_source_function[indp(ix,iy,iz,ip)] = T(particle.primary_source_function.d3[ix][iy][iz].s[ip]+particle.secondary_source_function.d3[ix][iy][iz].s[ip]);
                  cr_density[indp(ix,iy,iz,ip)] = T(particle.cr_density.d3[ix][iy][iz].s[ip]);
               }
            }
         }
      }

      //Need to look into these timestep modes Andy introduced
      //Ignore them for now
      
      //Storage and temporary reordering of alpha arrays and source function
      //Only need to store for x,y, and z propagation 
      //This could be avoided by overloading initilizeAlpha, but then we'd
      //need some mechanism for not copying the code because that is difficult to maintain properly

      //Do the reordering using NUMA aware looping
#pragma omp parallel for default(shared) schedule(static) 
      for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
         for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
            for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
               for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
                  malpha1_z[indz(ix,iy,iz,ip)] = alpha1_z[indp(ix,iy,iz,ip)];
                  malpha2_z[indz(ix,iy,iz,ip)] = alpha2_z[indp(ix,iy,iz,ip)];
                  malpha3_z[indz(ix,iy,iz,ip)] = alpha3_z[indp(ix,iy,iz,ip)];
                  mz_tsf[indz(ix,iy,iz,ip)] = total_source_function[indp(ix,iy,iz,ip)];
               }
            }
         }
      }
#pragma omp parallel for default(shared) schedule(static)
      for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
         for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
            for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
               for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
                  malpha1_y[indy(ix,iy,iz,ip)] = alpha1_y[indp(ix,iy,iz,ip)];
                  malpha2_y[indy(ix,iy,iz,ip)] = alpha2_y[indp(ix,iy,iz,ip)];
                  malpha3_y[indy(ix,iy,iz,ip)] = alpha3_y[indp(ix,iy,iz,ip)];
                  my_tsf[indy(ix,iy,iz,ip)] = total_source_function[indp(ix,iy,iz,ip)];
               }
            }
         }
      }

#pragma omp parallel for default(shared) schedule(static)
      for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
         for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
            for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
               for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
                  malpha1_x[indx(ix,iy,iz,ip)] = alpha1_x[indp(ix,iy,iz,ip)];
                  malpha2_x[indx(ix,iy,iz,ip)] = alpha2_x[indp(ix,iy,iz,ip)];
                  malpha3_x[indx(ix,iy,iz,ip)] = alpha3_x[indp(ix,iy,iz,ip)];
                  mx_tsf[indx(ix,iy,iz,ip)] = total_source_function[indp(ix,iy,iz,ip)];
               }
            }
         }
      }

      std::ostringstream tBuf;
      tBuf << "galdef_ID " << PropelBase::galdef.galdef_ID << " start timestep dt = " << dt/year2sec << " yr " << particle.name << endl;
      INFO(tBuf.str());

      while ( dt > end_timestep_sec*0.9999 ) {

         for (size_t irept = 1; irept <= size_t(PropelBase::galdef.timestep_repeat); ++irept) {
            // x propagation
            if (1 == PropelBase::galdef.prop_x) {
#pragma omp parallel default(shared) 
               {
                  T* Nx0 = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*sizeof(T)));
                  T* Nx1 = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*sizeof(T)));
                  T* Nx2 = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*sizeof(T)));
                  T* Nx3 = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*sizeof(T)));
                  T* Rx = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*sizeof(T)));
                  T* gamx = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*sizeof(T)));
#pragma omp for schedule(static)
                  for (size_t iy = 0; iy < size_t(particle.n_ygrid)-1; ++iy) {
                     if (iy == 0)
                        continue; //We do this with an if statment for NUMA purposes
                     for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
                        //There is no need to do the last slice, because the boundary condition forces it to 0
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           T* a1 = malpha1_x+indx(0,iy,iz,ip);
                           T* a2 = malpha2_x+indx(0,iy,iz,ip);
                           T* a3 = malpha3_x+indx(0,iy,iz,ip);
                           T* tsf = mx_tsf+indx(0,iy,iz,ip);
#if _OPENMP >= 201307
#pragma omp simd aligned(a1, a2, a3, Nx1, Nx2, Nx3, Nx0, tsf:64)
#endif
                           for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
                              Nx1[ix] =    - 0.5*a1[ix]*dt;
                              Nx2[ix] = 1. + 0.5*a2[ix]*dt;
                              Nx3[ix] =    - 0.5*a3[ix]*dt;

                              Nx0[ix] = tsf[ix]*f_use*dt;
                           }
                           for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
                              Nx0[ix] += (1. - 0.5*malpha2_x[indx(ix,iy,iz,ip)]*dt)*cr_density[indp(ix,iy,iz,ip)];
                           }
                           for (size_t ix = 1; ix < size_t(particle.n_xgrid); ++ix) {
                              Nx0[ix] += 0.5*malpha1_x[indx(ix,iy,iz,ip)]*dt*cr_density[indp(ix-1,iy,iz,ip)];
                           }
                           for (size_t ix = 0; ix < size_t(particle.n_xgrid)-1; ++ix) {
                              Nx0[ix] += 0.5*malpha3_x[indx(ix,iy,iz,ip)]*dt*cr_density[indp(ix+1,iy,iz,ip)];
                           }

                           PropelOperatorSplitting<T>::tridiag(Nx1, Nx2, Nx3, Nx0, Rx, gamx, particle.n_xgrid);

                           for (size_t ix = 1; ix < size_t(particle.n_xgrid)-1; ++ix) {
                              cr_density[indp(ix,iy,iz,ip)] = (Rx[ix] < 0.) ? 0. : Rx[ix];
                           }

                           //Boundary condition
                           cr_density[indp(0,iy,iz,ip)] = 0.;
                           cr_density[indp(particle.n_xgrid-1,iy,iz,ip)] = 0.;

                        }
                     }
                  }
                  mlc::FreeAligned(Nx0);
                  mlc::FreeAligned(Nx1);
                  mlc::FreeAligned(Nx2);
                  mlc::FreeAligned(Nx3);
                  mlc::FreeAligned(Rx);
                  mlc::FreeAligned(gamx);
               }
            } // x propagation

            // y propagation
            if (1 == PropelBase::galdef.prop_y) {
#pragma omp parallel default(shared) 
               {
                  T* Ny0 = static_cast<T*>(mlc::MallocAligned(particle.n_ygrid*sizeof(T)));
                  T* Ny1 = static_cast<T*>(mlc::MallocAligned(particle.n_ygrid*sizeof(T)));
                  T* Ny2 = static_cast<T*>(mlc::MallocAligned(particle.n_ygrid*sizeof(T)));
                  T* Ny3 = static_cast<T*>(mlc::MallocAligned(particle.n_ygrid*sizeof(T)));
                  T* Ry = static_cast<T*>(mlc::MallocAligned(particle.n_ygrid*sizeof(T)));
                  T* gamy = static_cast<T*>(mlc::MallocAligned(particle.n_ygrid*sizeof(T)));
#pragma omp for schedule(static)
                  for (size_t ix = 0; ix < size_t(particle.n_xgrid)-1; ++ix) {
                     if (ix == 0)
                        continue; //We do this with an if statment for NUMA purposes
                     for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
                        //There is no need to do the last slice, because the boundary condition forces it to 0
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           T* a1 = malpha1_y+indy(ix,0,iz,ip);
                           T* a2 = malpha2_y+indy(ix,0,iz,ip);
                           T* a3 = malpha3_y+indy(ix,0,iz,ip);
                           T* tsf = my_tsf+indy(ix,0,iz,ip);
                           //Cannot guarantee alignment of the ais at this point
#if _OPENMP >= 201307
#pragma omp simd aligned(Ny1, Ny2, Ny3, Ny0:64)
#endif
                           for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
                              Ny1[iy] =    - 0.5*a1[iy]*dt;
                              Ny2[iy] = 1. + 0.5*a2[iy]*dt;
                              Ny3[iy] =    - 0.5*a3[iy]*dt;

                              Ny0[iy] = tsf[iy]*f_use*dt;
                           }
                           for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
                              Ny0[iy] += (1. - 0.5*malpha2_y[indy(ix,iy,iz,ip)]*dt)*cr_density[indp(ix,iy,iz,ip)];
                           }
                           for (size_t iy = 1; iy < size_t(particle.n_ygrid); ++iy) {
                              Ny0[iy] += 0.5*malpha1_y[indy(ix,iy,iz,ip)]*dt*cr_density[indp(ix,iy-1,iz,ip)];
                           }
                           for (size_t iy = 0; iy < size_t(particle.n_ygrid)-1; ++iy) {
                              Ny0[iy] += 0.5*malpha3_y[indy(ix,iy,iz,ip)]*dt*cr_density[indp(ix,iy+1,iz,ip)];
                           }

                           PropelOperatorSplitting<T>::tridiag(Ny1, Ny2, Ny3, Ny0, Ry, gamy, particle.n_ygrid);

                           for (size_t iy = 1; iy < size_t(particle.n_ygrid)-1; ++iy) {
                              cr_density[indp(ix,iy,iz,ip)] = (Ry[iy] < 0.) ? 0. : Ry[iy];
                           }

                           //Boundary condition
                           cr_density[indp(ix,0,iz,ip)] = 0.;
                           cr_density[indp(ix,particle.n_ygrid-1,iz,ip)] = 0.;

                        }
                     }
                  }
                  mlc::FreeAligned(Ny0);
                  mlc::FreeAligned(Ny1);
                  mlc::FreeAligned(Ny2);
                  mlc::FreeAligned(Ny3);
                  mlc::FreeAligned(Ry);
                  mlc::FreeAligned(gamy);
               }
            } // y propagation

            // z propagation
            if (1 == PropelBase::galdef.prop_z) {
#pragma omp parallel default(shared) 
               {
                  T* Nz0 = static_cast<T*>(mlc::MallocAligned(particle.n_zgrid*sizeof(T)));
                  T* Nz1 = static_cast<T*>(mlc::MallocAligned(particle.n_zgrid*sizeof(T)));
                  T* Nz2 = static_cast<T*>(mlc::MallocAligned(particle.n_zgrid*sizeof(T)));
                  T* Nz3 = static_cast<T*>(mlc::MallocAligned(particle.n_zgrid*sizeof(T)));
                  T* Rz = static_cast<T*>(mlc::MallocAligned(particle.n_zgrid*sizeof(T)));
                  T* gamz = static_cast<T*>(mlc::MallocAligned(particle.n_zgrid*sizeof(T)));
#pragma omp for schedule(static)
                  for (size_t ix = 0; ix < size_t(particle.n_xgrid)-1; ++ix) {
                     if (ix == 0)
                        continue; //For NUMA, no need to do first slice
                     for (size_t iy = 1; iy < size_t(particle.n_ygrid)-1; ++iy) {
                        //There is no need to do the first and last slice, because the boundary condition forces it to 0
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           T* a1 = malpha1_z+indz(ix,iy,0,ip);
                           T* a2 = malpha2_z+indz(ix,iy,0,ip);
                           T* a3 = malpha3_z+indz(ix,iy,0,ip);
                           T* tsf = mz_tsf+indz(ix,iy,0,ip);
#if _OPENMP >= 201307
#pragma omp simd aligned(Nz1, Nz2, Nz3, Nz0:64)
#endif
                           for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
                              Nz1[iz] =    - 0.5*a1[iz]*dt;
                              Nz2[iz] = 1. + 0.5*a2[iz]*dt;
                              Nz3[iz] =    - 0.5*a3[iz]*dt;

                              Nz0[iz] = tsf[iz]*f_use*dt;
                           }
                           for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
                              Nz0[iz] += (1. - 0.5*malpha2_z[indz(ix,iy,iz,ip)]*dt)*cr_density[indp(ix,iy,iz,ip)];
                           }
                           for (size_t iz = 1; iz < size_t(particle.n_zgrid); ++iz) {
                              Nz0[iz] += 0.5*malpha1_z[indz(ix,iy,iz,ip)]*dt*cr_density[indp(ix,iy,iz-1,ip)];
                           }
                           for (size_t iz = 0; iz < size_t(particle.n_zgrid)-1; ++iz) {
                              Nz0[iz] += 0.5*malpha3_z[indz(ix,iy,iz,ip)]*dt*cr_density[indp(ix,iy,iz+1,ip)];
                           }

                           PropelOperatorSplitting<T>::tridiag(Nz1, Nz2, Nz3, Nz0, Rz, gamz, particle.n_zgrid);

                           for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
                              cr_density[indp(ix,iy,iz,ip)] = (Rz[iz] < 0.) ? 0. : Rz[iz];
                           }

                           //Boundary condition
                           cr_density[indp(ix,iy,0,ip)] = 0.;
                           cr_density[indp(ix,iy,particle.n_zgrid-1,ip)] = 0.;

                        }
                     }
                  }
                  mlc::FreeAligned(Nz0);
                  mlc::FreeAligned(Nz1);
                  mlc::FreeAligned(Nz2);
                  mlc::FreeAligned(Nz3);
                  mlc::FreeAligned(Rz);
                  mlc::FreeAligned(gamz);
               }
            } // z propagation

            // momentum propagation
            if (1 == PropelBase::galdef.prop_p) {
#pragma omp parallel default(shared) 
               {
                  T* Np0 = static_cast<T*>(mlc::MallocAligned(particle.n_pgrid*sizeof(T)));
                  T* Np1 = static_cast<T*>(mlc::MallocAligned(particle.n_pgrid*sizeof(T)));
                  T* Np2 = static_cast<T*>(mlc::MallocAligned(particle.n_pgrid*sizeof(T)));
                  T* Np3 = static_cast<T*>(mlc::MallocAligned(particle.n_pgrid*sizeof(T)));
                  T* Rp = static_cast<T*>(mlc::MallocAligned(particle.n_pgrid*sizeof(T)));
                  T* gamp = static_cast<T*>(mlc::MallocAligned(particle.n_pgrid*sizeof(T)));
#pragma omp for schedule(static)
                  for (size_t ix = 0; ix < size_t(particle.n_xgrid)-1; ++ix) {
                     if (ix == 0)
                        continue; //For NUMA, no need to do first slice
                     for (size_t iy = 1; iy < size_t(particle.n_ygrid)-1; ++iy) {
                        for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
                           //There is no need to do the first and last slices, because the boundary condition forces them to 0
                           T* a1 = alpha1_p+indp(ix,iy,iz,0);
                           T* a2 = alpha2_p+indp(ix,iy,iz,0);
                           T* a3 = alpha3_p+indp(ix,iy,iz,0);
                           T* tsf = total_source_function+indp(ix,iy,iz,0);
                           T* crd = cr_density+indp(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(a1, a2, a3, Np1, Np2, Np3, Np0, tsf, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                              Np1[ip] =    - 0.5*a1[ip]*dt;
                              Np2[ip] = 1. + 0.5*a2[ip]*dt;
                              Np3[ip] =    - 0.5*a3[ip]*dt;

                              Np0[ip] = tsf[ip]*f_use*dt + (1. - 0.5*a2[ip]*dt)*crd[ip];
                           }
#if _OPENMP >= 201307
#pragma omp simd aligned(a1, Np0, crd:64)
#endif
                           for (size_t ip = 1; ip < size_t(particle.n_pgrid); ++ip) {
                              Np0[ip] += 0.5*a1[ip]*dt*crd[ip-1];
                           }
#if _OPENMP >= 201307
#pragma omp simd aligned(a3, Np0, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid)-1; ++ip) {
                              Np0[ip] += 0.5*a3[ip]*dt*crd[ip+1];
                           }

                           PropelOperatorSplitting<T>::tridiag(Np1, Np2, Np3, Np0, Rp, gamp, particle.n_pgrid);

#if _OPENMP >= 201307
#pragma omp simd aligned(Rp, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                              crd[ip] = (Rp[ip] < 0.) ? 0. : Rp[ip];
                           }

                           //No boundary conditions here
                        }
                     }
                  }
                  mlc::FreeAligned(Np0);
                  mlc::FreeAligned(Np1);
                  mlc::FreeAligned(Np2);
                  mlc::FreeAligned(Np3);
                  mlc::FreeAligned(Rp);
                  mlc::FreeAligned(gamp);
               }
            } // momentum propagation


            //The boundary conditions are properly set above, no need to re-apply them

         } //repeat

         //Fix dt for the next timestep 

         dt *= timestepFactor;
         tBuf.str("");
         tBuf<<" galdef_ID "<<PropelBase::galdef.galdef_ID<<" new timestep dt="<<dt/year2sec<<" yr  "
            <<particle.name;
         INFO(tBuf.str());

      }//Propagation

      //Assign the cr density
#pragma omp parallel for default(shared) schedule(static)
      for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
         for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
            for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
#if _OPENMP >= 201307
#pragma omp simd 
#endif
               for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                  particle.cr_density.d3[ix][iy][iz].s[ip] = cr_density[indp(ix,iy,iz,ip)];
               }
            }
         }
      }
   } //3D
}

//Explicit instantiation
template class PropelCrankNicolson<float>;
template class PropelCrankNicolson<double>;


//We make the assumption here in the code that n_pgrid, n_rgrid, and n_xgrid are a whole multiple of 8
//so we are guaranteed to get the pointers in the innermost loop aligned on 64 bytes.

//For best performance the n_xgrid, n_ygrid, n_rgrid, and n_zgrid should be a whole multiple of number of threads as well

template <typename T>
PropelCrankNicolsonVector<T>::PropelCrankNicolsonVector(const Particle &particle, const Galdef &_galdef) :
   PropelOperatorSplitting<T>(particle, _galdef),
   total_source_function(nullptr),
   cr_density(nullptr),
   malpha1_p(nullptr), malpha1_r(nullptr), malpha1_x(nullptr), malpha1_y(nullptr),
   malpha2_p(nullptr), malpha2_r(nullptr), malpha2_x(nullptr), malpha2_y(nullptr),
   malpha3_p(nullptr), malpha3_r(nullptr), malpha3_x(nullptr), malpha3_y(nullptr),
   mp_tsf(nullptr), mr_tsf(nullptr), mx_tsf(nullptr), my_tsf(nullptr)
{
   //Allocate the memory
   if (2 == particle.n_spatial_dimensions) { // ==== 2D ====

      total_source_function = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      cr_density = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha1_p = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha2_p = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha3_p = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha1_r = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha2_r = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha3_r = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      mp_tsf = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      mr_tsf = static_cast<T*>(mlc::MallocAligned(particle.n_rgrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));

   }
   else if (3 == particle.n_spatial_dimensions) { // ==== 3D ====

      total_source_function = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      cr_density = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha1_p = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha2_p = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha3_p = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha1_x = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha2_x = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha3_x = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha1_y = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha2_y = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      malpha3_y = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      mp_tsf = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      mx_tsf = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));
      my_tsf = static_cast<T*>(mlc::MallocAligned(particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid*sizeof(T)));

   }
}

template <typename T>
PropelCrankNicolsonVector<T>::~PropelCrankNicolsonVector()
{
   mlc::FreeAligned(total_source_function);
   mlc::FreeAligned(cr_density);
   mlc::FreeAligned(malpha1_p);
   mlc::FreeAligned(malpha2_p);
   mlc::FreeAligned(malpha3_p);
   mlc::FreeAligned(malpha1_r);
   mlc::FreeAligned(malpha2_r);
   mlc::FreeAligned(malpha3_r);
   mlc::FreeAligned(malpha1_x);
   mlc::FreeAligned(malpha2_x);
   mlc::FreeAligned(malpha3_x);
   mlc::FreeAligned(malpha1_y);
   mlc::FreeAligned(malpha2_y);
   mlc::FreeAligned(malpha3_y);
   mlc::FreeAligned(mp_tsf);
   mlc::FreeAligned(mr_tsf);
   mlc::FreeAligned(mx_tsf);
   mlc::FreeAligned(my_tsf);
}

template <typename T>
void PropelCrankNicolsonVector<T>::operator() (Particle &particle) 
{
   //Initialize the alpha and abort if necessary
   if ( ! PropelOperatorSplitting<T>::initializeAlpha(particle) )
      return;

   //Initialize the timestep
   const T timestepFactor = T(PropelBase::galdef.timestep_factor);
   const T end_timestep_sec = PropelBase::galdef.end_timestep*year2sec;
   
   if (2 == particle.n_spatial_dimensions) { // ==== 2D ====

      //Needed by the method
      const T f_use = 1./(PropelBase::galdef.prop_r + PropelBase::galdef.prop_z + PropelBase::galdef.prop_p);

      //Create the index functions
      const typename PropelOperatorSplitting<T>::IndexConverter2D indz(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
      const IndexConverter2Dp indp(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
      const IndexConverter2Dr indr(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);

      //populate it
#pragma omp parallel for default(shared) schedule(static)
      for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
         for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
            T* tsf = total_source_function+indz(ir,iz,0);
            T* crd = cr_density+indz(ir,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(tsf, crd:64)
#endif
            for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
               tsf[ip] = T(particle.primary_source_function.d2[ir][iz].s[ip]+particle.secondary_source_function.d2[ir][iz].s[ip]);
               crd[ip] = T(particle.cr_density.d2[ir][iz].s[ip]);
            }
         }
      }

      //Need to look into these timestep modes Andy introduced
      //Ignore them for now
      
      //Storage and temporary reordering of alpha arrays and source function
      //Only need to store for z and r propagation 
      //This could be avoided by overloading initilizeAlpha, but then we'd
      //need some mechanism for not copying the code because that is difficult to maintain properly

      //Do the reordering using NUMA aware looping
#pragma omp parallel for default(shared) schedule(static) 
      for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
         for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
            //There is no benefits in simd here because the right side is not continuous
            for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
               malpha1_p[indp(ir,iz,ip)] = alpha1_p[indz(ir,iz,ip)];
               malpha2_p[indp(ir,iz,ip)] = alpha2_p[indz(ir,iz,ip)];
               malpha3_p[indp(ir,iz,ip)] = alpha3_p[indz(ir,iz,ip)];
               mp_tsf[indp(ir,iz,ip)] = total_source_function[indz(ir,iz,ip)];
            }
         }
      }
      //Do the reordering using NUMA aware looping
#pragma omp parallel for default(shared) schedule(static)
      for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
         for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
            T* ma1 = malpha1_r+indr(ir,iz,0);
            T* ma2 = malpha2_r+indr(ir,iz,0);
            T* ma3 = malpha3_r+indr(ir,iz,0);
            T* a1 = alpha1_r+indz(ir,iz,0);
            T* a2 = alpha2_r+indz(ir,iz,0);
            T* a3 = alpha3_r+indz(ir,iz,0);
            T* mtsf = mr_tsf+indr(ir,iz,0);
            T* tsf = total_source_function+indz(ir,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(ma1, ma2, ma3, mtsf, a1, a2, a3, tsf:64)
#endif
            for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
               ma1[ip] = a1[ip];
               ma2[ip] = a2[ip];
               ma3[ip] = a3[ip];
               mtsf[ip] = tsf[ip];
            }
         }
      }

      //Sort the grid sizes for easier allocations in the parallel region
      std::vector<size_t> sortedGridSizes({size_t(particle.n_rgrid), size_t(particle.n_zgrid), size_t(particle.n_pgrid)});
      std::sort(sortedGridSizes.begin(), sortedGridSizes.end());


      //Have the entire loop a big parallel region.  Almost anything is run in parallel within anyway, reduces the overhead
#pragma omp parallel default(shared) 
      {
         //Storage for tridag, use a single one just big enough for all variations
         T* N0_store = static_cast<T*>(mlc::MallocAligned(sortedGridSizes[2]*sortedGridSizes[1]*sizeof(T)));
         T* N1_store = static_cast<T*>(mlc::MallocAligned(sortedGridSizes[2]*sortedGridSizes[1]*sizeof(T)));
         T* N2_store = static_cast<T*>(mlc::MallocAligned(sortedGridSizes[2]*sortedGridSizes[1]*sizeof(T)));
         T* N3_store = static_cast<T*>(mlc::MallocAligned(sortedGridSizes[2]*sortedGridSizes[1]*sizeof(T)));
         T* R_store = static_cast<T*>(mlc::MallocAligned(sortedGridSizes[2]*sortedGridSizes[1]*sizeof(T)));
         T* gam_store = static_cast<T*>(mlc::MallocAligned(sortedGridSizes[2]*sortedGridSizes[1]*sizeof(T)));
         T* bet_store = static_cast<T*>(mlc::MallocAligned(sortedGridSizes[2]*sizeof(T)));

         for (T dt = PropelBase::galdef.start_timestep*year2sec; dt > end_timestep_sec*0.9999; dt*= timestepFactor) {


#pragma omp single
            {
               std::ostringstream tBuf;
               tBuf<<" galdef_ID "<<PropelBase::galdef.galdef_ID<<" doing timestep dt="<<dt/year2sec<<" yr  "
                  <<particle.name;
               INFO(tBuf.str());
            }


            for (size_t irept = 1; irept <= size_t(PropelBase::galdef.timestep_repeat); ++irept) {
               // z propagation
               if (1 == PropelBase::galdef.prop_z) {
#pragma omp for schedule(static) 
                  for (size_t ir = 0; ir < size_t(particle.n_rgrid)-1; ++ir) {
                     //There is no need to do the last slice, because the boundary condition forces it to 0
                     for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
                        T* N0 = N0_store+(iz*particle.n_pgrid);
                        T* N1 = N1_store+(iz*particle.n_pgrid);
                        T* N2 = N2_store+(iz*particle.n_pgrid);
                        T* N3 = N3_store+(iz*particle.n_pgrid);
                        T* a1 = alpha1_z+indz(ir,iz,0);
                        T* a2 = alpha2_z+indz(ir,iz,0);
                        T* a3 = alpha3_z+indz(ir,iz,0);
                        T* tsf = total_source_function+indz(ir,iz,0);
                        T* crd = cr_density+indz(ir,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N1, N2, N3, N0, a1, a2, a3, tsf, crd:64)
#endif
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           N1[ip] = -0.5*a1[ip]*dt;
                           N2[ip] = 1. + 0.5*a2[ip]*dt;
                           N3[ip] = -0.5*a3[ip]*dt;

                           N0[ip] = tsf[ip]*f_use*dt + (1. - 0.5*a2[ip]*dt)*crd[ip];
                        }
                     }
                     for (size_t iz = 1; iz < size_t(particle.n_zgrid); ++iz) {
                        T* N0 = N0_store+(iz*particle.n_pgrid);
                        T* a1 = alpha1_z+indz(ir,iz,0);
                        T* crd = cr_density+indz(ir,iz-1,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N0, a1, crd:64)
#endif
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           N0[ip] += 0.5*a1[ip]*dt*crd[ip];
                        }
                     }
                     for (size_t iz = 0; iz < size_t(particle.n_zgrid)-1; ++iz) {
                        T* N0 = N0_store+(iz*particle.n_pgrid);
                        T* a3 = alpha3_z+indz(ir,iz,0);
                        T* crd = cr_density+indz(ir,iz+1,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N0, a3, crd:64)
#endif
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           N0[ip] += 0.5*a3[ip]*dt*crd[ip];
                        }
                     }

                     PropelOperatorSplitting<T>::vectorTridiag(N1_store, N2_store, N3_store, N0_store, R_store, gam_store, bet_store, particle.n_zgrid, particle.n_pgrid);

                     for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
                        T* crd = cr_density+indz(ir,iz,0);
                        T* R = R_store+(iz*particle.n_pgrid);
#if _OPENMP >= 201307
#pragma omp simd aligned(R, crd:64)
#endif
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           crd[ip] = (R[ip] < 0.) ? 0. : R[ip];
                        }
                     }

                     T* crd1 = cr_density+indz(ir,0,0);
                     T* crd2 = cr_density+indz(ir,particle.n_zgrid-1,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(crd1, crd2:64)
#endif
                     for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                        //Boundary condition
                        crd1[ip] = 0.;
                        crd2[ip] = 0.;
                     }

                  } //OMP loop
               } // z propagation

               // momentum propagation
               if (1 == PropelBase::galdef.prop_p) {
#pragma omp for schedule(static) 
                  for (size_t iz = 0; iz < size_t(particle.n_zgrid)-1; ++iz) {
                     //There is no need to do the first and last slices, because the boundary condition forces them to 0
                     if (iz == 0)
                        continue;
                     for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                        T* N0 = N0_store+(ip*particle.n_rgrid);
                        T* N1 = N1_store+(ip*particle.n_rgrid);
                        T* N2 = N2_store+(ip*particle.n_rgrid);
                        T* N3 = N3_store+(ip*particle.n_rgrid);
                        T* a1 = malpha1_p+indp(0,iz,ip);
                        T* a2 = malpha2_p+indp(0,iz,ip);
                        T* a3 = malpha3_p+indp(0,iz,ip);
                        T* tsf = mp_tsf+indp(0,iz,ip);
#if _OPENMP >= 201307
#pragma omp simd aligned(N1, N2, N3, N0, a1, a2, a3, tsf:64)
#endif
                        for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
                           N1[ir] =    - 0.5*a1[ir]*dt;
                           N2[ir] = 1. + 0.5*a2[ir]*dt;
                           N3[ir] =    - 0.5*a3[ir]*dt;

                           N0[ir] = tsf[ir]*f_use*dt;
                        }
                        //CR density is not continuous, so no simd
                        for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
                              N0[ir] += (1. - 0.5*a2[ir]*dt)*cr_density[indz(ir,iz,ip)];
                        }
                     }
                     for (size_t ip = 1; ip < size_t(particle.n_pgrid); ++ip) {
                        T* N0 = N0_store+(ip*particle.n_rgrid);
                        T* a1 = malpha1_p+indp(0,iz,ip);
                        for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
                           N0[ir] += 0.5*a1[ir]*dt*cr_density[indz(ir,iz,ip-1)];
                        }
                     }
                     for (size_t ip = 0; ip < size_t(particle.n_pgrid)-1; ++ip) {
                        T* N0 = N0_store+(ip*particle.n_rgrid);
                        T* a3 = malpha3_p+indp(0,iz,ip);
                        for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
                           N0[ir] += 0.5*a3[ir]*dt*cr_density[indz(ir,iz,ip+1)];
                        }
                     }

                     PropelOperatorSplitting<T>::vectorTridiag(N1_store, N2_store, N3_store, N0_store, R_store, gam_store, bet_store, particle.n_pgrid, particle.n_rgrid);

                     for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                        T* R = R_store+(ip*particle.n_rgrid);
                        for (size_t ir = 0; ir < size_t(particle.n_rgrid)-1; ++ir) {
                           cr_density[indz(ir,iz,ip)] = (R[ir] < 0.) ? 0. : R[ir];
                        }
                     }

                     //No boundary conditions here
                  }
               } // momentum propagation

               // r propagation
               if (1 == PropelBase::galdef.prop_r) {
#pragma omp for schedule(static) 
                  for (size_t iz = 0; iz < size_t(particle.n_zgrid)-1; ++iz) {
                     //There is no need to do the last slice, because the boundary condition forces it to 0
                     if (iz == 0)
                        continue; //We do this with an if statment for NUMA purposes
                     for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
                        T* N0 = N0_store+(ir*particle.n_pgrid);
                        T* N1 = N1_store+(ir*particle.n_pgrid);
                        T* N2 = N2_store+(ir*particle.n_pgrid);
                        T* N3 = N3_store+(ir*particle.n_pgrid);
                        T* a1 = malpha1_r+indr(ir,iz,0);
                        T* a2 = malpha2_r+indr(ir,iz,0);
                        T* a3 = malpha3_r+indr(ir,iz,0);
                        T* tsf = mr_tsf+indr(ir,iz,0);
                        T* crd = cr_density+indz(ir,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N1, N2, N3, N0, a1, a2, a3, tsf, crd:64)
#endif
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           N1[ip] = -0.5*a1[ip]*dt;
                           N2[ip] = 1. + 0.5*a2[ip]*dt;
                           N3[ip] = -0.5*a3[ip]*dt;

                           N0[ip] = tsf[ip]*f_use*dt + (1. - 0.5*a2[ip]*dt)*crd[ip];
                        }
                     }
                     for (size_t ir = 1; ir < size_t(particle.n_rgrid); ++ir) {
                        T* N0 = N0_store+(ir*particle.n_pgrid);
                        T* a1 = malpha1_r+indr(ir,iz,0);
                        T* crd = cr_density+indz(ir-1,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N0, a1, crd:64)
#endif
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           N0[ip] += 0.5*a1[ip]*dt*crd[ip];
                        }
                     }
                     for (size_t ir = 0; ir < size_t(particle.n_rgrid)-1; ++ir) {
                        T* N0 = N0_store+(ir*particle.n_pgrid);
                        T* a3 = malpha3_r+indr(ir,iz,0);
                        T* crd = cr_density+indz(ir+1,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N0, a3, crd:64)
#endif
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           N0[ip] += 0.5*a3[ip]*dt*crd[ip];
                        }
                     }
                     PropelOperatorSplitting<T>::vectorTridiag(N1_store, N2_store, N3_store, N0_store, R_store, gam_store, bet_store, particle.n_rgrid, particle.n_pgrid);

                     for (size_t ir = 0; ir < size_t(particle.n_rgrid)-1; ++ir) {
                        T* R = R_store+(ir*particle.n_pgrid);
                        T* crd = cr_density+indz(ir,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(R, crd:64)
#endif
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           crd[ip] = (R[ip] < 0.) ? 0. : R[ip];
                        }
                     }

                     //Boundary condition
                     T* crd = cr_density+indz(particle.n_rgrid-1,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(crd:64)
#endif
                     for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                        crd[ip] = 0.;
                     }

                  }
               } // r propagation

               //The boundary conditions are properly set above, no need to re-apply them
               /*
                * This doesn't work very well
                bool converged = true;
                for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
                for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
                for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                if (fabs(cr_density[indz(ir,iz,ip)]/cr_density_old[indz(ir,iz,ip)] - 1) > PropelBase::galdef.solution_rel_accuracy) {
                converged = false;
                }
                cr_density_old[indz(ir,iz,ip)] = cr_density[indz(ir,iz,ip)];
                }
                }
                }

                if (converged) {
                std::ostringstream buf;
                buf << "Converged after "<<irept<<" iterations";
                INFO(buf.str());
                break;
                }
                */

            } //repeat

            //Fix dt for the next timestep 

         }//Propagation

         mlc::FreeAligned(N0_store);
         mlc::FreeAligned(N1_store);
         mlc::FreeAligned(N2_store);
         mlc::FreeAligned(N3_store);
         mlc::FreeAligned(R_store);
         mlc::FreeAligned(gam_store);
         mlc::FreeAligned(bet_store);
      }//parallel

      //Assign the cr density
#pragma omp parallel for default(shared) schedule(static)
      for (size_t ir = 0; ir < size_t(particle.n_rgrid); ++ir) {
         for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
            T* crd = cr_density+indz(ir,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(crd:64)
#endif
            for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
               particle.cr_density.d2[ir][iz].s[ip] = crd[ip];
            }
         }
      }
   } //2D

   else if (3 == particle.n_spatial_dimensions) { // ==== 3D ====

      //Needed by the method
      const T f_use = 1./(PropelBase::galdef.prop_x + PropelBase::galdef.prop_y + PropelBase::galdef.prop_z + PropelBase::galdef.prop_p);

      //Create the index functions, indz is the one corresponding to the base class
      const typename PropelOperatorSplitting<T>::IndexConverter3D indz(particle.n_xgrid, particle.n_ygrid, particle.n_zgrid, particle.n_pgrid);
      const IndexConverter3Dp indp(particle.n_xgrid, particle.n_ygrid, particle.n_zgrid, particle.n_pgrid);
      const IndexConverter3Dy indy(particle.n_xgrid, particle.n_ygrid, particle.n_zgrid, particle.n_pgrid);
      const IndexConverter3Dx indx(particle.n_xgrid, particle.n_ygrid, particle.n_zgrid, particle.n_pgrid);

      //populate the total source function and cr density arrays
#pragma omp parallel for default(shared) schedule(static)
      for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
         for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
            for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
               T* tsf = total_source_function+indz(ix,iy,iz,0);
               T* crd = cr_density+indz(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(tsf, crd:64)
#endif
               for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                  tsf[ip] = T(particle.primary_source_function.d3[ix][iy][iz].s[ip]+particle.secondary_source_function.d3[ix][iy][iz].s[ip]);
                  crd[ip] = T(particle.cr_density.d3[ix][iy][iz].s[ip]);
               }
            }
         }
      }

      //Need to look into these timestep modes Andy introduced
      //Ignore them for now
      
      //Storage and temporary reordering of alpha arrays and source function
      //Only need to store for x,y, and z propagation 
      //This could be avoided by overloading initilizeAlpha, but then we'd
      //need some mechanism for not copying the code because that is difficult to maintain properly

      //Do the reordering using NUMA aware looping
#pragma omp parallel for default(shared) schedule(static) 
      for (size_t iy = 1; iy < size_t(particle.n_ygrid)-1; ++iy) {
         for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
            for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
               for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
                  malpha1_p[indp(ix,iy,iz,ip)] = alpha1_p[indz(ix,iy,iz,ip)];
                  malpha2_p[indp(ix,iy,iz,ip)] = alpha2_p[indz(ix,iy,iz,ip)];
                  malpha3_p[indp(ix,iy,iz,ip)] = alpha3_p[indz(ix,iy,iz,ip)];
                  mp_tsf[indp(ix,iy,iz,ip)] = total_source_function[indz(ix,iy,iz,ip)];
               }
            }
         }
      }
#pragma omp parallel for default(shared) schedule(static)
      for (size_t ix = 1; ix < size_t(particle.n_xgrid)-1; ++ix) {
         for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
            for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
               for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                  malpha1_y[indy(ix,iy,iz,ip)] = alpha1_y[indz(ix,iy,iz,ip)];
                  malpha2_y[indy(ix,iy,iz,ip)] = alpha2_y[indz(ix,iy,iz,ip)];
                  malpha3_y[indy(ix,iy,iz,ip)] = alpha3_y[indz(ix,iy,iz,ip)];
                  my_tsf[indy(ix,iy,iz,ip)] = total_source_function[indz(ix,iy,iz,ip)];
               }
            }
         }
      }

#pragma omp parallel for default(shared) schedule(static)
      for (size_t iy = 1; iy < size_t(particle.n_ygrid)-1; ++iy) {
         for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
            for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
               for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                  malpha1_x[indx(ix,iy,iz,ip)] = alpha1_x[indz(ix,iy,iz,ip)];
                  malpha2_x[indx(ix,iy,iz,ip)] = alpha2_x[indz(ix,iy,iz,ip)];
                  malpha3_x[indx(ix,iy,iz,ip)] = alpha3_x[indz(ix,iy,iz,ip)];
                  mx_tsf[indx(ix,iy,iz,ip)] = total_source_function[indz(ix,iy,iz,ip)];
               }
            }
         }
      }

      //Sort the grid size for allocation below
      std::vector<size_t> sortedGridSizes({size_t(particle.n_xgrid),size_t(particle.n_ygrid), size_t(particle.n_zgrid), size_t(particle.n_pgrid)});
      std::sort(sortedGridSizes.begin(), sortedGridSizes.end());

      //Have the entire loop a big parallel region.  Almost anything is run in parallel within anyway, reduces the overhead
#pragma omp parallel default(shared) 
      {
         //Storage for tridag, use a single one just big enough for all variations
         T* N0_store = static_cast<T*>(mlc::MallocAligned(sortedGridSizes[3]*sortedGridSizes[2]*sizeof(T)));
         T* N1_store = static_cast<T*>(mlc::MallocAligned(sortedGridSizes[3]*sortedGridSizes[2]*sizeof(T)));
         T* N2_store = static_cast<T*>(mlc::MallocAligned(sortedGridSizes[3]*sortedGridSizes[2]*sizeof(T)));
         T* N3_store = static_cast<T*>(mlc::MallocAligned(sortedGridSizes[3]*sortedGridSizes[2]*sizeof(T)));
         T* R_store = static_cast<T*>(mlc::MallocAligned(sortedGridSizes[3]*sortedGridSizes[2]*sizeof(T)));
         T* gam_store = static_cast<T*>(mlc::MallocAligned(sortedGridSizes[3]*sortedGridSizes[2]*sizeof(T)));
         T* bet_store = static_cast<T*>(mlc::MallocAligned(sortedGridSizes[3]*sizeof(T)));

         for (T dt = PropelBase::galdef.start_timestep*year2sec; dt > end_timestep_sec*0.9999; dt*= timestepFactor) {


#pragma omp single
            {
               std::ostringstream tBuf;
               tBuf<<" galdef_ID "<<PropelBase::galdef.galdef_ID<<" doing timestep dt="<<dt/year2sec<<" yr  "
                  <<particle.name;
               INFO(tBuf.str());
            }

            for (size_t irept = 1; irept <= size_t(PropelBase::galdef.timestep_repeat); ++irept) {
               // x propagation
               if (1 == PropelBase::galdef.prop_x) {
#pragma omp for schedule(static) 
                  for (size_t iy = 1; iy < size_t(particle.n_ygrid)-1; ++iy) {
                     for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
                        //There is no need to do the last slice, because the boundary condition forces it to 0
                        for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
                           T* N0 = N0_store+(ix*particle.n_pgrid);
                           T* N1 = N1_store+(ix*particle.n_pgrid);
                           T* N2 = N2_store+(ix*particle.n_pgrid);
                           T* N3 = N3_store+(ix*particle.n_pgrid);
                           T* a1 = malpha1_x+indx(ix,iy,iz,0);
                           T* a2 = malpha2_x+indx(ix,iy,iz,0);
                           T* a3 = malpha3_x+indx(ix,iy,iz,0);
                           T* tsf = mx_tsf+indx(ix,iy,iz,0);
                           T* crd = cr_density+indz(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N1, N2, N3, N0, a1, a2, a3, tsf, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                              N1[ip] =    - 0.5*a1[ip]*dt;
                              N2[ip] = 1. + 0.5*a2[ip]*dt;
                              N3[ip] =    - 0.5*a3[ip]*dt;

                              N0[ip] = tsf[ip]*f_use*dt + (1. - 0.5*a2[ip]*dt)*crd[ip];
                           }
                        }
                        for (size_t ix = 1; ix < size_t(particle.n_xgrid); ++ix) {
                           T* N0 = N0_store+(ix*particle.n_pgrid);
                           T* a1 = malpha1_x+indx(ix,iy,iz,0);
                           T* crd = cr_density+indz(ix-1,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N0, a1, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                              N0[ip] += 0.5*a1[ip]*dt*crd[ip];
                           }
                        }
                        for (size_t ix = 0; ix < size_t(particle.n_xgrid)-1; ++ix) {
                           T* N0 = N0_store+(ix*particle.n_pgrid);
                           T* a3 = malpha3_x+indx(ix,iy,iz,0);
                           T* crd = cr_density+indz(ix+1,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N0, a3, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                              N0[ip] += 0.5*a3[ip]*dt*crd[ip];
                           }
                        }

                        PropelOperatorSplitting<T>::vectorTridiag(N1_store, N2_store, N3_store, N0_store, R_store, gam_store, bet_store, particle.n_xgrid, particle.n_pgrid);

                        for (size_t ix = 1; ix < size_t(particle.n_xgrid)-1; ++ix) {
                           T* R = R_store+(ix*particle.n_pgrid);
                           T* crd = cr_density+indz(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(R, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                              crd[ip] = (R[ip] < 0.) ? 0. : R[ip];
                           }
                        }

                        //Boundary condition
                        T* crd1 = cr_density+indz(0,iy,iz,0);
                        T* crd2 = cr_density+indz(particle.n_xgrid-1,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(crd1, crd2:64)
#endif
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           crd1[ip] = 0.;
                           crd2[ip] = 0.;
                        }

                     }
                  }
               } // x propagation

               // y propagation
               if (1 == PropelBase::galdef.prop_y) {
#pragma omp for schedule(static)
                  for (size_t ix = 1; ix < size_t(particle.n_xgrid)-1; ++ix) {
                     for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
                        //There is no need to do the last slice, because the boundary condition forces it to 0
                        for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
                           T* N0 = N0_store+(iy*particle.n_pgrid);
                           T* N1 = N1_store+(iy*particle.n_pgrid);
                           T* N2 = N2_store+(iy*particle.n_pgrid);
                           T* N3 = N3_store+(iy*particle.n_pgrid);
                           T* a1 = malpha1_y+indy(ix,iy,iz,0);
                           T* a2 = malpha2_y+indy(ix,iy,iz,0);
                           T* a3 = malpha3_y+indy(ix,iy,iz,0);
                           T* tsf = my_tsf+indy(ix,iy,iz,0);
                           T* crd = cr_density+indz(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N1, N2, N3, N0, a1, a2, a3, tsf, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                              N1[ip] =    - 0.5*a1[ip]*dt;
                              N2[ip] = 1. + 0.5*a2[ip]*dt;
                              N3[ip] =    - 0.5*a3[ip]*dt;

                              N0[ip] = tsf[ip]*f_use*dt + (1. - 0.5*a2[ip]*dt)*crd[ip];
                           }
                        }
                        for (size_t iy = 1; iy < size_t(particle.n_ygrid); ++iy) {
                           T* N0 = N0_store+(iy*particle.n_pgrid);
                           T* a1 = malpha1_y+indy(ix,iy,iz,0);
                           T* crd = cr_density+indz(ix,iy-1,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N0, a1, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                              N0[ip] += 0.5*a1[ip]*dt*crd[ip];
                           }
                        }
                        for (size_t iy = 0; iy < size_t(particle.n_ygrid)-1; ++iy) {
                           T* N0 = N0_store+(iy*particle.n_pgrid);
                           T* a3 = malpha3_y+indy(ix,iy,iz,0);
                           T* crd = cr_density+indz(ix,iy+1,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N0, a3, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                              N0[ip] += 0.5*a3[ip]*dt*crd[ip];
                           }
                        }

                        PropelOperatorSplitting<T>::vectorTridiag(N1_store, N2_store, N3_store, N0_store, R_store, gam_store, bet_store, particle.n_ygrid, particle.n_pgrid);

                        for (size_t iy = 1; iy < size_t(particle.n_ygrid)-1; ++iy) {
                           T* R = R_store+(iy*particle.n_pgrid);
                           T* crd = cr_density+indz(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(R, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                              crd[ip] = (R[ip] < 0.) ? 0. : R[ip];
                           }
                        }

                        //Boundary condition
                        T* crd1 = cr_density+indz(ix,0,iz,0);
                        T* crd2 = cr_density+indz(ix,particle.n_ygrid-1,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(crd1, crd2:64)
#endif
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           crd1[ip] = 0.;
                           crd2[ip] = 0.;
                        }

                     }
                  }
               } // y propagation

               // z propagation
               if (1 == PropelBase::galdef.prop_z) {
#pragma omp for schedule(static)
                  for (size_t ix = 0; ix < size_t(particle.n_xgrid)-1; ++ix) {
                     if (ix == 0)
                        continue; //For NUMA, no need to do first slice
                     for (size_t iy = 1; iy < size_t(particle.n_ygrid)-1; ++iy) {
                        //There is no need to do the first and last slice, because the boundary condition forces it to 0
                        for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
                           T* N0 = N0_store+(iz*particle.n_pgrid);
                           T* N1 = N1_store+(iz*particle.n_pgrid);
                           T* N2 = N2_store+(iz*particle.n_pgrid);
                           T* N3 = N3_store+(iz*particle.n_pgrid);
                           T* a1 = alpha1_z+indz(ix,iy,iz,0);
                           T* a2 = alpha2_z+indz(ix,iy,iz,0);
                           T* a3 = alpha3_z+indz(ix,iy,iz,0);
                           T* tsf = total_source_function+indz(ix,iy,iz,0);
                           T* crd = cr_density+indz(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N1, N2, N3, N0, a1, a2, a3, tsf, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                              N1[ip] =    - 0.5*a1[ip]*dt;
                              N2[ip] = 1. + 0.5*a2[ip]*dt;
                              N3[ip] =    - 0.5*a3[ip]*dt;

                              N0[ip] = tsf[ip]*f_use*dt + (1. - 0.5*a2[ip]*dt)*crd[ip];
                           }
                        }
                        for (size_t iz = 1; iz < size_t(particle.n_zgrid); ++iz) {
                           T* N0 = N0_store+(iz*particle.n_pgrid);
                           T* a1 = alpha1_z+indz(ix,iy,iz,0);
                           T* crd = cr_density+indz(ix,iy,iz-1,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N0, a1, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                              N0[ip] += 0.5*a1[ip]*dt*crd[ip];
                           }
                        }
                        for (size_t iz = 0; iz < size_t(particle.n_zgrid)-1; ++iz) {
                           T* N0 = N0_store+(iz*particle.n_pgrid);
                           T* a3 = alpha3_z+indz(ix,iy,iz,0);
                           T* crd = cr_density+indz(ix,iy,iz+1,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(N0, a3, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                              N0[ip] += 0.5*a3[ip]*dt*crd[ip];
                           }
                        }

                        PropelOperatorSplitting<T>::vectorTridiag(N1_store, N2_store, N3_store, N0_store, R_store, gam_store, bet_store, particle.n_zgrid, particle.n_pgrid);

                        for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
                           T* R = R_store+(iz*particle.n_pgrid);
                           T* crd = cr_density+indz(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(R, crd:64)
#endif
                           for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                              crd[ip] = (R[ip] < 0.) ? 0. : R[ip];
                           }
                        }

                        //Boundary condition
                        T* crd1 = cr_density+indz(ix,iy,0,0);
                        T* crd2 = cr_density+indz(ix,iy,particle.n_zgrid-1,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(crd1, crd2:64)
#endif
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           crd1[ip] = 0.;
                           crd2[ip] = 0.;
                        }

                     }
                  }
               } // z propagation

               // momentum propagation
               if (1 == PropelBase::galdef.prop_p) {
#pragma omp for schedule(static)
                  for (size_t iy = 1; iy < size_t(particle.n_ygrid)-1; ++iy) {
                     for (size_t iz = 1; iz < size_t(particle.n_zgrid)-1; ++iz) {
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           T* N0 = N0_store+(ip*particle.n_xgrid);
                           T* N1 = N1_store+(ip*particle.n_xgrid);
                           T* N2 = N2_store+(ip*particle.n_xgrid);
                           T* N3 = N3_store+(ip*particle.n_xgrid);
                           T* a1 = malpha1_p+indp(0,iy,iz,ip);
                           T* a2 = malpha2_p+indp(0,iy,iz,ip);
                           T* a3 = malpha3_p+indp(0,iy,iz,ip);
                           T* tsf = mp_tsf+indp(0,iy,iz,ip);
#if _OPENMP >= 201307
#pragma omp simd aligned(N1, N2, N3, N0, a1, a2, a3, tsf:64)
#endif
                           for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
                              N1[ix] =    - 0.5*a1[ix]*dt;
                              N2[ix] = 1. + 0.5*a2[ix]*dt;
                              N3[ix] =    - 0.5*a3[ix]*dt;

                              N0[ix] = tsf[ix]*f_use*dt;
                           }
                           for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
                              N0[ix] += (1. - 0.5*malpha2_p[indp(ix,iy,iz,ip)]*dt)*cr_density[indz(ix,iy,iz,ip)];
                           }
                        }
                        for (size_t ip = 1; ip < size_t(particle.n_pgrid); ++ip) {
                           T* N0 = N0_store+(ip*particle.n_xgrid);
                           for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
                              N0[ix] += 0.5*malpha1_p[indp(ix,iy,iz,ip)]*dt*cr_density[indz(ix,iy,iz,ip-1)];
                           }
                        }
                        for (size_t ip = 0; ip < size_t(particle.n_pgrid)-1; ++ip) {
                           T* N0 = N0_store+(ip*particle.n_xgrid);
                           for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
                              N0[ix] += 0.5*malpha3_p[indp(ix,iy,iz,ip)]*dt*cr_density[indz(ix,iy,iz,ip+1)];
                           }
                        }

                        PropelOperatorSplitting<T>::vectorTridiag(N1_store, N2_store, N3_store, N0_store, R_store, gam_store, bet_store, particle.n_pgrid, particle.n_xgrid);

                        for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                           T* R = R_store+(ip*particle.n_xgrid);
                           for (size_t ix = 1; ix < size_t(particle.n_xgrid)-1; ++ix) {
                              cr_density[indz(ix,iy,iz,ip)] = (R[ix] < 0.) ? 0. : R[ix];
                           }
                        }

                        //No boundary conditions here
                     }
                  }
               } // momentum propagation


               //The boundary conditions are properly set above, no need to re-apply them
               
               /*
                * This is not very useful and breaks in the faster big parallel region for some reason
                *
               //Do convergence criterion at each step
               bool converged = true;
               for (size_t i = 0; i < cr_density.size(); ++i) {
                  if (fabs(cr_density[i]/cr_density_old[i] - 1) > PropelBase::galdef.solution_rel_accuracy) {
                     converged = false;
                  }
                  cr_density_old[i] = cr_density[i];
               }

               if (converged) {
#pragma omp single
                  {
                     std::ostringstream buf;
                     buf << "Converged after "<<irept<<" iterations";
                     INFO(buf.str());
                  }
                  break;
               }
               */

            } //repeat

         }//Propagation
         mlc::FreeAligned(N0_store);
         mlc::FreeAligned(N1_store);
         mlc::FreeAligned(N2_store);
         mlc::FreeAligned(N3_store);
         mlc::FreeAligned(R_store);
         mlc::FreeAligned(gam_store);
         mlc::FreeAligned(bet_store);
      } // Parallel


      //Assign the cr density
#pragma omp parallel for default(shared) schedule(static)
      for (size_t ix = 0; ix < size_t(particle.n_xgrid); ++ix) {
         for (size_t iy = 0; iy < size_t(particle.n_ygrid); ++iy) {
            for (size_t iz = 0; iz < size_t(particle.n_zgrid); ++iz) {
               T* crd = cr_density+indz(ix,iy,iz,0);
#if _OPENMP >= 201307
#pragma omp simd aligned(crd:64)
#endif
               for (size_t ip = 0; ip < size_t(particle.n_pgrid); ++ip) {
                  particle.cr_density.d3[ix][iy][iz].s[ip] = crd[ip];
               }
            }
         }
      }
   } //3D
}

//Explicit instantiation
template class PropelCrankNicolsonVector<float>;
template class PropelCrankNicolsonVector<double>;
