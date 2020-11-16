/*
 * This is the old propel routine using the new PropelBase class
 *
 * It is now separated into three classes depending on the method 
 * and turbo settings.
 *
 * This separation should make it easier to modify and optimize parts
 * of the code and try new stuff.
 *
 * We template the class to have both a double and float versions.
 *
 * There is a single base class that handles the calculation of 
 * the alphas.
 */

#include "PropelBase.h"

#include <vector>

template <typename T>
class PropelOperatorSplitting : public PropelBase {
   protected:
      //Store the alphas as a single arry, it should be a more efficient storage.
      //Handle memory allocations ourselves and use aligned allocator from boost.
      //The indexing scheme is xyzp in 3D and rzp in 2D, from slowest to fastest.  p is always fastest.
      //It is assumed that temporary storage is created to re-order the elements more favorably.
      T *alpha1_p, *alpha1_z, *alpha1_r, *alpha1_x, *alpha1_y;
      T *alpha2_p, *alpha2_z, *alpha2_r, *alpha2_x, *alpha2_y;
      T *alpha3_p, *alpha3_z, *alpha3_r, *alpha3_x, *alpha3_y;

      //This will return false if nothing needs to be done
      bool initializeAlpha(Particle &particle);

      //Re-writing tridag in c++ with modern storage
      //This is mostly for the templating benefits for automatic selection of
      //double and float.
      static void tridiag(const T* a,
            const T* b,
            const T* c, //a, b, and c is the matrix
            const T* r, //The results.
            T* u, //The solution
            T* gam, //Temporary array
            const size_t n //Number of equations
            );

      //Vectorized version of tridag that solves many equations at once.
      //Should enable simd vectorization of the solution
      //a, b, c, r, u, gam should have size stride*n were stride is the 
      //number of systems and n is the size of each system.
      //bet should have size stride;
      static void vectorTridiag(const T* a,
            const T* b,
            const T* c, //a, b, and c is the matrix
            const T* r, //The results.
            T* u, //The solution
            T* gam, //Temporary array
            T* bet, //Temporary array
            const size_t n, //number of equations of each system
            const size_t stride //number of systems
            );

      //Class to convert ir, iz, and ip indices to compressed irzp index
      class IndexConverter2D {
         private:
            const size_t nr, nz, np;
         public:
            IndexConverter2D(size_t Nr, size_t Nz, size_t Np) :
               nr(Nr), nz(Nz), np(Np)
            {}
            inline size_t operator() (size_t ir, size_t iz, size_t ip) const 
            { return ip + np*(iz + nz*ir); }
      };

      //Class to convert ix, iy, iz, and ip indices to compressed ixyzp index
      class IndexConverter3D {
         private:
            const size_t nx, ny, nz, np;
         public:
            IndexConverter3D(size_t Nx, size_t Ny, size_t Nz, size_t Np) :
               nx(Nx), ny(Ny), nz(Nz), np(Np)
            {}
            inline size_t operator() (size_t ix, size_t iy, size_t iz, size_t ip) const 
            { return ip + np*(iz + nz*(iy + ny*ix)); }
         };

   public:
      PropelOperatorSplitting(const Particle &particle, const Galdef &_galdef);

      virtual ~PropelOperatorSplitting();

      virtual void operator() (Particle &particle) override = 0;
};


template <typename T>
class PropelCrankNicolson : public PropelOperatorSplitting<T> {
   protected:
      T* total_source_function;
      T* cr_density;
      //Storage for re-ordering of alphas to make things faster in the loops
      //Do not need p, they are correct order
      T* malpha1_z, *malpha1_r, *malpha1_x, *malpha1_y;
      T* malpha2_z, *malpha2_r, *malpha2_x, *malpha2_y;
      T* malpha3_z, *malpha3_r, *malpha3_x, *malpha3_y;
      //Same for the total source function 
      T* mz_tsf, *mr_tsf, *mx_tsf, *my_tsf;
      
      //Classes for index conversion in r and z cases
      class IndexConverter2Dr {
         private:
            const size_t nr, nz, np;
         public:
            IndexConverter2Dr(size_t Nr, size_t Nz, size_t Np) :
               nr(Nr), nz(Nz), np(Np)
            {}
            inline size_t operator() (size_t ir, size_t iz, size_t ip) const 
            { return ir + nr*(ip + np*iz); }
      };
      class IndexConverter2Dz {
         private:
            const size_t nr, nz, np;
         public:
            IndexConverter2Dz(size_t Nr, size_t Nz, size_t Np) :
               nr(Nr), nz(Nz), np(Np)
            {}
            inline size_t operator() (size_t ir, size_t iz, size_t ip) const 
            { return iz + nz*(ip + np*ir); }
      };
      //Classes for index conversion in x,y, and z cases
      class IndexConverter3Dx {
         private:
            const size_t nx, ny, nz, np;
         public:
            IndexConverter3Dx(size_t Nx, size_t Ny, size_t Nz, size_t Np) :
               nx(Nx), ny(Ny), nz(Nz), np(Np)
            {}
            inline size_t operator() (size_t ix, size_t iy, size_t iz, size_t ip) const 
            { return ix + nx*(ip + np*(iz + nz*iy)); }
      };
      class IndexConverter3Dy {
         private:
            const size_t nx, ny, nz, np;
         public:
            IndexConverter3Dy(size_t Nx, size_t Ny, size_t Nz, size_t Np) :
               nx(Nx), ny(Ny), nz(Nz), np(Np)
            {}
            inline size_t operator() (size_t ix, size_t iy, size_t iz, size_t ip) const 
            { return iy + ny*(ip + np*(iz + nz*ix)); }
      };
      class IndexConverter3Dz {
         private:
            const size_t nx, ny, nz, np;
         public:
            IndexConverter3Dz(size_t Nx, size_t Ny, size_t Nz, size_t Np) :
               nx(Nx), ny(Ny), nz(Nz), np(Np)
            {}
            inline size_t operator() (size_t ix, size_t iy, size_t iz, size_t ip) const 
            { return iz + nz*(ip + np*(iy + ny*ix)); }
      };

      using PropelOperatorSplitting<T>::alpha1_p;
      using PropelOperatorSplitting<T>::alpha2_p;
      using PropelOperatorSplitting<T>::alpha3_p;
      using PropelOperatorSplitting<T>::alpha1_r;
      using PropelOperatorSplitting<T>::alpha2_r;
      using PropelOperatorSplitting<T>::alpha3_r;
      using PropelOperatorSplitting<T>::alpha1_z;
      using PropelOperatorSplitting<T>::alpha2_z;
      using PropelOperatorSplitting<T>::alpha3_z;
      using PropelOperatorSplitting<T>::alpha1_x;
      using PropelOperatorSplitting<T>::alpha2_x;
      using PropelOperatorSplitting<T>::alpha3_x;
      using PropelOperatorSplitting<T>::alpha1_y;
      using PropelOperatorSplitting<T>::alpha2_y;
      using PropelOperatorSplitting<T>::alpha3_y;

   public:
      PropelCrankNicolson(const Particle &_particle, const Galdef &_galdef);

      virtual ~PropelCrankNicolson();

      virtual void operator() (Particle &particle) override;
};

using PropelCrankNicolsonF = PropelCrankNicolson<float>;
using PropelCrankNicolsonD = PropelCrankNicolson<double>;

template <typename T>
class PropelCrankNicolsonVector : public PropelOperatorSplitting<T> {
   protected:
      T* total_source_function;
      T* cr_density;
      //Storage for re-ordering of alphas to make things faster in the loops
      //Do not need z, they are correct order
      T* malpha1_p, *malpha1_r, *malpha1_x, *malpha1_y;
      T* malpha2_p, *malpha2_r, *malpha2_x, *malpha2_y;
      T* malpha3_p, *malpha3_r, *malpha3_x, *malpha3_y;
      //Same for the total source function 
      T* mp_tsf, *mr_tsf, *mx_tsf, *my_tsf;
      
      //Classes for index conversion in r and p cases
      class IndexConverter2Dr {
         private:
            const size_t nr, nz, np;
         public:
            IndexConverter2Dr(size_t Nr, size_t Nz, size_t Np) :
               nr(Nr), nz(Nz), np(Np)
            {}
            inline size_t operator() (size_t ir, size_t iz, size_t ip) const 
            { return ip + np*(ir + nr*iz); }
      };
      class IndexConverter2Dp {
         private:
            const size_t nr, nz, np;
         public:
            IndexConverter2Dp(size_t Nr, size_t Nz, size_t Np) :
               nr(Nr), nz(Nz), np(Np)
            {}
            inline size_t operator() (size_t ir, size_t iz, size_t ip) const 
            { return ir + nr*(ip + np*iz); }
      };
      //Classes for index conversion in x,y, and z cases
      class IndexConverter3Dx {
         private:
            const size_t nx, ny, nz, np;
         public:
            IndexConverter3Dx(size_t Nx, size_t Ny, size_t Nz, size_t Np) :
               nx(Nx), ny(Ny), nz(Nz), np(Np)
            {}
            inline size_t operator() (size_t ix, size_t iy, size_t iz, size_t ip) const 
            { return ip + np*(ix + nx*(iz + nz*iy)); }
      };
      class IndexConverter3Dy {
         private:
            const size_t nx, ny, nz, np;
         public:
            IndexConverter3Dy(size_t Nx, size_t Ny, size_t Nz, size_t Np) :
               nx(Nx), ny(Ny), nz(Nz), np(Np)
            {}
            inline size_t operator() (size_t ix, size_t iy, size_t iz, size_t ip) const 
            { return ip + np*(iy + ny*(iz + nz*ix)); }
      };
      class IndexConverter3Dp {
         private:
            const size_t nx, ny, nz, np;
         public:
            IndexConverter3Dp(size_t Nx, size_t Ny, size_t Nz, size_t Np) :
               nx(Nx), ny(Ny), nz(Nz), np(Np)
            {}
            inline size_t operator() (size_t ix, size_t iy, size_t iz, size_t ip) const 
            { return ix + nx*(ip + np*(iz + nz*iy)); }
      };

      using PropelOperatorSplitting<T>::alpha1_p;
      using PropelOperatorSplitting<T>::alpha2_p;
      using PropelOperatorSplitting<T>::alpha3_p;
      using PropelOperatorSplitting<T>::alpha1_r;
      using PropelOperatorSplitting<T>::alpha2_r;
      using PropelOperatorSplitting<T>::alpha3_r;
      using PropelOperatorSplitting<T>::alpha1_z;
      using PropelOperatorSplitting<T>::alpha2_z;
      using PropelOperatorSplitting<T>::alpha3_z;
      using PropelOperatorSplitting<T>::alpha1_x;
      using PropelOperatorSplitting<T>::alpha2_x;
      using PropelOperatorSplitting<T>::alpha3_x;
      using PropelOperatorSplitting<T>::alpha1_y;
      using PropelOperatorSplitting<T>::alpha2_y;
      using PropelOperatorSplitting<T>::alpha3_y;

   public:
      PropelCrankNicolsonVector(const Particle &_particle, const Galdef &_galdef);

      virtual ~PropelCrankNicolsonVector();

      virtual void operator() (Particle &particle) override;
};

