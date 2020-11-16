#ifndef FULL_SKY_H
#define FULL_SKY_H

#include "BaseSky.h"

#include <tuple>
#include <type_traits>

namespace SM {


/** Skymap that stores all pixels
 *
 * Implements the storage as a single vector.
 * Stores it packed such that each spectrum is continuous
 *
 * TODO: Have an option to pack the spectra map first.
 */

template <typename T>
class FullSky : public BaseSky<T>
{
   protected:
      std::vector<T> storage;

      /** \brief Iterator implementation that covers the entire sky. 
       */
      template < typename C >
      class FullSky_IteratorImplementation : public SM::BaseSky<T>::template IteratorImplementation<C> {
         protected:
            //Declare the iterator const_iterator if C is const
            //Basically use this for accessing the memory
            using iterator_type = typename std::conditional<
                  std::is_const<C>::value, 
                  typename std::vector<T>::const_iterator, 
                  typename std::vector<T>::iterator>::type;
            iterator_type it;
            using const_iterator_type = typename std::vector<T>::const_iterator;

            //Need a constant reference to the beginning
            const_iterator_type itBegin;

            //We must know the size of the spectrum
            size_t enSize;

            //Keep track of where we are
            size_t eIndex, hpIndex;

         public:
            FullSky_IteratorImplementation( iterator_type &&_it, 
                  const_iterator_type &&begin,
                  size_t enS ) 
               : it(std::forward<iterator_type>(_it)), itBegin(std::forward<const_iterator_type>(begin)), 
               enSize(enS) 
            {
               size_t index = it - itBegin;
               eIndex = index%enSize;
               hpIndex = index/enSize;
            }

            virtual void increment() override 
            { 
               ++it; 
               ++eIndex;
               if (eIndex == enSize) {
                  eIndex = 0;
                  ++hpIndex;
               }
            }
            virtual void decrement() override 
            { 
               --it; 
               if (eIndex == 0) {
                  eIndex = enSize-1;
                  --hpIndex;
               } else {
                  --eIndex;
               }
            }

            virtual bool equal(const typename SM::BaseSky<T>::template IteratorImplementation<C>& other) const override { 
               return it == dynamic_cast<const FullSky_IteratorImplementation<C>& >(other).it; 
            }
            
            virtual C& dereference() const override {
               return *it;
            }

            virtual size_t healpixIndex() const override {
               return hpIndex;
            }

            virtual size_t energyIndex() const override {
               return eIndex;
            }

            virtual FullSky_IteratorImplementation<C>* clone() const override {
               //Need new copies of the iterators and the size of the spectrum
               return new FullSky_IteratorImplementation<C>( iterator_type(it), const_iterator_type(itBegin), enSize );
            }

      }; 

      /** \brief Iterator for stepping over maps in a single energy bin
       */
      template< typename C >
      class Map_IteratorImplementation : public FullSky_IteratorImplementation<C> {
         protected:
            using iterator_type = typename FullSky_IteratorImplementation<C>::iterator_type;
            using const_iterator_type = typename FullSky_IteratorImplementation<C>::const_iterator_type;
            using FullSky_IteratorImplementation<C>::it;
            using FullSky_IteratorImplementation<C>::enSize;
            using FullSky_IteratorImplementation<C>::hpIndex;
            
         public:
            Map_IteratorImplementation( iterator_type &&_it, 
                  const_iterator_type &&begin,
                  size_t enS ) 
               : FullSky_IteratorImplementation<C>(
                     std::forward<iterator_type>(_it),
                     std::forward<const_iterator_type>(begin),
                    enS) {} 

            virtual void increment() override { 
               it += enSize; 
               ++hpIndex;
            }
            virtual void decrement() override { 
               it -= enSize; 
               --hpIndex;
            }
      };

      /** \brief Iterator for stepping over spectra in a single pixel
       */
      template< typename C >
      class Spectra_IteratorImplementation : public FullSky_IteratorImplementation<C> {
         protected:
            using iterator_type = typename FullSky_IteratorImplementation<C>::iterator_type;
            using const_iterator_type = typename FullSky_IteratorImplementation<C>::const_iterator_type;
            using FullSky_IteratorImplementation<C>::it;
            using FullSky_IteratorImplementation<C>::eIndex;

         public:
            Spectra_IteratorImplementation( iterator_type &&_it, 
                  const_iterator_type &&begin,
                  size_t enS ) 
               : FullSky_IteratorImplementation<C>(
                     std::forward<iterator_type>(_it),
                     std::forward<const_iterator_type>(begin),
                    enS) {} 

            virtual void increment() override { 
               ++it; 
               ++eIndex;
            }
            virtual void decrement() override { 
               --it; 
               --eIndex;
            }
      };

      using BaseSky<T>::spBin;
      using BaseSky<T>::empty;

      void setStorage()
      {
         storage.resize(HealpixBaseExtended::Npix()*spBin.GetSize());

         //Assign values with openmp for efficient NUMA access.
         //Use static so we always have the same CPU allocation
#if (_OPENMP >= 201307)
#pragma omp parallel for simd schedule(static) default(shared) 
#else
#pragma omp parallel for schedule(static) default(shared) 
#endif
         for ( size_t i=0; i < storage.size(); ++i )
            storage[i] = empty;
      }

      static std::string getName() { return "FullSky"; }

   public:
      using iterator = typename BaseSky<T>::iterator;
      using const_iterator = typename BaseSky<T>::const_iterator;

      using HealpixBaseExtended::Npix;
      using HealpixBaseExtended::Order;
      using HealpixBaseExtended::Scheme;
      using BaseSky<T>::GetBinning;
      using BaseSky<T>::GetCoordinateSystem;
      using BaseSky<T>::GetHealpixOrder;
      using BaseSky<T>::GetEmptyValue;
      using BaseSky<T>::SetEmptyValue;
      using BaseSky<T>::GetSpectrum;
      using BaseSky<T>::GetValue;
      using BaseSky<T>::SetValue;
      using BaseSky<T>::GetReference;
      using BaseSky<T>::operator[];
      using BaseSky<T>::GetSkyInterpolatedSpectra;
      using BaseSky<T>::Equivalent;

      FullSky(SpectralBinning&& sp, int order, Healpix_Ordering_Scheme scheme, CoordSys cs, const T& null = T(0)) :
         BaseSky<T>(std::forward<SpectralBinning>(sp), order, scheme, cs, null)
      { setStorage(); }
      FullSky(const SpectralBinning& sp, int order, Healpix_Ordering_Scheme scheme, CoordSys cs, const T& null = T(0)) :
         BaseSky<T>(sp, order, scheme, cs, null)
      { setStorage(); }
      FullSky(const BaseSkyStorageIndependent &templ, const T& null) :
         BaseSky<T>(templ, null)
      { setStorage(); }

      virtual std::unique_ptr< BaseSky<T> > clone() const override {
         return std::unique_ptr< BaseSky<T> >(new FullSky<T>(*this));
      }

      virtual void Reset(const SpectralBinning& sp, int order, Healpix_Ordering_Scheme scheme, CoordSys cs, const T& null = T(0)) override
      {
         Healpix_Base::Set(order, scheme);
         BaseSkyStorageIndependent::spBin = sp;
         BaseSkyStorageIndependent::cSys = cs;
         empty = null;
         setStorage();
      }

      virtual std::string name() const noexcept override
      { return "FullSky"; }

      virtual iterator begin() noexcept override;
      virtual const_iterator begin() const noexcept override;
      virtual iterator end() noexcept override;
      virtual const_iterator end() const noexcept override;

      virtual const_iterator cbegin() const noexcept override;
      virtual const_iterator cend() const noexcept override;

      virtual T Accumulate() const noexcept override;

      virtual const T& GetValue(const size_t hpIndex, const size_t spIndex) const noexcept override;

      virtual void SetValue(const size_t hpIndex, const size_t spIndex, const T& value) noexcept override;

      virtual T& GetReference(const size_t hpIndex, const size_t spIndex) noexcept override;

      virtual void GetSpectrum(const size_t hpIndex, std::vector<T> &sp) const noexcept override;

      virtual Healpix_Map<T> GetHealpixMap(const size_t spIndex) const noexcept override;

      virtual void ApplyFunction(const std::function< void(T&) > &func) noexcept override;
      virtual void ApplyFunction(const std::function< void(const T&) > &func) const noexcept override;
      virtual void ApplyFunction(const std::function< void(T&, size_t) > &func) noexcept override;
      virtual void ApplyFunction(const std::function< void(const T&, size_t) > &func) const noexcept override;

      virtual void ApplyFunction(const std::function< void(T&, const T&) > &func, const BaseSky<T> &map ) override;
      virtual void ApplyFunction(const std::function< void(const T&, const T&) > &func, const BaseSky<T> &map ) const override;
      virtual void ApplyFunction(const std::function< void(T&, const T&, size_t) > &func, const BaseSky<T> &map ) override;
      virtual void ApplyFunction(const std::function< void(const T&, const T&, size_t) > &func, const BaseSky<T> &map ) const override;

      virtual T ApplyFunctionAccumulate(const std::function< T(T&) > &func) noexcept override;
      virtual T ApplyFunctionAccumulate(const std::function< T(const T&) > &func) const noexcept override;
      virtual T ApplyFunctionAccumulate(const std::function< T(T&, size_t) > &func) noexcept override;
      virtual T ApplyFunctionAccumulate(const std::function< T(const T&, size_t) > &func) const noexcept override;

      virtual T ApplyFunctionAccumulate(const std::function< T(T&, const T&) > &func, const BaseSky<T> &map ) override;
      virtual T ApplyFunctionAccumulate(const std::function< T(const T&, const T&) > &func, const BaseSky<T> &map ) const override;
      virtual T ApplyFunctionAccumulate(const std::function< T(T&, const T&, size_t) > &func, const BaseSky<T> &map ) override;
      virtual T ApplyFunctionAccumulate(const std::function< T(const T&, const T&, size_t) > &func, const BaseSky<T> &map ) const override;

      virtual bool IsFullSky() const noexcept override { return true; }

      virtual std::set<size_t> GetActivePixels() const noexcept override
      {
         std::set<size_t> pixels;
         for (size_t i(0); i < size_t(Healpix_Base::Npix()); ++i)
            pixels.insert(i);
         return pixels;
      }

      // This should be deprecated
      template <typename F, typename O>
      friend void ApplyFunction(FullSky<F> &arg1, const FullSky<O> &other, const std::function< void(F&, const O&) > &func);
      template <typename F, typename O>
      friend void ApplyFunction(const FullSky<F> &arg1, const FullSky<O> &other, const std::function< void(const F&, const O&) > &func);

      template <typename F, typename O>
      friend F ApplyFunctionAccumulate(FullSky<F> &arg1, const FullSky<O> &other, const std::function< F(F&, const O&) > &func);
      template <typename F, typename O>
      friend F ApplyFunctionAccumulate(const FullSky<F> &arg1, const FullSky<O> &other, const std::function< F(const F&, const O&) > &func);

      virtual iterator spectraBegin(const size_t hpIndex) noexcept override;
      virtual iterator spectraEnd(const size_t hpIndex) noexcept override;
      virtual const_iterator spectraBegin(const size_t hpIndex) const noexcept override;
      virtual const_iterator spectraEnd(const size_t hpIndex) const noexcept override;

      virtual iterator mapBegin(const size_t spIndex) noexcept override;
      virtual iterator mapEnd(const size_t spIndex) noexcept override;
      virtual const_iterator mapBegin(const size_t spIndex) const noexcept override;
      virtual const_iterator mapEnd(const size_t spIndex) const noexcept override;

};


//Function implementation
template <typename T>
typename FullSky<T>::iterator FullSky<T>::begin() noexcept
{
   return iterator(new FullSky_IteratorImplementation<T>(storage.begin(), storage.begin(), spBin.GetSize()));
}
template <typename T>
typename FullSky<T>::const_iterator FullSky<T>::begin() const noexcept
{
   return const_iterator(new FullSky_IteratorImplementation<const T>(storage.begin(), storage.begin(), spBin.GetSize()));
}
template <typename T>
typename FullSky<T>::iterator FullSky<T>::end() noexcept
{
   return iterator(new FullSky_IteratorImplementation<T>(storage.end(), storage.begin(), spBin.GetSize()));
}
template <typename T>
typename FullSky<T>::const_iterator FullSky<T>::end() const noexcept
{
   return const_iterator(new FullSky_IteratorImplementation<const T>(storage.end(), storage.begin(), spBin.GetSize()));
}
template <typename T>
typename FullSky<T>::const_iterator FullSky<T>::cbegin() const noexcept
{
   return const_iterator(new FullSky_IteratorImplementation<const T>(storage.begin(), storage.begin(), spBin.GetSize()));
}
template <typename T>
typename FullSky<T>::const_iterator FullSky<T>::cend() const noexcept
{
   return const_iterator(new FullSky_IteratorImplementation<const T>(storage.end(), storage.begin(), spBin.GetSize()));
}

template <typename T>
const T& FullSky<T>::GetValue(const size_t hpIndex, const size_t spIndex) const noexcept
{
   return storage[hpIndex*spBin.GetSize()+spIndex];
}

template <typename T>
void FullSky<T>::SetValue(const size_t hpIndex, const size_t spIndex, const T& value) noexcept
{
   storage[hpIndex*spBin.GetSize()+spIndex] = value;
}

template <typename T>
T& FullSky<T>::GetReference(const size_t hpIndex, const size_t spIndex) noexcept
{
   return storage[hpIndex*spBin.GetSize()+spIndex];
}

template <typename T>
void FullSky<T>::GetSpectrum(const size_t hpIndex, std::vector<T> &sp) const noexcept
{
   sp.resize(spBin.GetSize());
   for (size_t i = 0; i < sp.size(); ++i)
      sp[i] = storage[i+hpIndex*spBin.GetSize()];
}

template <typename T>
Healpix_Map<T> FullSky<T>::GetHealpixMap(const size_t spIndex) const noexcept
{
   Healpix_Map<T> out(Order(), Scheme());

   const size_t enSize = spBin.GetSize();
#pragma omp parallel for schedule(static) default(shared)
   for (int i = 0; i < Npix(); ++i)
      out[i] = storage[i*enSize + spIndex];

   return out;
}

template <typename T>
typename FullSky<T>::iterator FullSky<T>::spectraBegin(const size_t hpIndex) noexcept
{
   return iterator(new Spectra_IteratorImplementation<T>(storage.begin()+hpIndex*spBin.GetSize(), storage.begin(), spBin.GetSize()));
}
template <typename T>
typename FullSky<T>::iterator FullSky<T>::spectraEnd(const size_t hpIndex) noexcept
{
   return iterator(new Spectra_IteratorImplementation<T>(storage.begin()+(hpIndex+1)*spBin.GetSize(), storage.begin(), spBin.GetSize()));
}
template <typename T>
typename FullSky<T>::const_iterator FullSky<T>::spectraBegin(const size_t hpIndex) const noexcept
{
   return const_iterator(new Spectra_IteratorImplementation<const T>(storage.begin()+hpIndex*spBin.GetSize(), storage.begin(), spBin.GetSize()));
}
template <typename T>
typename FullSky<T>::const_iterator FullSky<T>::spectraEnd(const size_t hpIndex) const noexcept
{
   return const_iterator(new Spectra_IteratorImplementation<const T>(storage.begin()+(hpIndex+1)*spBin.GetSize(), storage.begin(), spBin.GetSize()));
}

template <typename T>
typename FullSky<T>::iterator FullSky<T>::mapBegin(const size_t spIndex) noexcept
{
   return iterator(new Map_IteratorImplementation<T>(storage.begin()+spIndex, storage.begin(), spBin.GetSize()));
}
template <typename T>
typename FullSky<T>::iterator FullSky<T>::mapEnd(const size_t spIndex) noexcept
{
   return iterator(new Map_IteratorImplementation<T>(storage.end()+spIndex, storage.begin(), spBin.GetSize()));
}
template <typename T>
typename FullSky<T>::const_iterator FullSky<T>::mapBegin(const size_t spIndex) const noexcept
{
   return const_iterator(new Map_IteratorImplementation<const T>(storage.begin()+spIndex, storage.begin(), spBin.GetSize()));
}
template <typename T>
typename FullSky<T>::const_iterator FullSky<T>::mapEnd(const size_t spIndex) const noexcept
{
   return const_iterator(new Map_IteratorImplementation<const T>(storage.end()+spIndex, storage.begin(), spBin.GetSize()));
}

template <typename T>
T FullSky<T>::Accumulate() const noexcept
{
   //Use Kahan summing, even though it is slower for integers
   T sum(0);
   T c(0);
#if (_OPENMP >= 201307)
#pragma omp parallel for simd schedule(static) default(shared) reduction(+:sum) firstprivate(c)
#else
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
#endif
   for (size_t i = 0; i < storage.size(); ++i) {
      const T y = storage[i] - c;
      const T t = sum + y;
      c = t - sum;
      c -= y;
      sum = t;
   }

   return sum;
}

template <typename T>
void FullSky<T>::ApplyFunction(const std::function< void(T&) > &func) noexcept
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;

   //Try to get a function pointer for simd evaluation
   auto ptr = func.template target<void(*)(T&)>();

   if (ptr != nullptr) {
#if (_OPENMP >= 201307)
#pragma omp parallel for simd schedule(static) default(shared)
#else
#pragma omp parallel for schedule(static) default(shared)
#endif
      for (size_t i = 0; i < storage.size(); ++i)
         (*ptr)(storage[i]);
   }
   else
   {
#pragma omp parallel for schedule(static) default(shared)
      for (size_t i = 0; i < storage.size(); ++i)
         func(storage[i]);
   }
}

template <typename T>
void FullSky<T>::ApplyFunction(const std::function< void(const T&) > &func) const noexcept
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;

   //Try to get a function pointer for simd evaluation
   auto ptr = func.template target<void(*)(const T&)>();

   if (ptr != nullptr) {
#if (_OPENMP >= 201307)
#pragma omp parallel for simd schedule(static) default(shared)
#else
#pragma omp parallel for schedule(static) default(shared)
#endif
      for (size_t i = 0; i < storage.size(); ++i)
         (*ptr)(storage[i]);
   }
   else
   {
#pragma omp parallel for schedule(static) default(shared)
      for (size_t i = 0; i < storage.size(); ++i)
         func(storage[i]);
   }
}

template <typename T>
void FullSky<T>::ApplyFunction(const std::function< void(T&, size_t) > &func) noexcept
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;

   //Try to get a function pointer for simd evaluation
   auto ptr = func.template target<void(*)(T&,size_t)>();

   const auto enSize = spBin.GetSize();
   const auto mapSize = Npix();
   if (ptr != nullptr) {
#pragma omp parallel for schedule(static) default(shared)
      for (int i = 0; i < mapSize; ++i)
#if (_OPENMP >= 201307)
#pragma omp simd
#endif
         for (size_t j = 0; j < enSize; ++j)
            (*ptr)(storage[i*enSize + j], j);
   }

   else

   {
#pragma omp parallel for schedule(static) default(shared)
      for (int i = 0; i < mapSize; ++i)
         for (size_t j = 0; j < enSize; ++j)
            func(storage[i*enSize + j], j);
   }
}

template <typename T>
void FullSky<T>::ApplyFunction(const std::function< void(const T&, size_t) > &func) const noexcept
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;

   //Try to get a function pointer for simd evaluation
   auto ptr = func.template target<void(*)(const T&,size_t)>();

   const auto enSize = spBin.GetSize();
   const auto mapSize = Npix();
   if (ptr != nullptr) {
#pragma omp parallel for schedule(static) default(shared)
      for (int i = 0; i < mapSize; ++i)
#if (_OPENMP >= 201307)
#pragma omp simd
#endif
         for (size_t j = 0; j < enSize; ++j)
            (*ptr)(storage[i*enSize + j], j);
   }

   else

   {
#pragma omp parallel for schedule(static) default(shared)
      for (int i = 0; i < mapSize; ++i)
         for (size_t j = 0; j < enSize; ++j)
            func(storage[i*enSize + j], j);
   }
}

template <typename T>
T FullSky<T>::ApplyFunctionAccumulate(const std::function< T(T&) > &func) noexcept
{
   //Use Kahan summing
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;

   //Try to get a function pointer for simd evaluation
   auto ptr = func.template target<T(*)(T&)>();

   if (ptr != nullptr) {
#if (_OPENMP >= 201307)
#pragma omp parallel for simd schedule(static) default(shared) reduction(+:sum) firstprivate(c)
#else
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
#endif
      for (size_t i = 0; i < storage.size(); ++i) {
         const T y = (*ptr)(storage[i]) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }
   }
   else
   {
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
      for (size_t i = 0; i < storage.size(); ++i) {
         const T y = func(storage[i]) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }
   }
   return sum;
}

template <typename T>
T FullSky<T>::ApplyFunctionAccumulate(const std::function< T(const T&) > &func) const noexcept
{
   //Use Kahan summing
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;

   //Try to get a function pointer for simd evaluation
   auto ptr = func.template target<T(*)(const T&)>();

   if (ptr != nullptr) {
#if (_OPENMP >= 201307)
#pragma omp parallel for simd schedule(static) default(shared) reduction(+:sum) firstprivate(c)
#else
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
#endif
      for (size_t i = 0; i < storage.size(); ++i) {
         const T y = (*ptr)(storage[i]) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }
   }
   else
   {
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
      for (size_t i = 0; i < storage.size(); ++i) {
         const T y = func(storage[i]) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }
   }
   return sum;
}

template <typename T>
T FullSky<T>::ApplyFunctionAccumulate(const std::function< T(T&, size_t) > &func) noexcept
{
   //Use Kahan summing
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;

   //Try to get a function pointer for simd evaluation
   auto ptr = func.template target<T(*)(T&,size_t)>();

   const auto enSize = spBin.GetSize();
   const auto mapSize = Npix();

   if (ptr != nullptr) {
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
      for (int i = 0; i < mapSize; ++i) {
#if (_OPENMP >= 201307)
#pragma omp simd reduction(+:sum)
#endif
         for (size_t j = 0; j < enSize; ++j) {
            const T y = (*ptr)(storage[i*enSize+j], j) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }
      }
   }
   else
   {
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
      for (int i = 0; i < mapSize; ++i) {
         for (size_t j = 0; j < enSize; ++j) {
            const T y = func(storage[i*enSize+j], j) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }
      }
   }
   return sum;
}

template <typename T>
T FullSky<T>::ApplyFunctionAccumulate(const std::function< T(const T&, size_t) > &func) const noexcept
{
   //Use Kahan summing
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;

   //Try to get a function pointer for simd evaluation
   auto ptr = func.template target<T(*)(const T&,size_t)>();

   const auto enSize = spBin.GetSize();
   const auto mapSize = Npix();

   if (ptr != nullptr) {
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
      for (int i = 0; i < mapSize; ++i) {
#if (_OPENMP >= 201307)
#pragma omp simd reduction(+:sum)
#endif
         for (size_t j = 0; j < enSize; ++j) {
            const T y = (*ptr)(storage[i*enSize+j], j) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }
      }
   }
   else
   {
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
      for (int i = 0; i < mapSize; ++i) {
         for (size_t j = 0; j < enSize; ++j) {
            const T y = func(storage[i*enSize+j], j) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }
      }
   }
   return sum;
}


template <typename T>
void FullSky<T>::ApplyFunction(const std::function< void(T&, const T&) > &func, const BaseSky<T> &map)
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;

   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

   //Try to cast the reference to a FullSky instance, otherwise use the general interface
   try {
      const auto &fsMap = dynamic_cast< const FullSky<T>& > ( map );

      //FullSky interface with possible simd

      //Try to get a function pointer for simd evaluation
      auto ptr = func.template target<void(*)(T&, const T&)>();

      if (ptr != nullptr) {
#if (_OPENMP >= 201307)
#pragma omp parallel for simd schedule(static) default(shared)
#else
#pragma omp parallel for schedule(static) default(shared)
#endif
         for (size_t i = 0; i < storage.size(); ++i)
            (*ptr)(storage[i], fsMap.storage[i]);
      }
      else
      {
#pragma omp parallel for schedule(static) default(shared)
         for (size_t i = 0; i < storage.size(); ++i)
            func(storage[i], fsMap.storage[i]);
      }

   } catch (const std::bad_cast &e) {

      //Generic interface
      const auto enSize = spBin.GetSize();
      const auto mapSize = Npix();

#pragma omp parallel for schedule(static) default(shared)
      for (int i = 0; i < mapSize; ++i) 
         for (size_t j = 0; j < enSize; ++j) 
            func(storage[i*enSize+j], map.GetValue(i,j));

   } 
}

template <typename T>
void FullSky<T>::ApplyFunction(const std::function< void(const T&, const T&) > &func, const BaseSky<T> &map) const
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;

   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

   //Try to cast the reference to a FullSky instance, otherwise use the general interface
   try {
      const auto &fsMap = dynamic_cast< const FullSky<T>& > ( map );

      //FullSky interface with possible simd

      //Try to get a function pointer for simd evaluation
      auto ptr = func.template target<void(*)(const T&, const T&)>();

      if (ptr != nullptr) {
#if (_OPENMP >= 201307)
#pragma omp parallel for simd schedule(static) default(shared)
#else
#pragma omp parallel for schedule(static) default(shared)
#endif
         for (size_t i = 0; i < storage.size(); ++i)
            (*ptr)(storage[i], fsMap.storage[i]);
      }
      else
      {
#pragma omp parallel for schedule(static) default(shared)
         for (size_t i = 0; i < storage.size(); ++i)
            func(storage[i], fsMap.storage[i]);
      }

   } catch (const std::bad_cast &e) {
      
      //Generic interface
      const auto enSize = spBin.GetSize();
      const auto mapSize = Npix();

#pragma omp parallel for schedule(static) default(shared)
      for (int i = 0; i < mapSize; ++i) 
         for (size_t j = 0; j < enSize; ++j) 
            func(storage[i*enSize+j], map.GetValue(i,j));

   } 
}

template <typename T>
void FullSky<T>::ApplyFunction(const std::function< void(T&, const T&, size_t) > &func, const BaseSky<T> &map)
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;

   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

   const auto enSize = spBin.GetSize();
   const auto mapSize = Npix();

   //Try to cast the reference to a FullSky instance, otherwise use the general interface
   try {
      const auto &fsMap = dynamic_cast< const FullSky<T>& > ( map );

      //FullSky interface with possible simd

      //Try to get a function pointer for simd evaluation
      auto ptr = func.template target<void(*)(T&, const T&, size_t)>();

      if (ptr != nullptr) {
#pragma omp parallel for schedule(static) default(shared)
         for (int i = 0; i < mapSize; ++i) 
#if (_OPENMP >= 201307)
#pragma omp simd
#endif
            for (size_t j = 0; j < enSize; ++j) 
               (*ptr)(storage[i*enSize+j], fsMap.storage[i*enSize+j], j);
      }
      else
      {
#pragma omp parallel for schedule(static) default(shared)
         for (int i = 0; i < mapSize; ++i) 
            for (size_t j = 0; j < enSize; ++j) 
               func(storage[i*enSize+j], fsMap.storage[i*enSize+j], j);
      }

   } catch (const std::bad_cast &e) {

      //Generic interface
#pragma omp parallel for schedule(static) default(shared)
      for (int i = 0; i < mapSize; ++i) 
         for (size_t j = 0; j < enSize; ++j) 
            func(storage[i*enSize+j], map.GetValue(i,j), j);

   } 
}

template <typename T>
void FullSky<T>::ApplyFunction(const std::function< void(const T&, const T&, size_t) > &func, const BaseSky<T> &map) const
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;

   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

   const auto enSize = spBin.GetSize();
   const auto mapSize = Npix();

   //Try to cast the reference to a FullSky instance, otherwise use the general interface
   try {
      const auto &fsMap = dynamic_cast< const FullSky<T>& > ( map );

      //FullSky interface with possible simd

      //Try to get a function pointer for simd evaluation
      auto ptr = func.template target<void(*)(const T&, const T&, size_t)>();

      if (ptr != nullptr) {
#pragma omp parallel for schedule(static) default(shared)
         for (int i = 0; i < mapSize; ++i) 
#if (_OPENMP >= 201307)
#pragma omp simd
#endif
            for (size_t j = 0; j < enSize; ++j) 
               (*ptr)(storage[i*enSize+j], fsMap.storage[i*enSize+j], j);
      }
      else
      {
#pragma omp parallel for schedule(static) default(shared)
         for (int i = 0; i < mapSize; ++i) 
            for (size_t j = 0; j < enSize; ++j) 
               func(storage[i*enSize+j], fsMap.storage[i*enSize+j], j);
      }

   } catch (const std::bad_cast &e) {

      //Generic interface
#pragma omp parallel for schedule(static) default(shared)
      for (int i = 0; i < mapSize; ++i) 
         for (size_t j = 0; j < enSize; ++j) 
            func(storage[i*enSize+j], map.GetValue(i,j), j);

   } 
}

template <typename T>
T FullSky<T>::ApplyFunctionAccumulate(const std::function< T(T&, const T&) > &func, const BaseSky<T> &map)
{
   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

   //Use Kahan summing
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;

   //Try to cast the reference to a FullSky instance, otherwise use the general interface
   try {
      const auto &fsMap = dynamic_cast< const FullSky<T>& > ( map );

      //FullSky interface with possible simd

      //Try to get a function pointer for simd evaluation
      auto ptr = func.template target<T(*)(T&, const T&)>();

      if (ptr != nullptr) {
#if (_OPENMP >= 201307)
#pragma omp parallel for simd schedule(static) default(shared) reduction(+:sum) firstprivate(c)
#else
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
#endif
         for (size_t i = 0; i < storage.size(); ++i) {
            const T y = (*ptr)(storage[i], fsMap.storage[i]) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }
      }
      else
      {
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
         for (size_t i = 0; i < storage.size(); ++i) {
            const T y = func(storage[i], fsMap.storage[i]) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }
      }

   } catch (const std::bad_cast &e) {

      //Generic interface
      const auto enSize = spBin.GetSize();
      const auto mapSize = Npix();

#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
      for (int i = 0; i < mapSize; ++i) {
         for (size_t j = 0; j < enSize; ++j) {
            const T y = func(storage[i*enSize+j], map.GetValue(i,j)) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }
      }

   } 

   return sum;
}

template <typename T>
T FullSky<T>::ApplyFunctionAccumulate(const std::function< T(const T&, const T&) > &func, const BaseSky<T> &map) const
{
   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

   //Use Kahan summing
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;

   //Try to cast the reference to a FullSky instance, otherwise use the general interface
   try {
      const auto &fsMap = dynamic_cast< const FullSky<T>& > ( map );

      //FullSky interface with possible simd

      //Try to get a function pointer for simd evaluation
      auto ptr = func.template target<T(*)(const T&, const T&)>();

      if (ptr != nullptr) {
#if (_OPENMP >= 201307)
#pragma omp parallel for simd schedule(static) default(shared) reduction(+:sum) firstprivate(c)
#else
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
#endif
         for (size_t i = 0; i < storage.size(); ++i) {
            const T y = (*ptr)(storage[i], fsMap.storage[i]) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }
      }
      else
      {
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
         for (size_t i = 0; i < storage.size(); ++i) {
            const T y = func(storage[i], fsMap.storage[i]) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }
      }

   } catch (const std::bad_cast &e) {

      //Generic interface
      const auto enSize = spBin.GetSize();
      const auto mapSize = Npix();

#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
      for (int i = 0; i < mapSize; ++i) {
         for (size_t j = 0; j < enSize; ++j) {
            const T y = func(storage[i*enSize+j], map.GetValue(i,j)) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }
      }

   } 

   return sum;
}

template <typename T>
T FullSky<T>::ApplyFunctionAccumulate(const std::function< T(T&, const T&, size_t) > &func, const BaseSky<T> &map)
{
   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

   //Use Kahan summing
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;

   const auto enSize = spBin.GetSize();
   const auto mapSize = Npix();

   //Try to cast the reference to a FullSky instance, otherwise use the general interface
   try {
      const auto &fsMap = dynamic_cast< const FullSky<T>& > ( map );

      //FullSky interface with possible simd

      //Try to get a function pointer for simd evaluation
      auto ptr = func.template target<T(*)(T&, const T&, size_t)>();

      if (ptr != nullptr) {
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
         for (int i = 0; i < mapSize; ++i) 
#if (_OPENMP >= 201307)
#pragma omp simd
#endif
            for (size_t j = 0; j < enSize; ++j) {
               const T y = (*ptr)(storage[i*enSize + j], fsMap.storage[i*enSize+j], j) - c;
               const T t = sum + y;
               c = t - sum;
               c -= y;
               sum = t;
            }
      }
      else
      {
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
         for (int i = 0; i < mapSize; ++i) 
            for (size_t j = 0; j < enSize; ++j) {
               const T y = func(storage[i*enSize + j], fsMap.storage[i*enSize+j], j) - c;
               const T t = sum + y;
               c = t - sum;
               c -= y;
               sum = t;
            }
      }

   } catch (const std::bad_cast &e) {

      //Generic interface
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
      for (int i = 0; i < mapSize; ++i) {
         for (size_t j = 0; j < enSize; ++j) {
            const T y = func(storage[i*enSize+j], map.GetValue(i,j), j) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }
      }

   } 

   return sum;
}

template <typename T>
T FullSky<T>::ApplyFunctionAccumulate(const std::function< T(const T&, const T&, size_t) > &func, const BaseSky<T> &map) const
{
   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

   //Use Kahan summing
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;

   const auto enSize = spBin.GetSize();
   const auto mapSize = Npix();

   //Try to cast the reference to a FullSky instance, otherwise use the general interface
   try {
      const auto &fsMap = dynamic_cast< const FullSky<T>& > ( map );

      //FullSky interface with possible simd

      //Try to get a function pointer for simd evaluation
      auto ptr = func.template target<T(*)(const T&, const T&, size_t)>();

      if (ptr != nullptr) {
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
         for (int i = 0; i < mapSize; ++i) 
#if (_OPENMP >= 201307)
#pragma omp simd
#endif
            for (size_t j = 0; j < enSize; ++j) {
               const T y = (*ptr)(storage[i*enSize + j], fsMap.storage[i*enSize+j], j) - c;
               const T t = sum + y;
               c = t - sum;
               c -= y;
               sum = t;
            }
      }
      else
      {
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
         for (int i = 0; i < mapSize; ++i) 
            for (size_t j = 0; j < enSize; ++j) {
               const T y = func(storage[i*enSize + j], fsMap.storage[i*enSize+j], j) - c;
               const T t = sum + y;
               c = t - sum;
               c -= y;
               sum = t;
            }
      }

   } catch (const std::bad_cast &e) {

      //Generic interface
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
      for (int i = 0; i < mapSize; ++i) {
         for (size_t j = 0; j < enSize; ++j) {
            const T y = func(storage[i*enSize+j], map.GetValue(i,j), j) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }
      }

   } 

   return sum;
}




template <typename T, typename O>
void ApplyFunction(FullSky<T> &arg1, const FullSky<O> &arg2, const std::function< void(T&, const O&) > &func)
{
   if (! arg1.Equivalent(arg2)) {
      ERROR("Maps should be equivalent");
      throw(std::invalid_argument("Equivalence test failed"));
   }

   if ( func == nullptr ) {
      ERROR("Invalid function");
      throw(std::invalid_argument("Empty function"));
   }

   //Try to get a function pointer for simd evaluation
   auto ptr = func.template target<void(*)(T&, const O&)>();

   if (ptr != nullptr) {
#if (_OPENMP >= 201307)
#pragma omp parallel for simd schedule(static) default(shared)
#else
#pragma omp parallel for schedule(static) default(shared)
#endif
      for (size_t i = 0; i < arg1.storage.size(); ++i)
         (*ptr)(arg1.storage[i], arg2.storage[i]);
   }
   else
   {
#pragma omp parallel for schedule(static) default(shared)
      for (size_t i = 0; i < arg1.storage.size(); ++i)
         func(arg1.storage[i], arg2.storage[i]);
   }
}

/** \brief Apply function to all values in maps
 *
 * The maps have to be equivalent, otherwise we throw an exception.
 * \param func should be threadsafe
 *
 * Specialized implementation for two FullSky maps
 */
template <typename T, typename O>
void ApplyFunction(const FullSky<T> &arg1, const FullSky<O> &arg2, const std::function< void(const T&, const O&) > &func)
{
   if (! arg1.Equivalent(arg2)) {
      ERROR("Maps should be equivalent");
      throw(std::invalid_argument("Equivalence test failed"));
   }

   if ( func == nullptr ) {
      ERROR("Invalid function");
      throw(std::invalid_argument("Empty function"));
   }

   //Try to get a function pointer for simd evaluation
   auto ptr = func.template target<void(*)(const T&, const O&)>();

   if (ptr != nullptr) {
#if (_OPENMP >= 201307)
#pragma omp parallel for simd schedule(static) default(shared)
#else
#pragma omp parallel for schedule(static) default(shared)
#endif
      for (size_t i = 0; i < arg1.storage.size(); ++i)
         (*ptr)(arg1.storage[i], arg2.storage[i]);
   }
   else
   {
#pragma omp parallel for schedule(static) default(shared)
      for (size_t i = 0; i < arg1.storage.size(); ++i)
         func(arg1.storage[i], arg2.storage[i]);
   }
}

/** \brief Apply function to all values in maps
 *
 * For this to work the type T has to be usable in openmp reduce clause
 * \param func should be threadsafe
 *
 * Specialized implementation for two FullSky maps
 */
template <typename T, typename O>
T ApplyFunctionAccumulate(FullSky<T> &arg1, const FullSky<O> &arg2, const std::function< T(T&, const O&) > &func)
{
   if (! arg1.Equivalent(arg2)) {
      ERROR("Maps should be equivalent");
      throw(std::invalid_argument("Equivalence test failed"));
   }

   if ( func == nullptr ) {
      ERROR("Invalid function");
      throw(std::invalid_argument("Empty function"));
   }

   //Use Kahan summing
   T sum(0);
   T c(0);

   //Try to get a function pointer for simd evaluation
   auto ptr = func.template target<T(*)(T&, const O&)>();

   if (ptr != nullptr) {
#if (_OPENMP >= 201307)
#pragma omp parallel for simd schedule(static) default(shared) reduction(+:sum) firstprivate(c)
#else
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
#endif
      for (size_t i = 0; i < arg1.storage.size(); ++i) {
         const T y = (*ptr)(arg1.storage[i], arg2.storage[i]) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }
   }
   else {
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
      for (size_t i = 0; i < arg1.storage.size(); ++i) {
         const T y = func(arg1.storage[i], arg2.storage[i]) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }
   }

   return sum;
}

template <typename T, typename O>
T ApplyFunctionAccumulate(const FullSky<T> &arg1, const FullSky<O> &arg2, const std::function< T(const T&, const O&) > &func)
{
   if (! arg1.Equivalent(arg2)) {
      ERROR("Maps should be equivalent");
      throw(std::invalid_argument("Equivalence test failed"));
   }

   if ( func == nullptr ) {
      ERROR("Invalid function");
      throw(std::invalid_argument("Empty function"));
   }

   //Use Kahan summing
   T sum(0);
   T c(0);

   //Try to get a function pointer for simd evaluation
   auto ptr = func.template target<T(*)(const T&, const O&)>();

   if (ptr != nullptr) {
#if (_OPENMP >= 201307)
#pragma omp parallel for simd schedule(static) default(shared) reduction(+:sum) firstprivate(c)
#else
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
#endif
      for (size_t i = 0; i < arg1.storage.size(); ++i) {
         const T y = (*ptr)(arg1.storage[i], arg2.storage[i]) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }
   }
   else {
#pragma omp parallel for schedule(static) default(shared) reduction(+:sum) firstprivate(c)
      for (size_t i = 0; i < arg1.storage.size(); ++i) {
         const T y = func(arg1.storage[i], arg2.storage[i]) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }
   }

   return sum;
}

}

#endif
