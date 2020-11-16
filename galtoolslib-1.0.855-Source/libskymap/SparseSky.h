#ifndef SPARSE_SKY_H
#define SPARSE_SKY_H

#include "BaseSky.h"

#include <map>
#include <tuple>
#include <type_traits>

namespace SM {


/** Skymap that stores individual values
 *
 * Implements the storage as a vector of maps
 * Stores each healpix map in a vector, very suitable for point sources in Fermi.
 *
 */
template <typename T>
class SparseSky : public BaseSky<T>
{
   protected:
      std::vector< std::map<size_t,T> > storage;

      /** \brief Iterator implementation that covers the entire sky. 
       */
      template < typename C >
      class SparseSky_IteratorImplementation : public SM::BaseSky<T>::template IteratorImplementation<C> {
         protected:
            //Declare the iterator const_iterator if C is const
            //Use this for keeping track of our position within a sky
            using iterator_type = typename std::conditional<
                  std::is_const<C>::value, 
                  typename std::map<size_t,T>::const_iterator, 
                  typename std::map<size_t,T>::iterator>::type;
            iterator_type it;

            using data_type = typename std::conditional<
                  std::is_const<C>::value, 
                  const std::vector<std::map<size_t,T> >,
                  std::vector<std::map<size_t,T> > >::type;
            data_type &data_ref;

            //Store this for convenience in the decrement operator
            iterator_type firstIt;

            //Keep track of where we are in energy
            size_t eIndex;

         public:
            //Access first or last element
            SparseSky_IteratorImplementation( data_type &data, bool first ) 
               : data_ref(data)
            {
               //We always need to find firstIt
               eIndex = 0;
               it = data_ref[eIndex].begin();
               while ( eIndex+1 < data_ref.size() && it == data_ref[eIndex].end() ) {
                  ++eIndex;
                  it = data_ref[eIndex].begin();
               }
               firstIt = it;

               if (! first) {
                  //Simply point to the end of the last map
                  eIndex = data_ref.size()-1;
                  it = data_ref[eIndex].end();
               }
            }


            virtual void increment() override 
            { 
               ++it;
               while ( eIndex+1 < data_ref.size() && it == data_ref[eIndex].end() ) {
                  ++eIndex;
                  it = data_ref[eIndex].begin();
               }
            }
            virtual void decrement() override 
            { 
               //Make sure we do not pass firstIt
               if (it != firstIt) {
                  while (eIndex > 0 && it == data_ref[eIndex].begin()) {
                     --eIndex;
                     it = data_ref[eIndex].end();
                  }
                  --it; 
               }
            }

            virtual bool equal(const typename SM::BaseSky<T>::template IteratorImplementation<C>& other) const override { 
               return it == dynamic_cast<const SparseSky_IteratorImplementation<C>& >(other).it; 
            }
            
            virtual C& dereference() const override {
               return it->second;
            }

            virtual size_t healpixIndex() const override {
               return it->first;
            }

            virtual size_t energyIndex() const override {
               return eIndex;
            }

            virtual SparseSky_IteratorImplementation<C>* clone() const override {
               //Need new copies of the iterator
               auto copy = new SparseSky_IteratorImplementation<C>( data_ref, false );
               copy->it = iterator_type(it);
               copy->eIndex = eIndex;
               return copy;
            }

      }; 

      /** \brief Iterator for stepping over maps in a single energy bin
       */
      template< typename C >
      class Map_IteratorImplementation : public SparseSky_IteratorImplementation<C> {
         protected:
            using iterator_type = typename SparseSky_IteratorImplementation<C>::iterator_type;
            using data_type = typename SparseSky_IteratorImplementation<C>::data_type;
            using SparseSky_IteratorImplementation<C>::it;
            using SparseSky_IteratorImplementation<C>::data_ref;
            using SparseSky_IteratorImplementation<C>::eIndex;
            
         public:
            Map_IteratorImplementation( data_type &data, bool first, size_t ie ) 
               : SparseSky_IteratorImplementation<C>( data, first )
            {
               eIndex = ie;
               if (first) 
                  it = data_ref[eIndex].begin();
               else
                  it = data_ref[eIndex].end();
            }

            virtual void increment() override { 
               ++it;
            }
            virtual void decrement() override { 
               --it;
            }
      };

      /** \brief Iterator for stepping over spectra in a single pixel
       */
      template< typename C >
      class Spectra_IteratorImplementation : public SparseSky_IteratorImplementation<C> {
         protected:
            using iterator_type = typename SparseSky_IteratorImplementation<C>::iterator_type;
            using data_type = typename SparseSky_IteratorImplementation<C>::data_type;
            using SparseSky_IteratorImplementation<C>::it;
            using SparseSky_IteratorImplementation<C>::data_ref;
            using SparseSky_IteratorImplementation<C>::eIndex;
            using SparseSky_IteratorImplementation<C>::firstIt;

            const size_t hpIndex;

         public:
            Spectra_IteratorImplementation( data_type &data, bool first, size_t pix ) 
               : SparseSky_IteratorImplementation<C>( data, first )
               , hpIndex(pix)
            {
               //Need to set the firstIt
               size_t ei = 0;
               firstIt = data_ref[ei].find(hpIndex);
               while ( ei+1 < data_ref.size() && firstIt == data_ref[ei].end() ) {
                  ++ei;
                  firstIt = data_ref[ei].find(hpIndex);
               }

               if (first) {
                  it = firstIt;
                  eIndex = ei;
               }
               //The end is set properly in the original constructor
            }

            virtual void increment() override { 
               if ( eIndex+1 == data_ref.size() ) {
                  it = data_ref[eIndex].end();
               } else {
                  ++eIndex;
                  it = data_ref[eIndex].find(hpIndex);
                  while ( eIndex+1 < data_ref.size() && it == data_ref[eIndex].end() ) {
                     ++eIndex;
                     it = data_ref[eIndex].find(hpIndex);
                  }
               }
            }
            virtual void decrement() override { 
               //Need to be careful never to go beyond the first iterator
               if (it != firstIt) {
                  --eIndex;
                  it = data_ref[eIndex].find(hpIndex);
                  while ( eIndex > 0 && it == data_ref[eIndex].end() ) {
                     --eIndex;
                     it = data_ref[eIndex].find(hpIndex);
                  }
               }
            }
      };

      using BaseSky<T>::spBin;
      using BaseSky<T>::empty;

      void setStorage()
      {
         storage.resize(spBin.GetSize());
         for (auto &m : storage)
            m.clear();
      }

      static std::string getName() { return "SparseSky"; }

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
      using BaseSky<T>::GetValue;
      using BaseSky<T>::SetValue;
      using BaseSky<T>::GetReference;
      using BaseSky<T>::GetSpectrum;
      using BaseSky<T>::operator[];
      using BaseSky<T>::GetSkyInterpolatedSpectra;
      using BaseSky<T>::Equivalent;

      SparseSky(SpectralBinning&& sp, int order, Healpix_Ordering_Scheme scheme, CoordSys cs, const T& null = T(0)) :
         BaseSky<T>(std::forward<SpectralBinning>(sp), order, scheme, cs, null)
      { setStorage(); }
      SparseSky(const SpectralBinning& sp, int order, Healpix_Ordering_Scheme scheme, CoordSys cs, const T& null = T(0)) :
         BaseSky<T>(sp, order, scheme, cs, null)
      { setStorage(); }
      SparseSky(const BaseSkyStorageIndependent &templ, const T& null) :
         BaseSky<T>(templ, null)
      { setStorage(); }

      virtual std::unique_ptr< BaseSky<T> > clone() const override {
         return std::unique_ptr< BaseSky<T> >(new SparseSky<T>(*this));
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
      { return getName(); }

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

      virtual void GetSpectrum (const size_t hpIndex, std::vector<T> &sp) const noexcept override;

      virtual void ApplyFunction(const std::function< void(T&) > &func) noexcept override;
      virtual void ApplyFunction(const std::function< void(const T&) > &func) const noexcept override;
      virtual void ApplyFunction(const std::function< void(T&, size_t) > &func) noexcept override;
      virtual void ApplyFunction(const std::function< void(const T&, size_t) > &func) const noexcept override;

      virtual void ApplyFunction( const std::function< void(T&, const T&) > &func, const BaseSky<T> &map ) override;
      virtual void ApplyFunction( const std::function< void(const T&, const T&) > &func, const BaseSky<T> &map ) const override;
      virtual void ApplyFunction( const std::function< void(T&, const T&, size_t) > &func, const BaseSky<T> &map ) override;
      virtual void ApplyFunction( const std::function< void(const T&, const T&, size_t) > &func, const BaseSky<T> &map ) const override;

      virtual T ApplyFunctionAccumulate(const std::function< T(T&) > &func) noexcept override;
      virtual T ApplyFunctionAccumulate(const std::function< T(const T&) > &func) const noexcept override;
      virtual T ApplyFunctionAccumulate(const std::function< T(T&, size_t) > &func) noexcept override;
      virtual T ApplyFunctionAccumulate(const std::function< T(const T&, size_t) > &func) const noexcept override;

      virtual T ApplyFunctionAccumulate( const std::function< T(T&, const T&) > &func, const BaseSky<T> &map ) override;
      virtual T ApplyFunctionAccumulate( const std::function< T(const T&, const T&) > &func, const BaseSky<T> &map ) const override;
      virtual T ApplyFunctionAccumulate( const std::function< T(T&, const T&, size_t) > &func, const BaseSky<T> &map ) override;
      virtual T ApplyFunctionAccumulate( const std::function< T(const T&, const T&, size_t) > &func, const BaseSky<T> &map ) const override;

      virtual bool IsFullSky() const noexcept override { return false; }

      virtual std::set<size_t> GetActivePixels() const noexcept override;

      // Specialization so we only set as many pixels as needed and do OpenMP
      template <typename F, typename O>
      friend void ApplyFunction(SparseSky<F> &arg1, const SparseSky<O> &other, const std::function< void(F&, const O&) > &func);
      template <typename F, typename O>
      friend void ApplyFunction(const SparseSky<F> &arg1, const SparseSky<O> &other, const std::function< void(const F&, const O&) > &func);

      template <typename F, typename O>
      friend F ApplyFunctionAccumulate(SparseSky<F> &arg1, const SparseSky<O> &other, const std::function< F(F&, const O&) > &func);
      template <typename F, typename O>
      friend F ApplyFunctionAccumulate(const SparseSky<F> &arg1, const SparseSky<O> &other, const std::function< F(const F&, const O&) > &func);

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
typename SparseSky<T>::iterator SparseSky<T>::begin() noexcept
{
   return iterator(new SparseSky_IteratorImplementation<T>(storage, true));
}
template <typename T>
typename SparseSky<T>::const_iterator SparseSky<T>::begin() const noexcept
{
   return const_iterator(new SparseSky_IteratorImplementation<const T>(storage, true));
}
template <typename T>
typename SparseSky<T>::iterator SparseSky<T>::end() noexcept
{
   return iterator(new SparseSky_IteratorImplementation<T>(storage, false));
}
template <typename T>
typename SparseSky<T>::const_iterator SparseSky<T>::end() const noexcept
{
   return const_iterator(new SparseSky_IteratorImplementation<const T>(storage, false));
}
template <typename T>
typename SparseSky<T>::const_iterator SparseSky<T>::cbegin() const noexcept
{
   return const_iterator(new SparseSky_IteratorImplementation<const T>(storage, true));
}
template <typename T>
typename SparseSky<T>::const_iterator SparseSky<T>::cend() const noexcept
{
   return const_iterator(new SparseSky_IteratorImplementation<const T>(storage, false));
}

template <typename T>
const T& SparseSky<T>::GetValue(const size_t hpIndex, const size_t spIndex) const noexcept
{
   auto it = storage[spIndex].find(hpIndex);
   if (it != storage[spIndex].end())
      return it->second;
   else
      return empty;
}

template <typename T>
void SparseSky<T>::SetValue(const size_t hpIndex, const size_t spIndex, const T& value) noexcept
{
   storage[spIndex][hpIndex] = value;
}

template <typename T>
T& SparseSky<T>::GetReference(const size_t hpIndex, const size_t spIndex) noexcept
{
   auto it = storage[spIndex].find(hpIndex);
   if (it == storage[spIndex].end()) {
      T& val = storage[spIndex][hpIndex];
      val = empty;
      return val;
   } else {
      return it->second;
   }
}

template <typename T>
void SparseSky<T>::GetSpectrum(const size_t hpIndex, std::vector<T> &sp) const noexcept
{
   sp.resize(spBin.GetSize());

   for (size_t spIndex = 0; spIndex < spBin.GetSize(); ++spIndex) {
      
      auto it = storage[spIndex].find(hpIndex);
      if (it == storage[spIndex].end())
         sp[spIndex] = empty;
      else
         sp[spIndex] = it->second;

   }

}

template <typename T>
typename SparseSky<T>::iterator SparseSky<T>::spectraBegin(const size_t hpIndex) noexcept
{
   return iterator(new Spectra_IteratorImplementation<T>(storage, true, hpIndex));
}
template <typename T>
typename SparseSky<T>::const_iterator SparseSky<T>::spectraBegin(const size_t hpIndex) const noexcept
{
   return const_iterator(new Spectra_IteratorImplementation<const T>(storage, true, hpIndex));
}
template <typename T>
typename SparseSky<T>::iterator SparseSky<T>::spectraEnd(const size_t hpIndex) noexcept
{
   return iterator(new Spectra_IteratorImplementation<T>(storage, false, hpIndex));
}
template <typename T>
typename SparseSky<T>::const_iterator SparseSky<T>::spectraEnd(const size_t hpIndex) const noexcept
{
   return const_iterator(new Spectra_IteratorImplementation<const T>(storage, false, hpIndex));
}

template <typename T>
typename SparseSky<T>::iterator SparseSky<T>::mapBegin(const size_t spIndex) noexcept
{
   return iterator(new Map_IteratorImplementation<T>(storage, true, spIndex));
}
template <typename T>
typename SparseSky<T>::const_iterator SparseSky<T>::mapBegin(const size_t spIndex) const noexcept
{
   return const_iterator(new Map_IteratorImplementation<const T>(storage, true, spIndex));
}
template <typename T>
typename SparseSky<T>::iterator SparseSky<T>::mapEnd(const size_t spIndex) noexcept
{
   return iterator(new Map_IteratorImplementation<T>(storage, false, spIndex));
}
template <typename T>
typename SparseSky<T>::const_iterator SparseSky<T>::mapEnd(const size_t spIndex) const noexcept
{
   return const_iterator(new Map_IteratorImplementation<const T>(storage, false, spIndex));
}

template <typename T>
std::set<size_t> SparseSky<T>::GetActivePixels() const noexcept
{
   std::set<size_t> pixels;
   for (auto &m : storage)
      for (auto &val : m)
         pixels.insert(val.first);

   return pixels;
}


template <typename T>
T SparseSky<T>::Accumulate() const noexcept
{
   //Kahan sum
   T sum(0);
   T c(0);

#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:sum) firstprivate(c)
   for (size_t i = 0; i < storage.size(); ++i)
      for (auto &val : storage[i]) {
         const T y = val.second - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }

   //Add the empty pixels
   size_t nEmpty(HealpixBaseExtended::Npix()*spBin.GetSize());
   for (auto &m : storage)
      nEmpty -= m.size();
   sum += empty * nEmpty;

   return sum;
}

template <typename T>
void SparseSky<T>::ApplyFunction(const std::function< void(T&) > &func) noexcept
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;

#pragma omp parallel for schedule(dynamic) default(shared)
   for (size_t i = 0; i < storage.size(); ++i)
      for (auto &val : storage[i])
         func(val.second);

}

template <typename T>
void SparseSky<T>::ApplyFunction(const std::function< void(const T&) > &func) const noexcept
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;

#pragma omp parallel for schedule(dynamic) default(shared)
   for (size_t i = 0; i < storage.size(); ++i)
      for (const auto &val : storage[i])
         func(val.second);

}

template <typename T>
void SparseSky<T>::ApplyFunction(const std::function< void(T&, size_t) > &func) noexcept
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;

#pragma omp parallel for schedule(dynamic) default(shared)
   for (size_t i = 0; i < storage.size(); ++i)
      for (auto &val : storage[i])
         func(val.second, i);

}

template <typename T>
void SparseSky<T>::ApplyFunction(const std::function< void(const T&, size_t) > &func) const noexcept
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;

#pragma omp parallel for schedule(dynamic) default(shared)
   for (size_t i = 0; i < storage.size(); ++i)
      for (const auto &val : storage[i])
         func(val.second, i);

}

template <typename T>
T SparseSky<T>::ApplyFunctionAccumulate(const std::function< T(T&) > &func) noexcept
{
   //Kahan sum
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;


#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:sum) firstprivate(c)
   for (size_t i = 0; i < storage.size(); ++i)
      for (auto &val : storage[i]) {
         const T y = func(val.second) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }

   //Add the empty pixels
   size_t nEmpty(HealpixBaseExtended::Npix()*spBin.GetSize());
   for (auto &m : storage)
      nEmpty -= m.size();

   sum += nEmpty*func(empty);
   
   return sum;
}

template <typename T>
T SparseSky<T>::ApplyFunctionAccumulate(const std::function< T(const T&) > &func) const noexcept
{
   //Kahan sum
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;

#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:sum) firstprivate(c)
   for (size_t i = 0; i < storage.size(); ++i)
      for (auto &val : storage[i]) {
         const T y = func(val.second) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }

   //Add the empty pixels
   size_t nEmpty(HealpixBaseExtended::Npix()*spBin.GetSize());
   for (const auto &m : storage)
      nEmpty -= m.size();

   sum += nEmpty*func(empty);
   
   return sum;
}

template <typename T>
T SparseSky<T>::ApplyFunctionAccumulate(const std::function< T(T&, size_t) > &func) noexcept
{
   //Kahan sum
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;

#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:sum) firstprivate(c)
   for (size_t i = 0; i < storage.size(); ++i)
      for (auto &val : storage[i]) {
         const T y = func(val.second, i) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }

   //Add the empty pixels
   const auto mapSize = HealpixBaseExtended::Npix();
   for (size_t i = 0; i < storage.size(); ++i)
      sum += func(empty, i)*(mapSize - storage[i].size());

   return sum;
}

template <typename T>
T SparseSky<T>::ApplyFunctionAccumulate(const std::function< T(const T&, size_t) > &func) const noexcept
{
   //Kahan sum
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;

#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:sum) firstprivate(c)
   for (size_t i = 0; i < storage.size(); ++i)
      for (auto &val : storage[i]) {
         const T y = func(val.second, i) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }

   //Add the empty pixels
   const auto mapSize = HealpixBaseExtended::Npix();
   for (size_t i = 0; i < storage.size(); ++i)
      sum += func(empty, i)*(mapSize - storage[i].size());

   return sum;
}


template <typename T>
void SparseSky<T>::ApplyFunction(const std::function< void(T&, const T&) > &func, const BaseSky<T> &map)
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;
   
   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

#pragma omp parallel for schedule(dynamic) default(shared)
   for (size_t i = 0; i < storage.size(); ++i) {
      //Find out the active pixels in both maps
      if (map.IsFullSky()) {
         //Combining sparse sky maps with full sky maps is incredibly stupid.
         for (int pix = 0; pix < Npix(); ++pix)
            func(storage[i][pix], map.GetValue(pix, i));

      } else {
         std::set<size_t> pixels;
         auto pit = pixels.end();
         for (auto &val : storage[i])
            pit = pixels.insert(pit, val.first); 

         const auto endIt = map.mapEnd(i);
         pit = pixels.begin(); //The first element is a good guess at the initial position
         for (auto it = map.mapBegin(i); it != endIt; ++it)
            pit = pixels.insert(pit, it.healpixIndex());

         for (auto pix : pixels)
            func(storage[i][pix], map.GetValue(pix, i));
      }
   }

}


template <typename T>
void SparseSky<T>::ApplyFunction(const std::function< void(const T&, const T&) > &func, const BaseSky<T> &map) const 
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;
   
   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

#pragma omp parallel for schedule(dynamic) default(shared)
   for (size_t i = 0; i < storage.size(); ++i) {
      //Find out the active pixels in both maps
      if (map.IsFullSky()) {

         //Use empty value unless storage is set
         //Use the fact that the pixels are sorted in storage
         auto it = storage[i].begin();
         for (size_t pix = 0; pix < size_t(Npix()); ++pix) {
            if ( (it != storage[i].end()) && pix == it->first) {
               func(it->second, map.GetValue(pix,i));
               ++it;
            } else {
               func(empty, map.GetValue(pix, i));
            }
         }

      } else {

         std::set<size_t> pixels;
         auto pit = pixels.end();
         for (auto &val : storage[i])
            pit = pixels.insert(pixels.end(), val.first); 

         const auto endIt = map.mapEnd(i);
         pit = pixels.begin(); //The first element is a good guess at the initial position
         for (auto it = map.mapBegin(i); it != endIt; ++it)
            pit = pixels.insert(pit, it.healpixIndex());

         //Use the fact that the pixels are sorted in both cases
         auto it = storage[i].begin();
         for (auto pix : pixels) {
            if ( (it != storage[i].end()) && pix == it->first) {
               func(it->second, map.GetValue(pix, i));
               ++it;
            } else {
               func(empty, map.GetValue(pix,i));
            }
         }
      }
   }

}

template <typename T>
void SparseSky<T>::ApplyFunction(const std::function< void(T&, const T&, size_t) > &func, const BaseSky<T> &map) 
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;
   
   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

#pragma omp parallel for schedule(dynamic) default(shared)
   for (size_t i = 0; i < storage.size(); ++i) {
      //Find out the active pixels in both maps
      if (map.IsFullSky()) {
         //Combining sparse sky maps with full sky maps is incredibly stupid.
         for (int pix = 0; pix < Npix(); ++pix)
            func(storage[i][pix], map.GetValue(pix, i), i);

      } else {
         std::set<size_t> pixels;
         auto pit = pixels.end();
         for (auto &val : storage[i])
            pit = pixels.insert(pit, val.first); 

         const auto endIt = map.mapEnd(i);
         pit = pixels.begin(); //The first element is a good guess at the initial position
         for (auto it = map.mapBegin(i); it != endIt; ++it)
            pit = pixels.insert(pit, it.healpixIndex());

         for (auto pix : pixels)
            func(storage[i][pix], map.GetValue(pix, i), i);
      }
   }

}


template <typename T>
void SparseSky<T>::ApplyFunction(const std::function< void(const T&, const T&, size_t) > &func, const BaseSky<T> &map) const 
{
   //Make sure we have a valid function
   if ( func == nullptr ) 
      return;
   
   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

#pragma omp parallel for schedule(dynamic) default(shared)
   for (size_t i = 0; i < storage.size(); ++i) {
      //Find out the active pixels in both maps
      if (map.IsFullSky()) {

         //Use empty value unless storage is set
         //Use the fact that the pixels are sorted in storage
         auto it = storage[i].begin();
         for (size_t pix = 0; pix < size_t(Npix()); ++pix) {
            if ( (it != storage[i].end()) && pix == it->first) {
               func(it->second, map.GetValue(pix,i), i);
               ++it;
            } else {
               func(empty, map.GetValue(pix, i), i);
            }
         }

      } else {

         std::set<size_t> pixels;
         auto pit = pixels.end();
         for (auto &val : storage[i])
            pit = pixels.insert(pit, val.first); 

         const auto endIt = map.mapEnd(i);
         pit = pixels.begin(); //The first element is a good guess at the initial position
         for (auto it = map.mapBegin(i); it != endIt; ++it)
            pit = pixels.insert(pit, it.healpixIndex());

         //Use the fact that the pixels are sorted in both cases
         auto it = storage[i].begin();
         for (auto pix : pixels) {
            if ( (it != storage[i].end()) && pix == it->first) {
               func(it->second, map.GetValue(pix, i), i);
               ++it;
            } else {
               func(empty, map.GetValue(pix,i), i);
            }
         }
      }
   }

}

template <typename T>
T SparseSky<T>::ApplyFunctionAccumulate(const std::function< T(T&, const T&) > &func, const BaseSky<T> &map) 
{
   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

   //Kahan sum
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;
   
#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:sum) firstprivate(c)
   for (size_t i = 0; i < storage.size(); ++i) {
      //Find out the active pixels in both maps
      if (map.IsFullSky()) {
         //Combining sparse sky maps with full sky maps is incredibly stupid.
         for (int pix = 0; pix < Npix(); ++pix) {
            const T y = func(storage[i][pix], map.GetValue(pix,i)) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }

      } else {

         std::set<size_t> pixels;
         auto pit = pixels.end();
         for (auto &val : storage[i])
            pit = pixels.insert(pit, val.first); 

         const auto endIt = map.mapEnd(i);
         pit = pixels.begin(); //The first element is a good guess at the initial position
         for (auto it = map.mapBegin(i); it != endIt; ++it)
            pit = pixels.insert(pit, it.healpixIndex());

         for (auto pix : pixels) {
            const T y = func(storage[i][pix], map.GetValue(pix,i)) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }

         //Add the empty pixels
         size_t nEmpty(HealpixBaseExtended::Npix() - pixels.size());

         const T y = nEmpty*func(empty, map.GetEmptyValue()) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }
   }

   return sum;
}

template <typename T>
T SparseSky<T>::ApplyFunctionAccumulate(const std::function< T(const T&, const T&) > &func, const BaseSky<T> &map) const
{
   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

   //Kahan sum
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;

#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:sum) firstprivate(c)
   for (size_t i = 0; i < storage.size(); ++i) {
      //Find out the active pixels in both maps
      if (map.IsFullSky()) {
         //Use empty value unless storage is set
         //Use the fact that the pixels are sorted in storage
         auto it = storage[i].begin();
         for (size_t pix = 0; pix < size_t(Npix()); ++pix) {
            T y;
            if ( (it != storage[i].end()) && pix == it->first) {
               y = func(it->second, map.GetValue(pix,i)) - c;
               ++it;
            } else {
               y = func(empty, map.GetValue(pix,i)) - c;
            }
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }

      } else {

         std::set<size_t> pixels;
         auto pit = pixels.end();
         for (auto &val : storage[i])
            pit = pixels.insert(pit, val.first); 

         const auto endIt = map.mapEnd(i);
         pit = pixels.begin(); //The first element is a good guess at the initial position
         for (auto it = map.mapBegin(i); it != endIt; ++it)
            pit = pixels.insert(pit, it.healpixIndex());

         //Use the fact that the pixels are sorted in both cases
         auto it = storage[i].begin();
         for (auto pix : pixels) {
            T y;
            if ( (it != storage[i].end()) && pix == it->first) {
               y = func(it->second, map.GetValue(pix,i)) - c; 
               ++it;
            } else {
               y = func(empty, map.GetValue(pix,i)) - c; 
            }
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }

         //Add the empty pixels
         size_t nEmpty(HealpixBaseExtended::Npix() - pixels.size());

         const T y = nEmpty*func(empty, map.GetEmptyValue()) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }
   }

   return sum;
}

template <typename T>
T SparseSky<T>::ApplyFunctionAccumulate(const std::function< T(T&, const T&, size_t) > &func, const BaseSky<T> &map)
{
   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

   //Kahan sum
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;
   
#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:sum) firstprivate(c)
   for (size_t i = 0; i < storage.size(); ++i) {
      //Find out the active pixels in both maps
      if (map.IsFullSky()) {
         //Combining sparse sky maps with full sky maps is incredibly stupid.
         for (int pix = 0; pix < Npix(); ++pix) {
            const T y = func(storage[i][pix], map.GetValue(pix,i), i) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }

      } else {

         std::set<size_t> pixels;
         auto pit = pixels.end();
         for (auto &val : storage[i])
            pit = pixels.insert(pit, val.first); 

         const auto endIt = map.mapEnd(i);
         pit = pixels.begin(); //The first element is a good guess at the initial position
         for (auto it = map.mapBegin(i); it != endIt; ++it)
            pit = pixels.insert(pit, it.healpixIndex());

         for (auto pix : pixels) {
            const T y = func(storage[i][pix], map.GetValue(pix,i), i) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }

         //Add the empty pixels
         size_t nEmpty(HealpixBaseExtended::Npix() - pixels.size());

         const T y = nEmpty*func(empty, map.GetEmptyValue(), i) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }
   }

   return sum;
}

template <typename T>
T SparseSky<T>::ApplyFunctionAccumulate(const std::function< T(const T&, const T&, size_t) > &func, const BaseSky<T> &map) const
{
   //The maps should be equivalent
   if ( ! Equivalent(map) )
      throw(std::invalid_argument("maps are not equivalent"));

   //Kahan sum
   T sum(0);
   T c(0);

   //Make sure we have a valid function
   if ( func == nullptr ) 
      return sum;
   
#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:sum) firstprivate(c)
   for (size_t i = 0; i < storage.size(); ++i) {
      //Find out the active pixels in both maps
      if (map.IsFullSky()) {
         //Use empty value unless storage is set
         //Use the fact that the pixels are sorted in storage
         auto it = storage[i].begin();
         for (size_t pix = 0; pix < size_t(Npix()); ++pix) {
            T y;
            if ( (it != storage[i].end()) && pix == it->first) {
               y = func(it->second, map.GetValue(pix,i), i) - c;
               ++it;
            } else {
               y = func(empty, map.GetValue(pix,i), i) - c;
            }
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }

      } else {

         std::set<size_t> pixels;
         auto pit = pixels.end();
         for (auto &val : storage[i])
            pit = pixels.insert(pit, val.first); 

         const auto endIt = map.mapEnd(i);
         pit = pixels.begin(); //The first element is a good guess at the initial position
         for (auto it = map.mapBegin(i); it != endIt; ++it)
            pit = pixels.insert(pit, it.healpixIndex());

         //Use the fact that the pixels are sorted in both cases
         auto it = storage[i].begin();
         for (auto pix : pixels) {
            T y;
            if ( (it != storage[i].end()) && pix == it->first) {
               y = func(it->second, map.GetValue(pix,i), i) - c; 
               ++it;
            } else {
               y = func(empty, map.GetValue(pix,i), i) - c; 
            }
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }

         //Add the empty pixels
         size_t nEmpty(HealpixBaseExtended::Npix() - pixels.size());

         const T y = nEmpty*func(empty, map.GetEmptyValue(), i) - c;
         const T t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }
   }

   return sum;
}


template <typename F, typename O>
void ApplyFunction(SparseSky<F> &arg1, const SparseSky<O> &other, const std::function< void(F&, const O&) > &func) 
{
   if (! arg1.Equivalent(other)) {
      ERROR("Maps should be equivalent");
      throw(std::invalid_argument("Equivalence test failed"));
   }

   //Loop over the energy planes and merge the pixel sets
#pragma omp parallel for schedule(dynamic) default(shared)
   for (size_t i = 0; i < arg1.storage.size(); ++i) {
      std::set<size_t> pixList;

      for (const auto &p : arg1.storage[i])
         pixList.insert(p.first);

      for (const auto &p : other.storage[i])
         pixList.insert(p.first);

      for (const auto &p : pixList)
         func(arg1.GetReference(p,i), other.GetValue(p,i));
   }

   //Finally call it once with the empty value
   func(arg1.GetEmptyValue(), other.GetEmptyValue());
}

template <typename F, typename O>
void ApplyFunction(const SparseSky<F> &arg1, const SparseSky<O> &other, const std::function< void(const F&, const O&) > &func) 
{
   if (! arg1.Equivalent(other)) {
      ERROR("Maps should be equivalent");
      throw(std::invalid_argument("Equivalence test failed"));
   }

   //Loop over the energy planes and merge the pixel sets
#pragma omp parallel for schedule(dynamic) default(shared)
   for (size_t i = 0; i < arg1.storage.size(); ++i) {
      std::set<size_t> pixList;

      for (const auto &p : arg1.storage[i])
         pixList.insert(p.first);

      for (const auto &p : other.storage[i])
         pixList.insert(p.first);

      for (const auto &p : pixList)
         func(arg1.GetValue(p,i), other.GetValue(p,i));
   }

   //Finally call it once with the empty value
   func(arg1.GetEmptyValue(), other.GetEmptyValue());
}

template <typename F, typename O>
F ApplyFunctionAccumulate(SparseSky<F> &arg1, const SparseSky<O> &other, const std::function< F(F&, const O&) > &func) 
{
   if (! arg1.Equivalent(other)) {
      ERROR("Maps should be equivalent");
      throw(std::invalid_argument("Equivalence test failed"));
   }

   //Kahan sum
   F sum(0);
   F c(0);

   size_t nEmpty(0);

   //Loop over the energy planes and merge the pixel sets
#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:sum,nEmpty) firstprivate(c)
   for (size_t i = 0; i < arg1.storage.size(); ++i) {
      std::set<size_t> pixList;

      for (const auto &p : arg1.storage[i])
         pixList.insert(p.first);

      for (const auto &p : other.storage[i])
         pixList.insert(p.first);

      nEmpty += arg1.Npix() - pixList.size();

      for (const auto &p : pixList) {
         const F y = func(arg1.GetReference(p,i), other.GetValue(p,i)) - c;
         const F t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }
   }

   //Finally call it once with the empty value
   sum += nEmpty * func(arg1.GetEmptyValue(), other.GetEmptyValue());

   return sum;
}

template <typename F, typename O>
F ApplyFunctionAccumulate(const SparseSky<F> &arg1, const SparseSky<O> &other, const std::function< F(const F&, const O&) > &func) 
{
   if (! arg1.Equivalent(other)) {
      ERROR("Maps should be equivalent");
      throw(std::invalid_argument("Equivalence test failed"));
   }

   //Kahan sum
   F sum(0);
   F c(0);

   size_t nEmpty(0);

   //Loop over the energy planes and merge the pixel sets
#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:sum,nEmpty) firstprivate(c)
   for (size_t i = 0; i < arg1.storage.size(); ++i) {
      std::set<size_t> pixList;

      for (const auto &p : arg1.storage[i])
         pixList.insert(p.first);

      for (const auto &p : other.storage[i])
         pixList.insert(p.first);

      nEmpty += arg1.Npix() - pixList.size();

      for (const auto &p : pixList) {
         const F y = func(arg1.GetReference(p,i), other.GetValue(p,i)) - c;
         const F t = sum + y;
         c = t - sum;
         c -= y;
         sum = t;
      }
   }

   //Finally call it once with the empty value
   sum += nEmpty * func(arg1.GetEmptyValue(), other.GetEmptyValue());

   return sum;
}

}

#endif
