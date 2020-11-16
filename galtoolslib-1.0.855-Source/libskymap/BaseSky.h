#ifndef BASE_SKY_H
#define BASE_SKY_H

#include "HealpixBaseExtended.h"
#include "Coordinate.h"

#include <Interpolation.h>
#include <ErrorLogger.h>
#include <Registry.h>

#include <healpix_map.h>
#include <tuple>
#include <functional>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <map>

#include <md5.h>

//Convenince macro to register a class with all the basic types available.  We do not register classes with other types, this is for fitsio only
#define REGISTER_SKY_CLASS(classname, name) \
   static utl::RegistryN<BaseSky<char>, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, char>::Registrar<classname<char> > charRegistrar ## name(#name); \
   static utl::RegistryN<BaseSky<short>, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, short>::Registrar<classname<short> > shortRegistrar ## name(#name); \
   static utl::RegistryN<BaseSky<int>, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, int>::Registrar<classname<int> > intRegistrar ## name(#name); \
   static utl::RegistryN<BaseSky<long>, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, long>::Registrar<classname<long> > longRegistrar ## name(#name); \
   static utl::RegistryN<BaseSky<float>, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, float>::Registrar<classname<float> > floatRegistrar ## name(#name); \
   static utl::RegistryN<BaseSky<double>, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, double>::Registrar<classname<double> > doubleRegistrar ## name(#name); \
   static utl::RegistryN<BaseSky<unsigned char>, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, unsigned char>::Registrar<classname<unsigned char> > uCharRegistrar ## name(#name); \
   static utl::RegistryN<BaseSky<unsigned short>, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, unsigned short>::Registrar<classname<unsigned short> > uShortRegistrar ## name(#name); \
   static utl::RegistryN<BaseSky<unsigned int>, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, unsigned int>::Registrar<classname<unsigned int> > uIntRegistrar ## name(#name); \
   static utl::RegistryN<BaseSky<unsigned long>, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, unsigned long>::Registrar<classname<unsigned long> > uLongRegistrar ## name(#name)

namespace SM {


class SpectralBinning {

   public:
      /** \brief Discrete point values.  sp must be ascending */
      SpectralBinning(std::vector<double>&& sp);
      SpectralBinning(const std::vector<double>& sp);

      /** \brief Integrated spectrum.  spMin and spMax must be ascending.  The bins must be adjacent. */
      SpectralBinning(std::vector<double>&& spMin, std::vector<double>&& spMax);
      SpectralBinning(const std::vector<double>& spMin, const std::vector<double>& spMax);

      //Default copy constructor

      /** \brief Returns the index of the closest energy bin */
      size_t GetIndex (double value) const noexcept;

      /** \brief True if the spectrum is integrated rather than discrete */
      bool IsIntegrated() const noexcept {return spectraMin.size() > 0;}

      /** \brief Return the size of the spectrum */
      size_t GetSize() const noexcept {return IsIntegrated() ? spectraMin.size() : spectra.size();}

      /** \brief Returns 0 length vector if spectrum is integrated */
      const std::vector<double>& GetBins() const noexcept {return spectra;}
      /** \brief Returns 0 length vector if spectrum is discrete */
      const std::vector<double>& GetBinsMin() const noexcept {return spectraMin;}
      /** \brief Returns 0 length vector if spectrum is discrete */
      const std::vector<double>& GetBinsMax() const noexcept {return spectraMax;}

      /** \brief Set the bins for a discrete spectrum.  Changes the spectrum to discrete.  Throws an error if size differs. */
      void SetBinning( const std::vector<double>& sp );
      /** \brief Set the bins for an integrated spectrum.  Changes the spectrum to integrated.  Throws an error if size differs. */
      void SetBinning( const std::vector<double>& spMin, const std::vector<double>& spMax );

      /** \brief Test for identical spectra */
      bool IsIdentical( const SpectralBinning &other) const noexcept;

   private:
      void AssertSpectra();
      std::vector<double> spectra, spectraMin, spectraMax;

};

class BaseSkyStorageIndependent : public HealpixBaseExtended {
   protected:
      SpectralBinning spBin;
      CoordSys cSys;
   public:
      BaseSkyStorageIndependent(SpectralBinning&& sp, int order, Healpix_Ordering_Scheme sc, CoordSys cs) :
         HealpixBaseExtended(order,sc),
         spBin(std::forward<SpectralBinning>(sp)),
         cSys(cs)
      {}

      BaseSkyStorageIndependent(const SpectralBinning& sp, int order, Healpix_Ordering_Scheme sc, CoordSys cs) :
         HealpixBaseExtended(order,sc),
         spBin(sp),
         cSys(cs)
      {}

      //!Access to the spectral binning
      SpectralBinning& GetBinning() noexcept {return spBin;}
      const SpectralBinning& GetBinning() const noexcept {return spBin;}
      //!Access to the coordinate system 
      CoordSys GetCoordinateSystem() const noexcept {return cSys;}

      //!Return a coordinate given a pixel
      Coordinate GetCoordinate(size_t iPix) { return Coordinate(pix2ang(iPix), cSys); }
      //!Return a pixel given a coordinate
      size_t GetPixelIndex(const Coordinate &co) { return ang2pix(co.healpixAng(cSys)); }


      /** \brief Test for equivalence of maps
       *
       * If this evaluates to true the maps have the same indexing.
       */
      bool Equivalent(const BaseSkyStorageIndependent &other) const {
         bool equiv = spBin.IsIdentical(other.spBin);
         equiv &= conformable(other);
         equiv &= cSys == other.cSys;
         return equiv;
      }

      //! Return the healpix order of the map
      virtual size_t GetHealpixOrder() const noexcept {return Order();}

};


/** Base class defining the interface to Skymaps
 *
 * Some restrictions on the maps
 * * The HEALPix grid is the same for all energies
 * * The spectral grid is the same for all pixels
 * * The type T must have basic arithmetic (+,-,*,/) defined between its type and the basic types of integers and floats.
 *
 * Allows for access to arbitrary underlying storage mechanism, both continuous and discrete.
 * There is no guarantee that memory has been reserved for all pixels and spectral bins.
 *
 * Access through the iterators is preferred.
 * The iterators provide access to the HEALPix index
 * and the energy index.
 * Modifications through the function interface is even more preferred.
 * Pixels with no data (no iterator pointing to it) return
 * an empty value in spectra and maps, empty values are
 * implementation defined and set in the derived classes.
 *
 */
template < typename T >
class BaseSky : public BaseSkyStorageIndependent {
   protected:
      //We only support bidirectional iterators
      template < typename C > 
      class IteratorImplementation {
         public:
            virtual ~IteratorImplementation() {}
            virtual IteratorImplementation<C>* clone() const = 0;
            virtual void increment() = 0;
            virtual void decrement() = 0;
            virtual C& dereference() const = 0;
            virtual size_t healpixIndex() const = 0;
            virtual size_t energyIndex() const = 0;
            virtual bool equal(const IteratorImplementation<C>&) const = 0;
      };

      //Template away the constantness
      template < typename C > 
      class T_iterator {
            IteratorImplementation<C> *impl;
         public:
            //Memory management for the implementation
            T_iterator() : impl() {}
            T_iterator(IteratorImplementation<C> *im) : impl(im) {}
            T_iterator(T_iterator<C> const &right) : impl(right.impl->clone()) {}
            ~T_iterator() { delete impl; }
            T_iterator<C>& operator= (T_iterator<C> const& right) {
               if (&right != this) {
                  delete impl;
                  impl = right.impl->clone();
               }
               return *this;
            }

            //using difference_type = size_t; //Not used. 
            using value_type = T;
            using pointer = C*;
            using reference = C&;
            using iterator_category = std::bidirectional_iterator_tag;

            //Comparison operators
            bool operator == (const T_iterator<C> & other) const { return impl->equal(*other.impl); }
            bool operator != (const T_iterator<C> & other) const { return  !(*this == other); }

            //Indices
            size_t healpixIndex() const { return impl->healpixIndex(); }
            size_t energyIndex() const { return impl->energyIndex(); }

            //Dereference
            reference operator*() const { return impl->dereference(); }

            //Increment and decrement
            T_iterator<C>& operator++ () {
               impl->increment();
               return *this;
            }
            T_iterator<C>& operator-- () {
               impl->decrement();
               return *this;
            }
            T_iterator<C> operator++ (int) {
               auto newit = T_iterator<C>(*this);
               impl->increment();
               return *this;
            }
            T_iterator<C> operator-- (int) {
               auto newit = T_iterator<C>(*this);
               impl->decrement();
               return *this;
            }
      };

   public:

      //! Use the registry to create different types of maps.
      static std::unique_ptr< BaseSky<T> > create(const std::string &name, const SpectralBinning& sp, int order, Healpix_Ordering_Scheme sc, CoordSys cs, const T &null);
      static std::unique_ptr< BaseSky<T> > create(const std::string &name, const BaseSkyStorageIndependent &templ, const T &null);

      /** \brief Iterator definitions
       *
       * The iterators return a tuple when dereferenced with a 
       * reference to the pixel value, the HEALPix index, and the energy index
       * in that order.
       */
      using iterator = T_iterator<T>;
      using const_iterator = T_iterator<const T>;

      /** \brief Constructor specifying everything
       *
       * The spectral binning, order, scheme, and coordsys cannot be changed without
       * destroying the map contents.
       *
       * null is used as the default value for the map.
       */
      BaseSky(SpectralBinning&& sp, int order, Healpix_Ordering_Scheme sc, CoordSys cs, const T& null = T(0)) :
         BaseSkyStorageIndependent(std::forward<SpectralBinning>(sp), order, sc, cs),
         empty(null)
      {}
      BaseSky(const SpectralBinning& sp, int order, Healpix_Ordering_Scheme sc, CoordSys cs, const T& null = T(0)) :
         BaseSkyStorageIndependent(sp, order, sc, cs),
         empty(null)
      {}
      BaseSky(const BaseSkyStorageIndependent &templ, const T& null) :
         BaseSkyStorageIndependent(templ),
         empty(null)
      {}


      virtual ~BaseSky() {}

      //Use default copy constructors
      
      // Classes should be able to clone themselves,
      virtual std::unique_ptr< BaseSky<T> > clone() const = 0;

      /** \brief Reset the map to a new size and binning
       *
       * This destroys the contents of the map.
       */
      virtual void Reset(const SpectralBinning& sp, int order, Healpix_Ordering_Scheme sc, CoordSys cs, const T& null = T(0)) = 0;
      void Reset(const BaseSkyStorageIndependent &templ, const T& null);

      //! Return the name of the object used in the registry
      virtual std::string name() const noexcept = 0;

      //! Constant reference to the value of empty pixels, only relevant for maps that don't store the full sky
      const T& GetEmptyValue() const noexcept { return empty; }
      T& GetEmptyValue() noexcept { return empty; }

      //! Set the default empty value, only relevant for maps that don't store the full sky
      void SetEmptyValue(T value) noexcept { empty = value; }

      /** \brief Iterator to the beginning of the map
       *
       * Loops through all pixels and energy range (no order is guaranteed)
       *
       * Dereferencing an iterator gives a tuple with (value reference, energy index, healpix index)
       */
      virtual iterator begin() noexcept = 0;
      virtual const_iterator begin() const noexcept = 0;
      //! Iterator to the end of the map
      virtual iterator end() noexcept = 0;
      virtual const_iterator end() const noexcept = 0;

      //! Constant iterator to the beginning of the map
      virtual const_iterator cbegin() const noexcept = 0;
      //! Constant iterator to the end of the map
      virtual const_iterator cend() const noexcept = 0;

      /** \brief Return the value at a given coordinate and energy
       *
       * Returns \par empty if value not set
       * Uses nearest neighbor search
       */
      const T& GetValue(const SM::Coordinate &co, const double energy) const noexcept;

      /** \brief Return the value at a given healpix index and spectral index
       *
       * Returns \par empty if value not set
       */
      virtual const T& GetValue(const size_t hpIndex, const size_t spIndex) const noexcept = 0;

      /** \brief Set the value at a given coordinate and energy
       *
       * Uses nearest neighbor search
       */
      void SetValue(const SM::Coordinate &co, const double energy, const T& value) noexcept;

      /** \brief Set the value at a given healpix index and spectral index
       *
       */
      virtual void SetValue(const size_t hpIndex, const size_t spIndex, const T& value) noexcept = 0;

      /** \brief Return a reference to the value at a given coordinate and energy
       *
       * Creates the value if it is not set and initializes to \par empty
       * Uses nearest neighbor search
       */
      T& GetReference(const SM::Coordinate &co, const double energy) noexcept;

      /** \brief Return a reference to the value at a given healpix index and spectral index
       *
       * Creates the value if it is not set and initializes to \par empty
       */
      virtual T& GetReference(const size_t hpIndex, const size_t spIndex) noexcept = 0;

      /** \brief Return the spectrum for a given coordinate 
       * 
       * Selects the closest pixel only.
       */
      virtual std::vector<T> operator [] (const SM::Coordinate &co) const noexcept;

      /** \brief Return the spectrum for a given coordinate 
       *
       * Selects the closest pixel.
       * The sp vector is resized to match the spectrum of the map.
       */
      void GetSpectrum (const SM::Coordinate &co, std::vector<T> &sp) const noexcept;

      /** \brief Return the spectrum for a given healpix index
       * 
       */
      virtual std::vector<T> operator [] (const size_t hpIndex) const noexcept;

      /** \brief Return the spectrum for a given healpix index
       *
       * The sp vector is resized to match the spectrum of the map.
       */
      virtual void GetSpectrum (const size_t hpIndex, std::vector<T> &sp) const noexcept;

      /** \brief Return the spectra for a given coordinate with sky interpolation
       *
       * Uses linear interpolation on the sky.
       * Assumes the value of empty for pixels without value.
       */
      std::vector<T> GetSkyInterpolatedSpectra(const SM::Coordinate &co) const noexcept;

      /** \brief Return a HEALPix map at the given energy
       *
       * No interpolation is performed, selects the nearest neighbor.
       */
      Healpix_Map<T> GetHealpixMap(double energy) const noexcept;

      /** \brief Return a HEALPix map at the given spectral index
       *
       */
      virtual Healpix_Map<T> GetHealpixMap(const size_t spIndex) const noexcept;

      /** \brief Convert to another coordinate system
       *
       * Creates a new map in each case, even if coordinate system is the same.
       */
      std::unique_ptr< BaseSky<T> > CoordinateConversion(CoordSys newCoord) const;

      /** \brief Fill input map with interpolated values from current map
       *
       * Uses linear interpolation over sky coordinates and power-law interpolation
       * for the spectrum.
       *
       * If \param extrapolate is true we use extrapolation for spectral values that are
       * out of bounds, otherwise we set those to empty value.
       *
       * Uses geometric mean of spectral boundaries for maps with binned spectrum
       */
      void Interpolate( BaseSky<T> &map, bool extrapolate ) const;

      /** \brief Return a new map with new order but otherwise identical
       *
       * If solidAngleCorrection is false then map is assumed to be in units 
       * sr^-1 (intensity), otherwise we take into account the different solid
       * angle of the input and output.
       *
       * A new map is always created, even if newOrder is the same as current order.
       */
      std::unique_ptr<BaseSky<T> > ResetOrder( int newOrder, bool solidAngleCorrection ) const;

      /** \brief Sum up all the values
       *
       * The type T has to be appropriate for openmp reduction clause.
       * Empty value is added nEmpty times to the sum where nEmpty is
       * the number of empty pixels in the map.
       */
      virtual T Accumulate() const noexcept = 0;
      
      /** \brief Apply function to all values in map
       *
       * The function is only called for active values in the map.
       * It completely ignores pixels that are empty.
       *
       * \param func should be threadsafe
       */
      virtual void ApplyFunction( const std::function< void(T&) > &func ) noexcept = 0;
      virtual void ApplyFunction( const std::function< void(const T&) > &func ) const noexcept = 0;

      /** \brief Apply function to all values in map and passes the spectral index
       *
       * The function is only called for active values in the map.
       * It completely ignores pixels that are empty.
       *
       * \param func should be threadsafe.
       */
      virtual void ApplyFunction( const std::function< void(T&, size_t) > &func ) noexcept = 0;
      virtual void ApplyFunction( const std::function< void(const T&, size_t) > &func ) const noexcept = 0;

      /** \brief Apply function to all values in both maps
       *
       * The function is called only for active values in both maps.
       * empty pixels are ignored.
       *
       * The map should be equivalent to current map.
       *
       * \param func should be threadsafe
       */
      virtual void ApplyFunction( const std::function< void(T&, const T&) > &func, const BaseSky<T> &map ) = 0;
      virtual void ApplyFunction( const std::function< void(const T&, const T&) > &func, const BaseSky<T> &map ) const = 0;
      virtual void ApplyFunction( const std::function< void(T&, const T&, size_t) > &func, const BaseSky<T> &map ) = 0;
      virtual void ApplyFunction( const std::function< void(const T&, const T&, size_t) > &func, const BaseSky<T> &map ) const = 0;

      /** \brief Apply function to all values in map and accumulate
       *
       * For this to work the type T has to be usable in openmp reduce clause
       * \param func should be threadsafe
       *
       * The function is called once with the 
       * empty value and the result added nEmpty times to the sum 
       * where nEmpty is the number of empty pixels in the map.
       */
      virtual T ApplyFunctionAccumulate( const std::function< T(T&) > &func ) noexcept = 0;
      virtual T ApplyFunctionAccumulate( const std::function< T(const T&) > &func ) const noexcept = 0;
      
      /** \brief Similar as above but with spectral index
       *
       * The function is called once for each energy bin 
       * with the empty value, starting at the lowest energy bin.
       * The value in each energy bin is added nEmpty_i times to
       * the sum where nEmpty_i is the number of empty pixels 
       * in each energy bin.
       */
      virtual T ApplyFunctionAccumulate( const std::function< T(T&, size_t) > &func ) noexcept = 0;
      virtual T ApplyFunctionAccumulate( const std::function< T(const T&, size_t) > &func ) const noexcept = 0;

      /** \brief Similar as above, but include other map
       *
       * The map should be equivalent to the current map.
       *
       */
      virtual T ApplyFunctionAccumulate( const std::function< T(T&, const T&) > &func, const BaseSky<T> &map ) = 0;
      virtual T ApplyFunctionAccumulate( const std::function< T(const T&, const T&) > &func, const BaseSky<T> &map ) const = 0;
      virtual T ApplyFunctionAccumulate( const std::function< T(T&, const T&, size_t) > &func, const BaseSky<T> &map ) = 0;
      virtual T ApplyFunctionAccumulate( const std::function< T(const T&, const T&, size_t) > &func, const BaseSky<T> &map ) const = 0;

      //! Returns true if the underlying storage is full-sky.
      virtual bool IsFullSky() const noexcept = 0;

      /** \brief Returns the indices of the active pixels in the map
       * 
       * This should list all pixels that have non-empty values in at least
       * one of the energy planes. 
       */
      virtual std::set<size_t> GetActivePixels() const noexcept = 0;

      //! Calculate an md5 checksum, using the data and binning information
      std::string md5() const noexcept;

      //! Return iterators for the spectrum given the index
      virtual iterator spectraBegin(const size_t hpIndex) noexcept = 0;
      virtual iterator spectraEnd(const size_t hpIndex) noexcept = 0;
      virtual const_iterator spectraBegin(const size_t hpIndex) const noexcept = 0;
      virtual const_iterator spectraEnd(const size_t hpIndex) const noexcept = 0;

      //! Return iterators for a map given a spectral index
      virtual iterator mapBegin(const size_t spIndex) noexcept = 0;
      virtual iterator mapEnd(const size_t spIndex) noexcept = 0;
      virtual const_iterator mapBegin(const size_t spIndex) const noexcept = 0;
      virtual const_iterator mapEnd(const size_t spIndex) const noexcept = 0;

   protected:
      T empty;
};


//Implementations of functions follow
template < typename T>
std::unique_ptr<BaseSky<T> > BaseSky<T>::create(const std::string &name, const SpectralBinning& sp, int order, Healpix_Ordering_Scheme sc, CoordSys cs, const T &null)
{
   return utl::RegistryN<BaseSky<T>, const SpectralBinning&, int, Healpix_Ordering_Scheme, CoordSys, T>::create(name, sp, order, sc, cs, null);
}

template < typename T>
std::unique_ptr< BaseSky<T> > BaseSky<T>::create(const std::string &name, const BaseSkyStorageIndependent &templ, const T &null) {
   return create(name, templ.GetBinning(), templ.Order(), templ.Scheme(), templ.GetCoordinateSystem(), null);
}

template < typename T >
void BaseSky<T>::Reset(const BaseSkyStorageIndependent &templ, const T &null) {
   Reset(templ.GetBinning(), templ.Order(), templ.Scheme(), templ.GetCoordinateSystem(), null);
}

template < typename T >
const T& BaseSky<T>::GetValue(const SM::Coordinate &co, const double energy) const noexcept
{
   return GetValue(coord2pix(co.healpixAng(cSys)), spBin.GetIndex(energy));
}

template < typename T >
void BaseSky<T>::SetValue(const SM::Coordinate &co, const double energy, const T& value) noexcept
{
   SetValue(coord2pix(co.healpixAng(cSys)), spBin.GetIndex(energy), value);
}

template < typename T >
T& BaseSky<T>::GetReference(const SM::Coordinate &co, const double energy) noexcept
{
   return GetReference(coord2pix(co.healpixAng(cSys)), spBin.GetIndex(energy));
}

template < typename T >
std::vector<T> BaseSky<T>::operator [] (const SM::Coordinate &co) const noexcept
{
   return operator[] (coord2pix(co.healpixAng(cSys)));
}

template < typename T >
void BaseSky<T>::GetSpectrum (const SM::Coordinate &co, std::vector<T> &sp) const noexcept
{
   GetSpectrum(coord2pix(co.healpixAng(cSys)), sp);
}

template < typename T >
std::vector<T> BaseSky<T>::operator [] (const size_t hpIndex) const noexcept
{
   std::vector<T> sp;

   GetSpectrum(hpIndex, sp);

   return sp;
}

template < typename T >
void BaseSky<T>::GetSpectrum (const size_t hpIndex, std::vector<T> &sp) const noexcept
{
   sp.assign(spBin.GetSize(), empty);

   const auto endIt = spectraEnd(hpIndex);
   for (auto it = spectraBegin(hpIndex); it != endIt; ++it) 
      sp[it.energyIndex()] = *it;
}

template < typename T >
std::vector<T> BaseSky<T>::GetSkyInterpolatedSpectra(const SM::Coordinate &co) const noexcept
{
   //For storing the output
   std::vector<T> sp(spBin.GetSize(), 0);
   //For storing the values in each pixel
   std::vector<T> spPix(spBin.GetSize());

   //Use the interpol function from Healpix
   fix_arr<int, 4> pixels;
   fix_arr<double, 4> weight;

   get_interpol(co.healpixAng(cSys), pixels, weight);

   for ( size_t i(0); i < 4; ++i ) {
      GetSpectrum(pixels[i], spPix);
      for ( size_t j(0); j < spPix.size(); ++j)
         sp[j] += spPix[j]*weight[i];
   }

   return sp;

}

template < typename T >
Healpix_Map<T> BaseSky<T>::GetHealpixMap(double energy) const noexcept
{
   const size_t index = spBin.GetIndex(energy);
   return GetHealpixMap(index);
}

template < typename T >
Healpix_Map<T> BaseSky<T>::GetHealpixMap(size_t spind) const noexcept
{
   Healpix_Map<T> out(Order(), Scheme());

   if (! IsFullSky())
      out.fill(empty);

   auto it = mapBegin(spind);
   const auto endit = mapEnd(spind); 
   for ( ; it != endit; ++it) {
      out[it.healpixIndex()] = *it;
   }

   return out;
}

template < typename T >
std::unique_ptr< BaseSky<T> > BaseSky<T>::CoordinateConversion(CoordSys newCoord) const 
{

   //Create the new map, use the same type of map
   //Have empty value as 0
   auto newMap = BaseSky<T>::create(name(), GetBinning(), Order(), Scheme(), newCoord, T(0));

   if (newCoord == cSys) {

      //Simply copy the values to the new map
      auto it = begin();
      const auto endit = end();

      while (it != endit) {
         newMap->SetValue(it.healpixIndex(), it.energyIndex(), *it);
         ++it;
      }

   } else {

      //Use brute force conversion, loop over the active pixels at higher resolution and
      //insert into the new map.  Use Coordinate for the transformation

      //Use 5 orders finer map (or up to 13)
      const int dOrder = std::min(13-order_, 5);
      const int nPix = (1<<dOrder)*(1<<dOrder);
      
      //For the pointing to use in the coordinate transform
      Healpix_Base hpb(order_+dOrder, NEST);

      std::vector<size_t> vecPixels;

      if (IsFullSky()) {
         //All pixels
         vecPixels.resize(Npix());
         for (size_t hpi = 0; hpi < size_t(Npix()); ++hpi)
            vecPixels[hpi] = hpi;
         
      } else {

         // Need to find the effected pixel list
         std::set<size_t> pixels;

         //Loop over pixels in current map and find neighbours of that direction in input
         const auto active = GetActivePixels();

         const double radius = resolution()*utl::kConvertDegreesToRadians;
         for ( auto pix : active ) {
            //Use coordinates to swap between coordinate systems
            Coordinate co(pix2ang(pix), cSys);

            rangeset<int> newPixels;
            newMap->query_disc_inclusive(co.healpixAng(newMap->GetCoordinateSystem()), radius, newPixels);

            for ( int np : newPixels.toVector() )
               pixels.insert(np);
         }

         //To allow for parallelization, create all the pixels with 0 value
         for (auto pix : pixels)
            for (size_t j = 0; j < spBin.GetSize(); ++j) 
               newMap->GetReference(pix, j) = 0;

         //Turn pixel list into vectors so we can use openmp
         vecPixels.reserve(pixels.size());
         vecPixels.insert(vecPixels.begin(), pixels.begin(), pixels.end());

      }

      //Loop over the affected pixels and calculate the interpolation

#pragma omp parallel for schedule(static) default(shared)
      for (size_t i = 0; i < vecPixels.size(); ++i) {

         const int nestPix = Scheme() == NEST ? vecPixels[i] : ring2nest(vecPixels[i]);

         for (int p = nestPix*nPix; p < (nestPix+1)*nPix; ++p) {

            Coordinate co(hpb.pix2ang(p),newCoord);
            auto newPix = ang2pix(co.healpixAng(cSys));

            for (size_t j = 0; j < spBin.GetSize(); ++j) {
               newMap->GetReference(vecPixels[i], j) += GetValue(newPix, j);
            }

         }

         for (size_t j = 0; j < spBin.GetSize(); ++j) {
            newMap->GetReference(vecPixels[i], j) /= nPix;
         }
      }

   }

   newMap->SetEmptyValue(GetEmptyValue());

   return newMap;
}


template < typename T >
std::string BaseSky<T>::md5() const noexcept
{
   MD5 cs;
   for (auto val : *this)
      cs.update(reinterpret_cast<const unsigned char*>(&val), sizeof(T)/sizeof(char));

   cs.update(reinterpret_cast<const unsigned char*>(&GetBinning().GetBins()[0]), GetBinning().GetBins().size()*sizeof(double)/sizeof(char));
   cs.update(reinterpret_cast<const unsigned char*>(&GetBinning().GetBinsMin()[0]), GetBinning().GetBinsMin().size()*sizeof(double)/sizeof(char));
   cs.update(reinterpret_cast<const unsigned char*>(&GetBinning().GetBinsMax()[0]), GetBinning().GetBinsMax().size()*sizeof(double)/sizeof(char));
   cs.finalize();
   return cs.hexdigest();
}


template < typename T >
void BaseSky<T>::Interpolate(BaseSky<T> &map, bool extrapolate) const 
{
   // Fill map with empty value
   map.ApplyFunction([=](T&a){ a = empty; });

   // Need to find the affected pixel list
   std::set<size_t> pixels;

   if ( IsFullSky() ) {
      //All pixels
      for (size_t i=0; i < size_t(map.Npix()); ++i)
         pixels.insert(i);

   } else {
      //Loop over pixels in current map and find neighbours of that direction in input
      const auto active = GetActivePixels();

      const double radius = resolution()*utl::kConvertDegreesToRadians;
      for ( auto pix : active ) {
         //Use coordinates to swap between coordinate systems
         Coordinate co(pix2ang(pix), cSys);

         rangeset<int> newPixels;
         map.query_disc_inclusive(co.healpixAng(map.GetCoordinateSystem()), radius, newPixels);

         for ( auto np : newPixels.toVector() )
            pixels.insert(np);
      }

      //Insert values to the pixels and therefore create them so we can do openmp later
      for ( auto pix : pixels )
         for ( size_t i = 0; i < map.GetBinning().GetSize(); ++i )
            map.SetValue(pix, i, empty);
   }

   //Loop over the pixels and interpolate
   //Turn the set into vector for openmp parallelization
   std::vector<size_t> vecPixels(pixels.begin(), pixels.end());

   std::vector<double> inEnergies(spBin.GetSize()), outEnergies(map.GetBinning().GetSize());

   if (spBin.GetBinsMin().size() == 0) {
      inEnergies = spBin.GetBins();
   } else {
      for (size_t i = 0; i < spBin.GetSize(); ++i)
         inEnergies[i] = sqrt(spBin.GetBinsMin()[i]*spBin.GetBinsMax()[i]);
   }

   if (map.GetBinning().GetBinsMin().size() == 0) {
      outEnergies = map.GetBinning().GetBins();
   } else {
      for (size_t i = 0; i < map.GetBinning().GetSize(); ++i)
         outEnergies[i] = sqrt(map.GetBinning().GetBinsMin()[i]*map.GetBinning().GetBinsMax()[i]);
   }

   // Only do interpolation if the spectral bins do not match
   bool doInterpolation = true;
   if (inEnergies.size() == outEnergies.size()) {
      doInterpolation = false;
      for (size_t i = 0; i < inEnergies.size(); ++i) {
         if ( fabs(inEnergies[i] + outEnergies[i]) > 0 ) {
            if ( fabs( (inEnergies[i] - outEnergies[i])*2/(inEnergies[i]+outEnergies[i]) ) > 1e-8 ) {
               doInterpolation = true;
               break;
            }
         }
         else
         {
            if ( fabs(inEnergies[i] - outEnergies[i]) > 1e-10 ) {
               doInterpolation = true;
               break;
            }
         }
      }
   }

   if ( doInterpolation ) {

      // Pre-calculate the valid range for interpolation if extrapolation is turned off.
      size_t iBegin = 0, iEnd = outEnergies.size();
      if ( !extrapolate ) {
         size_t i(0);
         while (i < outEnergies.size() && outEnergies[i] < inEnergies.front())
            ++i;
         iBegin = i;

         i = outEnergies.size();
         while (i > 0 && outEnergies[i-1] > inEnergies.back())
            --i;
         iEnd = i;
      }

      // Do power-law interpolation by default, but revert to linear if 
      // the spectrum (unlikely) or any value is 0
      bool linearInterpolation = false;
      for (auto e : inEnergies)
         if ( e <= 0 )
            linearInterpolation = true;
      for (auto e : outEnergies)
         if ( e <= 0 )
            linearInterpolation = true;

      //Need to convert all types to double before interpolation
      std::vector<double> inSpectra(inEnergies.size());

      //Do interpolation (and extrapolation)
#pragma omp parallel for schedule(static) default(shared) firstprivate(inSpectra)
      for (size_t i = 0; i < vecPixels.size(); ++i) {
         auto intSpectra = GetSkyInterpolatedSpectra(Coordinate(map.pix2ang(vecPixels[i]), map.cSys));
         for ( size_t ie = 0; ie < intSpectra.size(); ++ie)
            inSpectra[ie] = double(intSpectra[ie]);

         bool doLin = linearInterpolation;
         if ( !linearInterpolation ) {
            for ( auto e : inSpectra ) {
               if ( e <= 0) {
                  doLin = true;
                  break;
               }
            }
         }

         std::unique_ptr<utl::Interpolation1D> interpolator;
         if (doLin)
            interpolator.reset( new utl::LinearInterpolation(inEnergies, inSpectra));
         else
            interpolator.reset( new utl::PowerLawInterpolation(inEnergies, inSpectra));

         for ( size_t ie = iBegin; ie < iEnd; ++ie ) 
            map.SetValue(vecPixels[i], ie, (*interpolator)(outEnergies[ie]));

      }
      
   } else {

#pragma omp parallel for schedule(static) default(shared)
      for (size_t i = 0; i < vecPixels.size(); ++i) {
         auto intSpectra = GetSkyInterpolatedSpectra(Coordinate(map.pix2ang(vecPixels[i]), map.cSys));

         for ( size_t ie = 0; ie < intSpectra.size(); ++ie ) 
            map.SetValue(vecPixels[i], ie, intSpectra[ie]);
      }
   }
}


template < typename T >
std::unique_ptr< BaseSky<T> > BaseSky<T>::ResetOrder( int newOrder, bool solidAngleCorrection ) const 
{

  //Create the new map, use the same type of map
  //Have empty value as 0
  auto newMap = BaseSky<T>::create(name(), GetBinning(), newOrder, Scheme(), GetCoordinateSystem(), T(0));
  
  //What happens depends on the order of the map
  if (Order() == newOrder) {
    //Copy all set pixels to the new map
    newMap->ApplyFunction([](T& a, const T& b){ a = b; }, *this);
  }
  
  if (Order() < newOrder) {
    //Set the value of pixels to the underlying pixel
    //Loop over the pixels in the current map and pick the pixels underneath
    const auto endit = end();
    for (auto it = begin(); it != endit; ++it) {
      std::vector<int> pixels;
      newMap->query_pixel(*this, it.healpixIndex(), pixels);
      
      //Apply solid angle correction
      const T value = solidAngleCorrection ? T(*it/double(pixels.size())) : *it;
      for (auto pix : pixels)
	newMap->SetValue(pix, it.energyIndex(), value);
    }
    
  }
  
  if (Order() > newOrder) {
    //Need to find which pixels are affected
    std::set<size_t> pixels;
    
    if (IsFullSky()) {
      for (size_t i = 0; i < size_t(newMap->Npix()); ++i)
	pixels.insert(i);
      
    } else {
      //Loop over the set pixels and add the one associated with it
      for (auto iPix : GetActivePixels())
	pixels.insert(newMap->ang2pix(pix2ang(iPix)));
    }
    
    //Loop over the pixel and use the average value underneath
    for (auto pix : pixels) {
      std::vector<int> listPix;
      query_pixel(*newMap, pix, listPix);
      for (auto iPix : listPix)
	for (size_t ie = 0; ie < spBin.GetSize(); ++ie) 
	  newMap->GetReference(pix,ie) += GetValue(iPix, ie);
      //Average, if needed
      if (! solidAngleCorrection) 
	for (size_t ie = 0; ie < spBin.GetSize(); ++ie) 
	  newMap->GetReference(pix,ie) /= double(listPix.size());
    }
    
  }

  //Set the correct empty value
  newMap->SetEmptyValue(empty);
    
  return newMap;
  
}

/** \brief Apply function to all values in maps
 *
 * The maps have to be equivalent, otherwise we throw an exception.
 * \param func should be threadsafe
 *
 * Default implementation uses GetReference for arg1 and GetValue for arg2, specialized 
 * implementation may be faster.
 *
 * Uses only the common subset of active pixels and calls func once with 
 * empty value of both if not all pixels are used
 */
template <typename T, typename O>
void ApplyFunction(BaseSky<T> &arg1, const BaseSky<O> &arg2, const std::function< void(T&, const O&) > &func)
{
   if (! arg1.Equivalent(arg2)) {
      ERROR("Maps should be equivalent");
      throw(std::invalid_argument("Equivalence test failed"));
   }

   //This cannot be parallel because GetReference is not guaranteed to be thread safe
   if ( arg1.IsFullSky() || arg2.IsFullSky() ) {
      for (size_t hpi = 0; hpi < size_t(arg1.Npix()); ++hpi)
         for (size_t ei = 0; ei < arg1.GetBinning().GetSize(); ++ei)
            func(arg1.GetReference(hpi,ei), arg2.GetValue(hpi,ei));
   } else {
      auto pix1 = arg1.GetActivePixels();
      auto pix2 = arg2.GetActivePixels();

      pix1.insert(pix2.begin(),pix2.end());
      for (auto hpi : pix1)
         for (size_t ei = 0; ei < arg1.GetBinning().GetSize(); ++ei)
            func(arg1.GetReference(hpi,ei), arg2.GetValue(hpi,ei));

      if (pix1.size() < size_t(arg1.Npix()))
         func(arg1.GetEmptyValue(), arg2.GetEmptyValue());
   }

}

/** \brief Apply function to all values in maps
 *
 * The maps have to be equivalent, otherwise we throw an exception.
 * \param func should be threadsafe
 *
 * Default implementation uses GetValue for both maps, specialized 
 * implementation may be faster.
 *
 * Uses only the common subset of active pixels and calls func once with 
 * empty value of both if not all pixels are used
 */
template <typename T, typename O>
void ApplyFunction(const BaseSky<T> &arg1, const BaseSky<O> &arg2, const std::function< void(const T&, const O&) > &func)
{
   if (! arg1.Equivalent(arg2)) {
      ERROR("Maps should be equivalent");
      throw(std::invalid_argument("Equivalence test failed"));
   }

   if ( arg1.IsFullSky() || arg2.IsFullSky() ) {
#pragma omp parallel for default(shared) schedule(static)
      for (size_t hpi = 0; hpi < size_t(arg1.Npix()); ++hpi)
         for (size_t ei = 0; ei < arg1.GetBinning().GetSize(); ++ei)
            func(arg1.GetValue(hpi,ei), arg2.GetValue(hpi,ei));
   } else {
      auto pix1 = arg1.GetActivePixels();
      auto pix2 = arg2.GetActivePixels();

      pix1.insert(pix2.begin(),pix2.end());
#pragma omp parallel default(shared)
      for (auto hpi : pix1)
#pragma omp single nowait
         for (size_t ei = 0; ei < arg1.GetBinning().GetSize(); ++ei)
            func(arg1.GetValue(hpi,ei), arg2.GetValue(hpi,ei));

      if (pix1.size() < size_t(arg1.Npix()))
         func(arg1.GetEmptyValue(), arg2.GetEmptyValue());
   }

}

/** \brief Apply function to all values in maps
 *
 * For this to work the type T has to be usable in openmp reduce clause
 * \param func should be threadsafe
 *
 * Default implementation uses GetReference for arg1 and GetValue for arg2. Specialized 
 * implementation may be faster.
 *
 * Uses only the common subset of active pixels and calls func once with 
 * empty value of both if not all pixels are used and adds the returned 
 * value multiplied with the number of pixels.
 */
template <typename T, typename O>
T ApplyFunctionAccumulate(BaseSky<T> &arg1, const BaseSky<O> &arg2, const std::function< T(T&, const O&) > &func)
{
   if (! arg1.Equivalent(arg2)) {
      ERROR("Maps should be equivalent");
      throw(std::invalid_argument("Equivalence test failed"));
   }

   //Use kahan summing
   T sum(0);
   T c(0);

   const bool isFull1 = arg1.IsFullSky();
   if ( isFull1 || arg2.IsFullSky() ) {
      //This cannot be parallel unless arg1 is fullsky because GetReference is not thread safe otherwise
#pragma omp parallel for default(shared) schedule(static) firstprivate(c), reduction(+:sum) if (isFull1)
      for (size_t hpi = 0; hpi < size_t(arg1.Npix()); ++hpi)
         for (size_t ei = 0; ei < arg1.GetBinning().GetSize(); ++ei) {
            const T y = func(arg1.GetReference(hpi,ei), arg2.GetValue(hpi,ei)) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }
   } else {
      auto pix1 = arg1.GetActivePixels();
      auto pix2 = arg2.GetActivePixels();

      pix1.insert(pix2.begin(),pix2.end());
      for (auto hpi : pix1)
         for (size_t ei = 0; ei < arg1.GetBinning().GetSize(); ++ei) {
            const T y = func(arg1.GetReference(hpi,ei), arg2.GetValue(hpi,ei)) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }

      const size_t nEmpty = arg1.Npix() - pix1.size();
      if (nEmpty > 0)
         sum += nEmpty * arg1.GetBinning().GetSize() * func(arg1.GetEmptyValue(), arg2.GetEmptyValue());
   }

   return sum;
}

/** \brief Apply function to all values in maps
 *
 * For this to work the type T has to be usable in openmp reduce clause
 * \param func should be threadsafe
 *
 * Default implementation uses GetValue for both maps, specialized 
 * implementation may be faster and more convenient.
 */
template <typename T, typename O>
T ApplyFunctionAccumulate(const BaseSky<T> &arg1, const BaseSky<O> &arg2, const std::function< T(const T&, const O&) > &func)
{
   if (! arg1.Equivalent(arg2)) {
      ERROR("Maps should be equivalent");
      throw(std::invalid_argument("Equivalence test failed"));
   }

   //Use kahan summing
   T sum(0);
   T c(0);

   if ( arg1.IsFullSky() || arg2.IsFullSky() ) {
#pragma omp parallel for default(shared) schedule(static) firstprivate(c), reduction(+:sum)
      for (size_t hpi = 0; hpi < size_t(arg1.Npix()); ++hpi)
         for (size_t ei = 0; ei < arg1.GetBinning().GetSize(); ++ei) {
            const T y = func(arg1.GetValue(hpi,ei), arg2.GetValue(hpi,ei)) - c;
            const T t = sum + y;
            c = t - sum;
            c -= y;
            sum = t;
         }
   } else {
      auto pix1 = arg1.GetActivePixels();
      auto pix2 = arg2.GetActivePixels();

      pix1.insert(pix2.begin(),pix2.end());
#pragma omp parallel default(shared) firstprivate(c), reduction(+:sum)
      {
         for (auto hpi : pix1)
#pragma omp single nowait
            for (size_t ei = 0; ei < arg1.GetBinning().GetSize(); ++ei) {
               const T y = func(arg1.GetValue(hpi,ei), arg2.GetValue(hpi,ei)) - c;
               const T t = sum + y;
               c = t - sum;
               c -= y;
               sum = t;
            }
      }

      const size_t nEmpty = arg1.Npix() - pix1.size();
      if (nEmpty > 0)
         sum += nEmpty * arg1.GetBinning().GetSize() * func(arg1.GetEmptyValue(), arg2.GetEmptyValue());
   }

   return sum;
}

}

#endif
