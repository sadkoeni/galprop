#ifndef SparseSkymap_h
#define SparseSkymap_h

#include "Skymap.h"
#include <map>
#include <algorithm>

typedef std::pair<long,int> indexType; //!< This is the index for each value in the sparse map, the pixel index and the spectral index

/** \brief Store only set values of pixels.
 *
 * This class was built for counts map for sources, where the psf excludes any emission in a large portion of the map.
 * Stores pixel and spectra indices as well as the pixel value.
 *
 * Only possible to initialize from Skymap object, where 0 values are eliminated.
 */
template <typename T>
class SparseSkymap : public HealpixBaseExtended, public std::map<indexType, T>{
   private:
      SkymapDim fdimension;
      typedef std::pair<const indexType,T> valueType; //!< Typedef for the pixel storage
      typedef typename std::map<indexType, T>::iterator iteratorType; //!< Typedef for the iterator

      /** \brief Helper class to handle copying mapTypes */
      class copyFunctor{
         private:
            SparseSkymap<T> *fmap;
         public:
            copyFunctor( SparseSkymap<T> * map) : fmap(map) {}
            void operator () (const valueType & pixel){
               fmap->insert(pixel);
            }
      };
      /** \brief Helper class to turn mapTypes to Skymaps */
      class fillFunctor{
         private:
            Skymap<T> &fmap; //!<Reference to the skymap to fill
         public:
            fillFunctor(Skymap<T> &map): fmap(map){}
            void operator () (const valueType & pixel){
               fmap[pixel.first.first][pixel.first.second] = pixel.second;
            }
      };


   public:
      /** \brief Default constructor, sets size to 0 */
      SparseSkymap() {Resize(0,0);}
      /** \brief Construct a map of given order and spectral size
       *
       * \param order is the healpix order of the map
       * \param nSpectra is the spectral size of the map
       * \param scheme is the healpix ordering scheme (default value is RING)
       */
      SparseSkymap(int order, int nSpectra, Healpix_Ordering_Scheme scheme=RING){
         Resize(order,nSpectra,scheme);
      }
      /** \brief Construct a map of given order and values at which spectra is evaluated
       *
       * \param order is the healpix order of the map
       * \param spectra is the values at which the spectra is evaluated
       * \param scheme is the healpix ordering scheme (default value is RING)
       */
      SparseSkymap(int order, const std::valarray<double> & spectra, Healpix_Ordering_Scheme scheme=RING){
         Resize(order,spectra,scheme);
      }

      /** \brief Construct a map of given order and values at which spectra is evaluated
       *
       * \param order is the healpix order of the map
       * \param specMin are the lower boundaries for binned maps
       * \param specMax are the lower boundaries for binned maps
       * \param scheme is the healpix ordering scheme (default value is RING)
       */
      SparseSkymap(int order, const std::valarray<double> & specMin, const std::valarray<double> & specMax, Healpix_Ordering_Scheme scheme=RING){
         Resize(order,specMin,specMax,scheme);
      }

      /** \brief Construct a map of given SkymapDim
       *
       * \param dimension describes the dimension of the map
       */
      SparseSkymap(const SkymapDim &dimension){
         Resize(dimension);
      }

      /** \brief Construct a map from a Skymap object.
       *
       * \param Skymap is the skymap to convert to sparse skymap.  All 0 value
       * pixels are removed.
       */
      SparseSkymap(const Skymap<T> & map){
         fromSkymap(map);
      }

      //Copy constructor
      /* Not working at the moment, don't know why, some error about too many
       * arguments to function from for_each
       SparseSkymap(const SparseSkymap & oldMap){
       Resize(oldMap.order_, oldMap.fdimension.spectra, oldMap.scheme_);
      //Copy all the pixels from the old map
      copyFunctor f();
      std::for_each(oldMap.begin(), oldMap.end(), f);
      }*/

      /** \brief Resize the map
       *
       * \param order is the new healpix order
       * \param nSpectra is the new spectral size for the map
       * \param scheme is the new healpix ordering scheme.
       *
       * \note This does not clear the data in the map and it could be rendered invalid with this method
       * if the order is changed.
       */
      void Resize(int order, int nSpectra, Healpix_Ordering_Scheme scheme=RING){
         std::valarray<double> spectra(1.0,nSpectra);
         Resize(order,spectra,scheme);
      }
      /** \brief Resize the map
       *
       * \param dimension is the new dimension of the map
       *
       * \note This does not clear the data in the map and it could be rendered invalid with this method.
       */
      void Resize(const SkymapDim &dimension) {
	Set(dimension.GetOrder(),dimension.GetScheme());
         fdimension = dimension;
      }	
      /** \brief Resize the map
       *
       * \param order is the new healpix order
       * \param spectra is the new values at which the spectra is evaluated
       * \param scheme is the new healpix ordering scheme.
       *
       * \note This does not clear the data in the map and it could be rendered invalid with this method.
       */
      void Resize(int order, const std::valarray<double> & spectra, Healpix_Ordering_Scheme scheme=RING){
         Set(order,scheme);
         fdimension.SetOrder(order);
         fdimension.SetScheme(scheme);
         fdimension.Spectra().resize(spectra.size());
         fdimension.Spectra()=spectra;
         fdimension.SetBinned(false);
      }	
      /** \brief Resize the map
       *
       * \param order is the new healpix order
       * \param specMin are the lower boundaries for binned maps
       * \param specMax are the lower boundaries for binned maps
       * \param scheme is the new healpix ordering scheme.
       *
       * \note This does not clear the data in the map and it could be rendered invalid with this method.
       */
      void Resize(int order, const std::valarray<double> & specMin, const std::valarray<double> & specMax, Healpix_Ordering_Scheme scheme=RING){
         if (specMin.size() != specMax.size()) {
            std::cerr<<"Cannot resize SparseSkymap since minimum and maximum boundaries don't have the same size"<<std::endl;
            throw(1);
         }
         Set(order,scheme);
         fdimension.SetOrder(order);
         fdimension.SetScheme(scheme);
         fdimension.Spectra().resize(specMin.size());
         fdimension.SpectraMin().resize(specMin.size());
         fdimension.SpectraMax().resize(specMin.size());
         fdimension.Spectra()=0.5*(specMin+specMax);
         fdimension.SpectraMin()=specMin;
         fdimension.SpectraMax()=specMax;
         fdimension.SetBinned(true);
      }	

      /** \brief Return the dimension of the skymap */
      const SkymapDim & getDimension() const { return fdimension; }

      /** \brief Return the number of bins in the spectra */
      int nSpectra() const { return fdimension.Spectra().size(); }

      /** \brief Return a constant reference to the values at which the spectra is evaluated */
      const std::valarray<double> & getSpectra() const{
	return fdimension.Spectra();
      }
      /**\brief Return the boundaries for binned maps.
       *
       * \return false if the map is not binned
       */
      bool getBoundaries(std::valarray<double> & specMin, std::valarray<double> & specMax) const{
	if (fdimension.IsBinned()) {
	  specMin.resize(fdimension.SpectraMin().size());
	  specMax.resize(fdimension.SpectraMax().size());
	  specMin = fdimension.SpectraMin();
	  specMax = fdimension.SpectraMax();
         }
	return fdimension.IsBinned();
      }

      /** \brief Convert to a Skymap */
      Skymap<T> toSkymap() const{
         Skymap<T> output(fdimension);
         fillFunctor f(output);
         std::for_each(this->begin(), this->end(), f);
         return output;
      }

      /** \brief Add a healpix map to the sparse map
       *
       * \param map is the input Healpix map.  All non-zero values
       * are added to the current map.
       * \param eBin is the bin number to add the map to
       */
      void addHealpixMap(const Healpix_Map<T> &map, size_t eBin) {
         //Check the order
         if (map.Order() != Order() ) {
            std::cerr<<"Order does not match when adding Healpix map to SparseSkymap.  Nothing added."<<std::endl;
            return;
         }
         //Check the energy bin within bounds
         if (eBin >= fdimension.Spectra().size()){
            std::cerr<<"Energy bin out of bounds when adding Healpix map to SparseSkymap.  Nothing added."<<std::endl;
            return;
         }
         for (long i = 0; i < map.Npix(); ++i){
            if (map[i] != T(0)) {
               indexType ind(i,eBin);
               iteratorType it = this->find(ind);
               if ( it == this->end() ) {
                  it = this->insert(valueType(ind, T(0))).first;
               }
               it->second += map[i];
            }
         }
      }


      /** \brief Convert from a Skymap
       *
       * \param map is the input skymap.  All non-zero values are inserted into the sparse skymap.
       */
      void fromSkymap(const Skymap<T> & map){
         Resize(map.getDimension());
         //Clear the map
         this->clear();
         //Insert all values that are not zero;
         for (long i = 0; i < map.Npix(); ++i){
            for (int j = 0; j < map.nSpectra(); ++j){
               if (map[i][j] != T(0)){
                  indexType index(i,j);
                  (*this)[index] = map[i][j];
               }
            }
         }
      }

      SparseSkymap<T> & operator = (const Skymap<T> &skymap) {
         fromSkymap(skymap);
         return (*this);
      }

      SparseSkymap<T> & operator += (const Skymap<T> &skymap) {
         if ( fdimension != skymap.getDimension() ) {
            std::cerr<<"Skymaps are not equivalent in addition"<<std::endl;
            std::cerr<<"Their sizes are not equal or their spectra is not the same"<<std::endl;
            return (*this);
         }
         //Add all values that are not zero;
         for (long i = 0; i < skymap.Npix(); ++i){
            for (int j = 0; j < skymap.nSpectra(); ++j){
               if (skymap[i][j] != T(0)){
                  //Check for existence of the pixel
                  indexType index(i,j);
                  iteratorType it = this->find(index);

                  if ( it == this->end() ) 
                     (*this)[index] = skymap[i][j];
                  else
                     it->second += skymap[i][j];
               }
            }
         }
         return (*this);
      }

      SparseSkymap<T> & operator -= (const Skymap<T> &skymap) {
         if ( fdimension != skymap.getDimension() ) {
            std::cerr<<"Skymaps are not equivalent in addition"<<std::endl;
            std::cerr<<"Their sizes are not equal or their spectra is not the same"<<std::endl;
            return (*this);
         }
         //Add all values that are not zero;
         for (long i = 0; i < skymap.Npix(); ++i){
            for (int j = 0; j < skymap.nSpectra(); ++j){
               if (skymap[i][j] != T(0)){
                  //Check for existence of the pixel
                  indexType index(i,j);
                  iteratorType it = this->find(index);

                  if ( it == this->end() ) 
                     (*this)[index] = -skymap[i][j];
                  else
                     it->second -= skymap[i][j];
               }
            }
         }
         return (*this);
      }

      SparseSkymap<T> & operator *= (const Skymap<T> &skymap) {
         if ( fdimension != skymap.getDimension() ) {
            std::cerr<<"Skymaps are not equivalent in addition"<<std::endl;
            std::cerr<<"Their sizes are not equal or their spectra is not the same"<<std::endl;
            return (*this);
         }
         //Add all values that are not zero;
         for (long i = 0; i < skymap.Npix(); ++i){
            for (int j = 0; j < skymap.nSpectra(); ++j){
               if (skymap[i][j] != T(0)){
                  //Check for existence of the pixel
                  indexType index(i,j);
                  iteratorType it = this->find(index);

                  if ( it != this->end() ) 
                     it->second *= skymap[i][j];
               }
            }
         }
         return (*this);
      }

      SparseSkymap<T> & operator /= (const Skymap<T> &skymap) {
         if ( fdimension != skymap.getDimension() ) {
            std::cerr<<"Skymaps are not equivalent in addition"<<std::endl;
            std::cerr<<"Their sizes are not equal or their spectra is not the same"<<std::endl;
            return (*this);
         }
         //Add all values that are not zero;
         for (long i = 0; i < skymap.Npix(); ++i){
            for (int j = 0; j < skymap.nSpectra(); ++j){
               if (skymap[i][j] != T(0)){
                  //Check for existence of the pixel
                  indexType index(i,j);
                  iteratorType it = this->find(index);

                  if ( it != this->end() ) 
                     it->second /= skymap[i][j];
               }
            }
         }
         return (*this);
      }

      //Assignment operator
      SparseSkymap<T> & operator = (const SparseSkymap<T> & oldMap){
         //Avoid self-assignment
         if (this != &oldMap){
            Resize(oldMap.fdimension);
            this->clear();
            copyFunctor f(this);
            std::for_each(oldMap.begin(), oldMap.end(), f);
         }
         return(*this);
      }
};

#endif
