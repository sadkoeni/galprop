/**\class Skymap
 * \brief Skymaps for gamma ray calculation in galprop
 *
 * The skymap class is used to store skymaps calculated in gardian
 * The information is stored in a healpix format.
 * Internally we use a std::valarray to store the spectra for each pixel
 * \author Gudlaugur Johannesson
 */

#ifndef SKYMAP_H
#define SKYMAP_H

#include "HealpixBaseExtended.h"
#include "Coordinate.h"
#include "PhysicalConstants.h"
#include "ArraySlice.h"

#include <healpix_map.h>
#include <arr.h>

#include <CCfits/CCfits>

#include <ErrorLogger.h>

#include <fitsio.h>

#include <iostream>
#include <vector>
#include <valarray>
#include <map>
#include <iterator>
#include <algorithm>
#include <functional>
#include <cmath>
#include <typeinfo>
#include <limits>

#include <md5.h>

//Need WCS to read any type of skymaps
#ifdef HAVE_WCS
extern "C" {
#include <wcslib/wcs.h>
#include <wcslib/wcshdr.h>
}
#endif

void empty(double mj, double x, double y, double *z, double *w);

class SkymapDim {
   public:
  //      bool binned;
  //  std::valarray<double> spectra, spectraMin, spectraMax;
  //  int healpixOrder;
  //  Healpix_Ordering_Scheme healpixScheme;

  // Used in the compare function
  enum compareBits {
     ORDERBIT = 1,
     SCHEMEBIT = 2,
     SPSIZEBIT = 4,
     SPBINNEDBIT = 8,
     SPVALUEBIT = 16
  };



  SkymapDim () : binned(false), healpixOrder(-1), healpixScheme(RING) {}
  SkymapDim (const SkymapDim &other);

  ~SkymapDim() {}
  
  SkymapDim & operator = (const SkymapDim & other);
  bool operator == (const SkymapDim & other) const;
  bool operator != (const SkymapDim & other) const;

  // Returns a bitmask that contains information about what is comparable in the two
  // NOTE: SPVALUEBIT is not set unless SPSIZEBIT and SPBINNEDBIT are also set.
  unsigned int compare (const SkymapDim & other) const;
  
  void SetOrder(const int order) { healpixOrder = order; }
  void SetScheme(const Healpix_Ordering_Scheme scheme) { healpixScheme = scheme; }
  void SetBinned(const bool state) { binned = state; }
  
  int GetOrder() const { return healpixOrder; }
  Healpix_Ordering_Scheme GetScheme() const { return healpixScheme; }
  bool IsBinned() const { return binned; }

  const std::valarray<double>& Spectra() const { return spectra; }
  std::valarray<double>& Spectra() { return spectra; }
  
  const std::valarray<double>& SpectraMin() const { return spectraMin; }
  std::valarray<double>& SpectraMin() { return spectraMin; }
  
  const std::valarray<double>& SpectraMax() const { return spectraMax; }
  std::valarray<double>& SpectraMax() { return spectraMax; }

 private:

  bool binned;
  std::valarray<double> spectra, spectraMin, spectraMax;
  int healpixOrder;
  Healpix_Ordering_Scheme healpixScheme;
  
};


/**\brief Store a full sky map in a healpix pixelation.
 *
 * This class is templated to allow for a flexible storage.  Inherits from Healpix_Base so all its methods
 * are available.  Can only be initialized with a whole number order, i.e. nside = 2^order.
 */
template <typename T>
class Skymap : public HealpixBaseExtended {
   private:
      std::valarray<T> fMap; //!< A valarray to store the values, we use slices to get single spectra
      SkymapDim fdimension;

      //!Mask for the skymap basic operators to allow adding maps that don't have the same spectra
      static unsigned int mathMask;

      //!Store the keywords of the FITS file loaded as strings
      std::map<std::string, std::string> fKeywords;

      //Resizes the storage if needed according to fdimension
      void Resize(const T& defaultValue) {
         //Only resize if needed
	const size_t newS = Npix()*fdimension.Spectra().size();
         if (newS == fMap.size()) {
            fMap = defaultValue;
         } else {
	   fMap.resize(Npix()*fdimension.Spectra().size(), defaultValue);
         }

	 //std::cout << "Resize: " << Npix() << " " << fdimension.Spectra().size() << " " << fMap.size() << std::endl;
      }

      /**\brief Convert a CAR map to healpix format
       *
       * Given arrays for CRVAL, CDELT and CRPIX (at the center of the pixel),
       * turn the mapcube into pixels.  The skymap has to be resized properly
       * before entering this method.  The \a image is assumed to be three
       * dimensional in the standard FITS way, longitude evolving fastest, then
       * latitude and last spectra.
       *
       * \param image is the valarray of image data, as from ccfits
       * \param crval is the CRVAL array from the fits header
       * \param crpix is the CRPIX array from the fits header
       * \param cdelt is the CDELT array from the fits header
       */
      void fillMapcube(const std::valarray<T> &image, const std::valarray<long> &axis, const std::valarray<double> &crval, const std::valarray<double> &cdelt, const std::valarray<double> &crpix, bool correctSA = false){
         //Divide each HEALPix pixel into 256 smaller ones and fill those with the
         //value from the map pixel directly beneeth its center.  The HEALPix
         //pixel value is the average of those 256 smaller ones.
         const int dOrder = std::min(5,13-order_); //Order can not exceed 13
         const int nPix = (1<<dOrder)*(1<<dOrder);
         //Create a nested healpix base for us to work with
         const Healpix_Base finer(order_+dOrder, NEST);
         //We must multiply with the solid angle if requested, since the method
         //does not work unless the value of the pixel is independent of its
         //solid angle
         double hpSolidAngle = 1;
         if (correctSA) {
            hpSolidAngle = solidAngle();
         }
#pragma omp parallel for default(shared) schedule(static)
         for (int p = 0; p<npix_; ++p){
            //Use pixel numbering from the nested scheme
            const int pp = (scheme_ == NEST) ? p : ring2nest(p);
            //To store the total value for all the pixels
            std::valarray<T> totSpectra(T(0), fdimension.Spectra().size());
            //Loop over all of the subpixels
            for (int sp = pp*nPix; sp<(pp+1)*nPix; ++sp){
               const pointing pnt = finer.pix2ang(sp);

               //Find the correct index in the array, first the longitude which can
               //loop
               double dil = (pnt.phi*180/utl::kPi-crval[0])/cdelt[0] + crpix[0] + 0.5;
               //If the pixel value is negative, we must loop it over the 360
               //degrees
               if (dil < 1) { 
                  dil += fabs(360./cdelt[0]); 
                  //Pixel values above the boundaries could also indicate looping
               } else if (dil > axis[0]) {
                  dil -= fabs(360./cdelt[0]);
               }
               int il = int(dil - 1); //Since the pixels are 1 based and the array is 0 based

               //Then the latitude, which can't loop
               int ib = int( ( (90 - pnt.theta*180/utl::kPi)-crval[1])/cdelt[1] + crpix[1] + 0.5) - 1;


               //They must both be within bounds to do something useful
               if (ib >= 0 && ib < axis[1] && il >= 0 && il < axis[0]) {
                  const int ind1 = il + ib*axis[0];
                  //We must divide with the solid angle if requested, since the
                  //method does not work unless the value of the pixel is independent
                  //of the solid angle of the pixel
                  double lbSolidAngle = 1;
                  if (correctSA) {
                     double bMiddle = crval[1] + (ib+1 - crpix[1])*cdelt[1];
                     //Note that sin(x-pi/2) = -cos(x)
                     lbSolidAngle = cdelt[0]*utl::kPi/180.*(cos((bMiddle-cdelt[1]/2.)*utl::kPi/180.) - cos((bMiddle+cdelt[1]/2.)*utl::kPi/180.));
                  }
                  for (size_t is=0; is<fdimension.Spectra().size(); ++is){
                     totSpectra[is] += T(image[ind1+is*axis[0]*axis[1]]/lbSolidAngle);
                  }
               } else {
                  //std::cerr<<"Pixels fall outside of boundaries in l and b conversion to healpix"<<std::endl;
                  //std::cerr<<"ib: "<<ib<<", il: "<<il<<std::endl;
                  //std::cerr<<"Phi: "<<pnt.phi*180/utl::kPi<<", theta: "<<90 - pnt.theta*180/utl::kPi<<std::endl;
               }
            }
            size_t first = p*fdimension.Spectra().size();
            for (int l = 0; l < nSpectra(); ++l){
               fMap[first+l] = T(totSpectra[l]/double(nPix)*hpSolidAngle);
            }
         }
      }

      /**\brief Helper class to create a conversion function for coordinate
       * transformation
       */
      class TransformFunction {
         private:
            void (*tfunc)(double,double,double,double*,double*);
         public:
            TransformFunction(SM::CoordSys from, SM::CoordSys to){
               if (from == to)
                  tfunc = &empty;
               else if (from == SM::CoordSys::GAL) {
                  if (to == SM::CoordSys::EQ)
                     tfunc = &SM::gal_eq2;
                  else
                     tfunc = &SM::gal_ecl2;
               } else if (from == SM::CoordSys::EQ) {
                  if (to == SM::CoordSys::ECL)
                     tfunc = &SM::eq_ecl2;
                  else
                     tfunc = &SM::eq_gal2;
               } else {
                  if (to == SM::CoordSys::GAL)
                     tfunc = &SM::ecl_gal2;
                  else
                     tfunc = &SM::ecl_eq2;
               }
            }
            void operator () (double mj, const pointing &from, pointing &to){
               const double pi_=3.141592653589793238462643383279502884197;
               const double fl = from.phi;
               const double fb = pi_/2. - from.theta;
               double tl, tb;
               (*tfunc)(mj,fb,fl,&tb,&tl);
               to.phi = tl;
               to.theta = pi_/2. - tb;
               //Make sure theta is within range
               if (to.theta < 0)
                  to.theta = 0;
               if (to.theta > pi_)
                  to.theta = pi_;
            }
      };


   public:
      /**\brief Default constructor
       *
       * Sets the size to 0 and scheme to RING
       */
      Skymap(){ Resize(0,0); }

      /**\brief Constructor that takes an order of the skymap, the value at which the spectra is evaluated
       * an ordering scheme and a default value for the map.
       *
       * \param order is the order of the healpix coordinate system.  The number
       * of sides per base pixel,  nside, is 2^order.
       * \param spectra is the values at which the spectra is evaluated in each pixel
       * \param scheme is the ordering scheme of the healpix array
       * \param default_value is the default value for the skymap
       */
      Skymap(const int order, const std::valarray<double> spectra, const Healpix_Ordering_Scheme scheme=RING, const T & default_value = T(0)){
         Resize(order,spectra,scheme,default_value);
      }

      /**\brief Constructor that takes an order of the skymap, the
       * boundaries of the spectral bins, an ordering scheme and a default value for the map.
       *
       * \param order is the order of the healpix coordinate system.  The number
       * of sides per base pixel,  nside, is 2^order.
       * \param specMin are the lower boundaries of the spectral bins
       * \param specMax are the upper boundaries of the spectral bins
       * \param scheme is the ordering scheme of the healpix array
       * \param default_value is the default value for the skymap
       */
      Skymap(const int order, const std::valarray<double> specMin, const std::valarray<double> specMax, const Healpix_Ordering_Scheme scheme=RING, const T & default_value = T(0)){
         Resize(order,specMin,specMax,scheme,default_value);
      }

      /**\brief Constructor that takes an order of the skymap, the size of the spectra,
       * an ordering scheme and a default value for the map.
       *
       * \param order is the order of the healpix coordinate system.  The number
       * of sides per base pixel,  nside, is 2^order.
       * \param nSpectra is the size of the spectra
       * \param scheme is the ordering scheme of the healpix array (defaults to RING)
       * \param default_value is the default value for the skymap (defaults to T(0))
       *
       * The values at which the spectra is located is set to 1 in all cases.  This is done so we can take the
       * logarithm of the values, but it is expected that the values at which the spectra is evaluated is in most
       * cases logarithmically distributed.
       */
      Skymap(const int order, const int nSpectra, const Healpix_Ordering_Scheme scheme=RING, const T & default_value = T(0)){
         Resize(order,nSpectra,scheme,default_value);
      }

      /**\brief Constructor that takes a SkymapDim struct
       * and a default value for the map.
       *
       * \param dimensions gives the dimension of the map
       * \param default_value is the default value for the skymap (defaults to T(0))
       *
       */
      Skymap(const SkymapDim &dimension, const T& defaultValue = T(0)) {
         Resize(dimension, defaultValue);
      }
      /**\brief Construct a skymap from file, either a healpix fits file or CAR
       * projected fits image.
       *
       * \param fileName is the name of the file to be opened
       */
      Skymap(const std::string & fileName, int order=-1, bool CorrectSA=false) {
         load(fileName, order, CorrectSA);
      }

      /**\brief Copy constructor */
      Skymap(const Skymap<T> & oldMap){
         Resize(oldMap.fdimension);
         fMap = oldMap.fMap;
      }

      /**\brief Resize the map to order and size of spectra
       *
       * \param order is the order of the healpix coordinate system
       * \param nSpectra is the size of the spectral array
       * \param scheme is the ordering scheme of the healpix map (defaults to RING)
       * \param defaultValue is the default value for the map (defaults to T(0))
       *
       * This method destroys the data in the skymap.  The value at which the spectra is evaluated is set to 1.
       */
      void Resize(const int order, const int nSpectra, const Healpix_Ordering_Scheme scheme=RING, const T & defaultValue = T(0)){
         std::valarray<double> spectra(1.0,nSpectra);
         Resize(order, spectra, scheme, defaultValue);
      }

      /**\brief Resize the map to order and size of spectra
       *
       * \param order is the order of the healpix coordinate system
       * \param spectra is the values at which the spectra is evaluated
       * \param scheme is the ordering scheme of the healpix map (defaults to RING)
       * \param defaultValue is the default value for the map (defaults to T(0))
       *
       * This method destroys the data in the skymap.
       */
      void Resize(const int order, const std::valarray<double> &spectra, const Healpix_Ordering_Scheme scheme=RING, const T & defaultValue = T(0)){
         Set(order, scheme); // Call the construcor for healpix with correct arguments
         fdimension.SetOrder(order);
         fdimension.SetScheme(scheme);
         fdimension.Spectra().resize(spectra.size());
         fdimension.Spectra() = spectra;
         fdimension.SetBinned(false);

         Resize(defaultValue);
      
	 //std::cout << "Skymap: " << order << " " << spectra.size() << " " << fdimension.Spectra().size() << std::endl;
      }

      /**\brief Resize the map to order and size of spectra in a binned fashion
       *
       * \param order is the order of the healpix coordinate system
       * \param specMin are the lower boundaries of the spectral bins
       * \param specMax are the upper boundaries of the spectral bins
       * \param scheme is the ordering scheme of the healpix map (defaults to RING)
       * \param defaultValue is the default value for the map (defaults to T(0))
       *
       * This method destroys the data in the skymap.
       */
      void Resize(const int order, const std::valarray<double> &specMin, const std::valarray<double> &specMax, const Healpix_Ordering_Scheme scheme=RING, const T & defaultValue = T(0)){
         //The spectral sizes must be the same
         if (specMin.size() != specMax.size()){
            std::cerr<<"Spectral sizes not equal for boundary arrays"<<std::endl;
            std::cerr<<specMin.size()<<" != "<<specMax.size()<<std::endl;
            throw(1);
         }
         Set(order, scheme); // Call the construcor for healpix with correct arguments
         fdimension.SetOrder(order);
         fdimension.SetScheme(scheme);
         fdimension.Spectra().resize(specMin.size());
         fdimension.SpectraMax().resize(specMin.size());
         fdimension.SpectraMin().resize(specMin.size());
         fdimension.Spectra() = 0.5*(specMin+specMax);
         fdimension.SpectraMin() = specMin;
         fdimension.SpectraMax() = specMax;
         fdimension.SetBinned(true);

         Resize(defaultValue);
      }

      /**\brief Resize the map from a SkymapDim class 
       *
       * \param dimension is the skymap dimension
       * \param defaultValue is the default value for the map (defaults to T(0))
       *
       * This method destroys the data in the skymap.
       */
      void Resize(const SkymapDim &dimension, const T& defaultValue = T(0)) {
	Set(dimension.GetOrder(),dimension.GetScheme());
         fdimension = dimension;
         Resize(defaultValue);
      }

      /**\brief Load a skymap from a file
       *
       * \param fileName is the name of the file to load.
       * \param order is the order of the skymap if we are rebinning
       * CAR maps.  If -1, it is automatically determined from the
       * CAR binning.  On return it is set to the order of the skymap
       *
       * Tries to be smart and first checks if it is a healpix file, and then
       * tries to load a standard FITS image.  Uses WCS library to read
       * FITS images if available, otherwise assumes CAR projection.  Converts
       * coordinates to GAL if astro is installed, otherwise assumes the
       * input is GAL.
       *
       * CorrectSA starts by dividing with the SA of the pixels before correction
       * and then converts back in the final HEALPix map.  Only relevant when reading
       * FITS images, ignored for healpix fits files.
       */
      void load(const std::string &fileName, int order=-1, bool correctSA=false){
         //Empty the keywords object
         fKeywords.clear();

         //Create the ccfits object, read only, finding the
         //extension with the pixtype HEALPIX and emin keyword
         try {
            std::vector< std::string > keywords(1,"");
            std::vector< std::string > values(1,"");
            keywords[0] = "PIXTYPE";
            values[0] = "HEALPIX";
            CCfits::FITS fits(fileName, CCfits::Read, keywords, values);

            //A reference to the table containing the skymap data
            CCfits::ExtHDU &skymapTable = fits.currentExtension(); 

            //Read the keywords to set up the skymap object
            int nRows, nSpectra, nSide;
            std::string ordering;
            skymapTable.readKey("NAXIS2", nRows);
            skymapTable.readKey("NBRBINS", nSpectra);
            skymapTable.readKey("NSIDE", nSide);
            skymapTable.readKey("ORDERING", ordering);
            //Calculate the order
            int hporder = int(log(double(nSide))/log(2.0)+0.1);// original: int order = int(log(nSide)/log(2)+0.1); AWS20080519

            //Try to find the EBOUNDS or ENERGIES extensions, either must exist
            bool foundEnExt = false;
            //First try the energies table
            try {
               fits.read(std::vector<std::string> (1,"ENERGIES"));
               CCfits::ExtHDU &energyTable = fits.extension("ENERGIES");
               std::valarray<double> energies;
               energyTable.column(1).read(energies, 1, nSpectra);

               if (ordering == "NEST" || ordering == "NESTED") {
		 if (Order() != hporder || nSpectra != int(fdimension.Spectra().size()) || Scheme() != NEST ) {
                     Resize(hporder, energies, NEST);
                  } else {
                     setSpectra(energies);
                  }
               } else {
		 if (Order() != hporder || nSpectra != int(fdimension.Spectra().size()) || Scheme() != RING ) {
		   Resize(hporder, energies, RING);
		 } else {
                     setSpectra(energies);
		 }
               }
               foundEnExt=true;
            } catch (CCfits::FITS::NoSuchHDU) {
               try { //Then EBOUNDS table
                  fits.read(std::vector<std::string> (1,"EBOUNDS"));
                  CCfits::ExtHDU &energyTable = fits.extension("EBOUNDS");
                  std::valarray<double> eMin, eMax;
                  energyTable.column(2).read(eMin, 1, nSpectra);
                  energyTable.column(3).read(eMax, 1, nSpectra);

                  if (ordering == "NEST"|| ordering == "NESTED") {
		    if (Order() != hporder || nSpectra != int(fdimension.Spectra().size()) || Scheme() != NEST ) {
                        Resize(hporder, eMin, eMax, NEST);
                     } else {
                        setSpectra(eMin, eMax);
                     }
                  } else {
		    if (Order() != hporder || nSpectra != int(fdimension.Spectra().size()) || Scheme() != RING ) {
                        Resize(hporder, eMin, eMax, RING);
                     } else {
                        setSpectra(eMin, eMax);
                     }
                  }
                  foundEnExt=true;
               } catch (CCfits::FITS::NoSuchHDU) { }
            }

            //If we found an energy extension, read in the data
            if (foundEnExt) {
               //Read in the data from the table
               //Skip the unnecessary overhead of using CCfits this
               //time and use cfitsio routines
               //Create a map of typeid's to cfitsio datatype
               std::map<const char*, int> formatMap;
               formatMap[typeid(char).name()] = TSBYTE;
               formatMap[typeid(short).name()] = TSHORT;
               formatMap[typeid(int).name()] = TINT;
               formatMap[typeid(long).name()] = TLONG;
               formatMap[typeid(float).name()] = TFLOAT;
               formatMap[typeid(double).name()] = TDOUBLE;
               formatMap[typeid(unsigned char).name()] = TBYTE;
               formatMap[typeid(unsigned short).name()] = TUSHORT;
               formatMap[typeid(unsigned int).name()] = TUINT;
               formatMap[typeid(unsigned long).name()] = TULONG;
               //Select the appropriate datatype
               int dataType = formatMap[typeid(T).name()];

               //Check for old skymap format with column vectors
               int ncols = skymapTable.numCols();
               if (ncols == nSpectra) {
                  //New format
                  long nrows(10);
                  int status(0);
                  fits_get_rowsize(fits.fitsPointer(),&nrows,&status);
                  //default to old method if this doesn't work for some reason
                  if (nrows < 1)
                     nrows = Npix();

                  for (size_t ipix(0); ipix < Npix(); ipix += nrows) {
                     const size_t start = ipix;
                     const size_t stop = std::min(ipix+nrows,size_t(Npix()));
                     for (int icol(1); icol <= nSpectra; ++icol) {
                        std::valarray<T> binMap(stop-ipix);
                        skymapTable.column(icol).read(binMap,ipix+1,stop);
                        for (size_t j(0); j < binMap.size(); ++j)
                           fMap[(j+ipix)*nSpectra + (icol-1)] = binMap[j];
                     }
                  }
               } else {
                  //Old format
                  //Point the fits file pointer to the correct extension
                  skymapTable.makeThisCurrent();
                  //Get the fits pointer
                  fitsfile* ff = fits.fitsPointer();
                  //Read the data
                  int status(0), anynul(0);
                  T null(0);
                  fits_read_col(ff, dataType, 1, 1, 1, fMap.size(), &null, &fMap[0], &anynul, &status);
                  if (status != 0) {
                     fits_report_error(stderr, status);
                     throw(1);
                  }
               }
            } else {
               std::cerr<<"Not a compatible fits file, did not find an EBOUNDS or ENERGIES extension"<<std::endl;
               throw(std::string("Not a compatible fits file, did not find an EBOUNDS or ENERGIES extension"));
            }

            //Read in all the keywords
            skymapTable.readAllKeys();
            std::map<std::string, CCfits::Keyword*>::iterator kit;
            for (kit = skymapTable.keyWord().begin(); kit != skymapTable.keyWord().end(); ++kit) {
               //Find the value type
               CCfits::ValueType vtype = kit->second->keytype();
               std::string value;
               std::ostringstream os;
               switch(vtype) {
                  case CCfits::Tstring:
                     kit->second->value(value);
                     break;

                  case CCfits::Tlogical:
                     bool btmp;
                     kit->second->value(btmp);
                     os << btmp;
                     value = os.str();
                     break;

                  case CCfits::Tbyte:
                     char ctmp;
                     kit->second->value(ctmp);
                     os << ctmp;
                     value = os.str();
                     break;

                  case CCfits::Tshort:
                     short stmp;
                     kit->second->value(stmp);
                     os << stmp;
                     value = os.str();
                     break;

                  case CCfits::Tushort:
                     unsigned short ustmp;
                     kit->second->value(ustmp);
                     os << ustmp;
                     value = os.str();
                     break;

                  case CCfits::Tint:
                     int itmp;
                     kit->second->value(itmp);
                     os << itmp;
                     value = os.str();
                     break;

                  case CCfits::Tuint:
                     unsigned int uitmp;
                     kit->second->value(uitmp);
                     os << uitmp;
                     value = os.str();
                     break;

                  case CCfits::Tlong:
                     long ltmp;
                     kit->second->value(ltmp);
                     os << ltmp;
                     value = os.str();
                     break;

                  case CCfits::Tulong:
                     unsigned long ultmp;
                     kit->second->value(ultmp);
                     os << ultmp;
                     value = os.str();
                     break;

                  case CCfits::Tlonglong:
                     long long lltmp;
                     kit->second->value(lltmp);
                     os << lltmp;
                     value = os.str();
                     break;

                  case CCfits::Tfloat:
                     float ftmp;
                     kit->second->value(ftmp);
                     os << ftmp;
                     value = os.str();
                     break;
                  case CCfits::Tdouble:
                     double dtmp;
                     kit->second->value(dtmp);
                     os << dtmp;
                     value = os.str();
                     break;

                  default:
                     value = "";
               }

               fKeywords[kit->first] = value;
            }

         } catch (CCfits::FITS::NoSuchHDU) {

            //Read from the primary image
            CCfits::FITS fits(fileName);
            CCfits::PHDU &mapCube = fits.pHDU();

#ifdef HAVE_WCS
            std::ostringstream oss;

            INFO("Using WCS to read the fits file");
            //Read the file using the wcslib
            char *header;
            int nkeys, status(0);

            mapCube.makeThisCurrent();

            fits_convert_hdr2str(mapCube.fitsPointer(), 1, NULL, 0, &header, &nkeys, &status);
            if (status != 0) {
               ERROR("Could not read header of file");
               throw(std::runtime_error("Broken file"));
            }

            //Read the WCS keywords
            wcsprm *wcsData;
            int ctrl=0, nreject, nwcs;
            status = wcspih(header, nkeys, 0, ctrl, &nreject, &nwcs, &wcsData);
            if (status != 0) {
               ERROR("Could not read WCS information");
               throw(std::runtime_error("Broken file"));
            }

            fits_free_memory(header, &status);
            if (status != 0) {
               ERROR("Could not free memory");
               throw(std::runtime_error("Major error"));
            }

            //Assert that we have at least 1 projection and the first has 2 axis, nothing complicated here
            if (nwcs < 1 || wcsData[0].naxis < 2) {
               ERROR("Number of axis less than 2, cannot continue");
               throw(std::runtime_error("Number of axis less than 2, cannot continue"));
            }

            //Make sure the files do not have too many axis
            if (wcsData[0].naxis > 3) {
               ERROR("More than 3 axis in the image, cannot continue");
               throw(std::runtime_error("Can at most handle a single spectrum"));
            }

            wcsprm *wcs = &wcsData[0];

            //Find the longitude and latitude range of the input map
            //Loop over all pixels and find the range
            double lmin = std::numeric_limits<double>::max();
            double lmax = -std::numeric_limits<double>::max();
            double bmin = std::numeric_limits<double>::max();
            double bmax = -std::numeric_limits<double>::max();

            //In case there are three axes, we don't know before hand which ones are the
            //sky coordinates.  
            //Pixel coordinates are 1 based
            if (mapCube.axes() == 2) {

               //Here it is clear and easy
               size_t mapPix = mapCube.axis(1);
               double worldTmp[mapPix][3];
               double phiTmp[mapPix], thetaTmp[mapPix];
               double imgcrdTmp[mapPix][3], pixcrdTmp[mapPix][3];
               int statTmp[mapPix];

               for (size_t i0 = 0; i0 < mapCube.axis(0); ++i0) {
                  size_t count(0);
                  for (size_t i1 = 0; i1 < mapCube.axis(1); ++i1) {
                     pixcrdTmp[count][0] = i0+1;
                     pixcrdTmp[count][1] = i1+1;
                     ++count;
                  }

                  wcsp2s(wcs, count, 3, pixcrdTmp[0], imgcrdTmp[0], phiTmp, thetaTmp, worldTmp[0], statTmp);

                  for (size_t i = 0; i < count; ++i) {
                     lmin = std::min(lmin, worldTmp[i][wcs->lng] - fabs(wcs->cdelt[wcs->lng]));
                     lmax = std::max(lmax, worldTmp[i][wcs->lng] + fabs(wcs->cdelt[wcs->lng]));
                     bmin = std::min(bmin, worldTmp[i][wcs->lat] - fabs(wcs->cdelt[wcs->lat]));
                     bmax = std::max(bmax, worldTmp[i][wcs->lat] + fabs(wcs->cdelt[wcs->lat]));
                  }
               }


            } else {

               //Here it is not so clear.  Could do a test run over the axis to see which one is spectrum
               //For now just loop over them all
               size_t mapPix = mapCube.axis(1);
               double worldTmp[mapPix][3];
               double phiTmp[mapPix], thetaTmp[mapPix];
               double imgcrdTmp[mapPix][3], pixcrdTmp[mapPix][3];
               int statTmp[mapPix];

               for (size_t i2 = 0; i2 < mapCube.axis(2); ++i2) {

                  for (size_t i0 = 0; i0 < mapCube.axis(0); ++i0) {
                     size_t count(0);
                     for (size_t i1 = 0; i1 < mapCube.axis(1); ++i1) {
                        pixcrdTmp[count][0] = i0+1;
                        pixcrdTmp[count][1] = i1+1;
                        pixcrdTmp[count][2] = i2+1;
                        ++count;
                     }

                     wcsp2s(wcs, count, 3, pixcrdTmp[0], imgcrdTmp[0], phiTmp, thetaTmp, worldTmp[0], statTmp);

                     for (size_t i = 0; i < count; ++i) {
                        lmin = std::min(lmin, worldTmp[i][wcs->lng] - fabs(wcs->cdelt[wcs->lng]));
                        lmax = std::max(lmax, worldTmp[i][wcs->lng] + fabs(wcs->cdelt[wcs->lng]));
                        bmin = std::min(bmin, worldTmp[i][wcs->lat] - fabs(wcs->cdelt[wcs->lat]));
                        bmax = std::max(bmax, worldTmp[i][wcs->lat] + fabs(wcs->cdelt[wcs->lat]));
                     }

                  }
               }
            }

            if (mapCube.axes() != wcs->naxis) {
               ERROR("Dimensions for wcs and FITS do not match");
               throw(std::runtime_error("This just shouldn't happen"));
            }

            //wcs->lat and wcs->lng are not set until after the first call to wcsp2s or wcss2p
            if (wcs->lat < 0 || wcs->lng < 0) {
               ERROR("Need two celestial axis");
               throw(std::runtime_error("Cannot convert this map"));
            }

            oss.str("");
            oss<<"Axis for longitude: "<<wcs->lng+1;
            DEBUGLOG(oss.str());
            oss.str("");
            oss<<"Axis for latitude: "<<wcs->lat+1;
            DEBUGLOG(oss.str());

            int spec = wcs->spec;
            if (spec < 0) {
               size_t i = 0;
               while (spec < 0 && i < mapCube.axes()) {
                  if (i != wcs->lat && i != wcs->lng) 
                     spec = i;
                  ++i;
               }
            }

            //We project everything to GAL and can only accept Equatorial and Ecliptic
            void (*tfunc)(double,double,double,double*,double*);
            std::string Coord = wcs->lngtyp;
            //Try to use ctype if Coord is empty, handles old gasmaps for galprop
            if (Coord == "    ")
               Coord = wcs->ctype[wcs->lng];

            oss.str("");
            oss<<"Coordinates of input map is";
            if (Coord == "RA") {
               tfunc = &SM::gal_eq2;
               oss<<" Equatorial";
            } else if (Coord == "ELON") {
               tfunc = &SM::gal_ecl2;
               oss<<" Ecliptic";
            } else if (Coord == "GLON") {
               tfunc = &empty;
               oss<<" Galactic";
            } else {
               ERROR("Unknown coordinate system \""+Coord+"\", assuming Galactic");
               oss<<" Galactic";
            }
            DEBUGLOG(oss.str());
            oss.str("");

            oss<<"Maximum dimensions of input map (longitude) : (latitude) : ("<<lmin<<", "<<lmax<<") : ("<<bmin<<", "<<bmax<<")";
            DEBUGLOG(oss.str());
            oss.str("");

            oss<<"Projection name: "<<wcs->cel.prj.name;
            DEBUGLOG(oss.str());

            for ( size_t i(0); i < wcs->naxis; ++i) {
               oss.str("");
               oss<<"CRVAL"<<i+1<<": "<<wcs->crval[i];
               DEBUGLOG(oss.str());
               oss.str("");
               oss<<"CRPIX"<<i+1<<": "<<wcs->crpix[i];
               DEBUGLOG(oss.str());
               oss.str("");
               oss<<"CDELT"<<i+1<<": "<<wcs->cdelt[i];
               DEBUGLOG(oss.str());
               oss.str("");
               oss<<"CUNIT"<<i+1<<": "<<wcs->cunit[i];
               DEBUGLOG(oss.str());
               oss.str("");
               oss<<"CTYPE"<<i+1<<": "<<wcs->ctype[i];
               DEBUGLOG(oss.str());
               for (size_t j(0); j < wcs->naxis; ++j) {
                  oss.str("");
                  oss<<"PC"<<i+1<<"_"<<j+1<<": "<<wcs->pc[wcs->naxis*i+j];
                  DEBUGLOG(oss.str());
               }
               oss.str("");
               oss<<"CROTA"<<i+1<<": "<<wcs->crota[i];
               DEBUGLOG(oss.str());
               for (size_t j(0); j < wcs->naxis; ++j) {
                  oss.str("");
                  oss<<"CD"<<i+1<<"_"<<j+1<<": "<<wcs->cd[wcs->naxis*i+j];
                  DEBUGLOG(oss.str());
               }
            }

            oss.str("");
            oss<<"cel.flag: "<<wcs->cel.flag;
            DEBUGLOG(oss.str());
            oss.str("");
            oss<<"cel.offset: "<<wcs->cel.offset;
            DEBUGLOG(oss.str());
            oss.str("");
            oss<<"cel.phi0: "<<wcs->cel.phi0;
            DEBUGLOG(oss.str());
            oss.str("");
            oss<<"cel.theta0: "<<wcs->cel.theta0;
            DEBUGLOG(oss.str());
            oss.str("");
            oss<<"cel.ref[0]: "<<wcs->cel.ref[0];
            DEBUGLOG(oss.str());
            oss.str("");
            oss<<"cel.ref[1]: "<<wcs->cel.ref[1];
            DEBUGLOG(oss.str());
            oss.str("");
            oss<<"cel.ref[2]: "<<wcs->cel.ref[2];
            DEBUGLOG(oss.str());
            oss.str("");
            oss<<"cel.ref[3]: "<<wcs->cel.ref[3];
            DEBUGLOG(oss.str());

            //If this is a proper mapcube, read in the energies value
            std::valarray<double> energies, emin, emax;
            try {
               //Try energies extension first
               CCfits::ExtHDU & energyTable = fits.extension("ENERGIES");
               int nSpectra;
               energyTable.readKey("NAXIS2", nSpectra);
               energyTable.column(1).read(energies, 1, nSpectra);
               DEBUGLOG("Energies from ENERGIES extension");
            } catch (CCfits::FITS::NoSuchHDU) {
               try{
                  //Then ebounds
                  CCfits::ExtHDU & energyTable = fits.extension("EBOUNDS");
                  int nSpectra;
                  energyTable.readKey("NAXIS2", nSpectra);
                  energyTable.column(2).read(emin, 1, nSpectra);
                  energyTable.column(3).read(emax, 1, nSpectra);
                  energies.resize(emin.size());
                  DEBUGLOG("Energies from EBOUNDS extension");
               } catch (CCfits::FITS::NoSuchHDU) {}
            }

            //If no energies extension is found, create the values from CRVAL and CDELT values
            //Assume the GALPROP convention for storing the spectral information axis
            if (energies.size() == 0) {
               if (wcs->naxis < 3) {
                  DEBUGLOG("No energy axis in map");
                  energies.resize(1);
                  energies[0] = 1;
               } else {
                  energies.resize(mapCube.axis(spec));
                  //Assume CDELT represents logarithmic values
                  for (int i = 0; i < energies.size(); ++i) {
                     energies[i] = pow(10,wcs->crval[spec] + i*wcs->cdelt[spec]);
                  }
                  oss.str("");
                  oss<<"Energy axes created from axis "<<spec<<" from input map assuming CDELT is log10";
                  DEBUGLOG(oss.str());
               }
            }

            //The resolution determines the order for the map, which is always
            //greater or equal to the resolution of the map
            double res = std::min(fabs(wcs->cdelt[wcs->lng]),fabs(wcs->cdelt[wcs->lat]));
            if (order == -1) {
               order = int(log(sqrt(3./utl::kPi)*60/res)/log(2.0)) + 1; //log(nside)/log(2) + 1      log(2)->log(2.0) AWS20080519
               oss.str("");
               oss<<"Setting order of map to "<<order;
               DEBUGLOG(oss.str());
            } else {
               oss.str("");
               oss<<"Using user provided order of "<<order;
               DEBUGLOG(oss.str());
            }

            //Now we can set up the map
            if (emin.size() == 0){
	      if (Order() != order || energies.size() != fdimension.Spectra().size() || Scheme() != RING ) {
                  Resize(order, energies, RING);
               }else{
                  setSpectra(energies);
               }
            } else {
	      if (Order() != order || emin.size() != fdimension.Spectra().size() || Scheme() != RING ) {
                  Resize(order, emin, emax, RING);
               }else{
                  setSpectra(emin, emax);
               }
            }

            //Read the image data and fill the skymap
            std::valarray<T> image;
            long npix(1);
            for (int i = 0; i < wcs->naxis; ++i){
               npix *= mapCube.axis(i);
            }
            mapCube.read(image,1,npix);

            //Correct the map for solid angle if necessary
            if (correctSA) {
               //One more than number of pixels
               size_t nvert = (mapCube.axis(wcs->lng)+1)*(mapCube.axis(wcs->lat)+1);
               
               //Arrays for conversion
               double world[nvert][3];
               double phi[nvert];  //Also used for solid angle
               double theta[nvert];
               double imgcrd[nvert][3];
               double pixcrd[nvert][3];
               int stat[nvert];

               size_t ind(0);
               for (int ib = 0; ib < mapCube.axis(wcs->lat)+1; ++ib) {
                  for (int il = 0; il < mapCube.axis(wcs->lng)+1; ++il) {
                     //Reset the values, in case spec is not defined
                     pixcrd[ind][0] = pixcrd[ind][1]= pixcrd[ind][2] = 1;
                     pixcrd[ind][wcs->lng] = il+0.5;
                     pixcrd[ind][wcs->lat] = ib+0.5;
                     ++ind;
                  }
               }
               wcsp2s(wcs, ind, 3, pixcrd[0], imgcrd[0], phi, theta, world[0], stat);

               //For conversion between il, ib, ispec to map index
               size_t ilmul = 1;
               if (wcs->lng == 1)
                  ilmul *= mapCube.axis(0);
               else if (wcs->lng == 2)
                  ilmul *= mapCube.axis(0)*mapCube.axis(1);
               size_t ibmul = 1;
               if (wcs->lat == 1)
                  ibmul *= mapCube.axis(0);
               else if (wcs->lat == 2)
                  ibmul *= mapCube.axis(0)*mapCube.axis(1);
               size_t ispecmul = 1;
               if (spec == 1)
                  ispecmul *= mapCube.axis(0);
               else if (spec == 2)
                  ispecmul *= mapCube.axis(0)*mapCube.axis(1);

               //Loop through the pixels and apply the solid angle correction
               //Need to specialize according to spec
               //Store the solid angle in the phi vector;
               size_t pind(0);
               for (int ib = 0; ib < mapCube.axis(wcs->lat); ++ib) {
                  for (int il = 0; il < mapCube.axis(wcs->lng); ++il) {
                     //The indices to the pixel corners
                     const size_t i1 = ib*(mapCube.axis(wcs->lng)+1)+il;
                     const size_t i2 = i1+1;
                     const size_t i4 = i2+mapCube.axis(wcs->lng);
                     const size_t i3 = i4+1;

                     //Need the distances between the corners, cache sin and cos of the lat
                     const double sinlat1 = std::sin(world[i1][wcs->lat]*utl::kConvertDegreesToRadians);
                     const double coslat1 = std::cos(world[i1][wcs->lat]*utl::kConvertDegreesToRadians);
                     const double sinlat2 = std::sin(world[i2][wcs->lat]*utl::kConvertDegreesToRadians);
                     const double coslat2 = std::cos(world[i2][wcs->lat]*utl::kConvertDegreesToRadians);
                     const double sinlat3 = std::sin(world[i3][wcs->lat]*utl::kConvertDegreesToRadians);
                     const double coslat3 = std::cos(world[i3][wcs->lat]*utl::kConvertDegreesToRadians);
                     const double sinlat4 = std::sin(world[i4][wcs->lat]*utl::kConvertDegreesToRadians);
                     const double coslat4 = std::cos(world[i4][wcs->lat]*utl::kConvertDegreesToRadians);

                     //Distances between adjacent corners
                     const double d12 = SM::chacos(sinlat1*sinlat2+coslat1*coslat2*std::cos((world[i1][wcs->lng]-world[i2][wcs->lng])*utl::kConvertDegreesToRadians));
                     const double d14 = SM::chacos(sinlat1*sinlat4+coslat1*coslat4*std::cos((world[i1][wcs->lng]-world[i4][wcs->lng])*utl::kConvertDegreesToRadians));
                     const double d32 = SM::chacos(sinlat3*sinlat2+coslat3*coslat2*std::cos((world[i3][wcs->lng]-world[i2][wcs->lng])*utl::kConvertDegreesToRadians));
                     const double d34 = SM::chacos(sinlat3*sinlat4+coslat3*coslat4*std::cos((world[i3][wcs->lng]-world[i4][wcs->lng])*utl::kConvertDegreesToRadians));

                     //In case of triangles
                     if (d12 <= 0 || d14 <= 0) {
                        //Diagonal
                        const double d24 = SM::chacos(sinlat2*sinlat4+coslat2*coslat4*std::cos((world[i2][wcs->lng]-world[i4][wcs->lng])*utl::kConvertDegreesToRadians));

                        //Cache for sines and cosines
                        const double sind32 = std::sin(d32);
                        const double sind24 = std::sin(d24);
                        const double sind34 = std::sin(d34);
                        const double cosd32 = std::cos(d32);
                        const double cosd24 = std::cos(d24);
                        const double cosd34 = std::cos(d34);

                        //The angles
                        const double a2 = SM::chacos((cosd34-cosd24*cosd32)/(sind24*sind32));
                        const double a3 = SM::chacos((cosd24-cosd34*cosd32)/(sind34*sind32));
                        const double a4 = SM::chacos((cosd32-cosd24*cosd34)/(sind24*sind34));

                        phi[pind] = (a2+a3+a4) - utl::kPi;

                     } else if ( d32 <= 0 || d34 <= 0) {
                        //Diagonal
                        const double d24 = SM::chacos(sinlat2*sinlat4+coslat2*coslat4*std::cos((world[i2][wcs->lng]-world[i4][wcs->lng])*utl::kConvertDegreesToRadians));

                        //Cache for sines and cosines
                        const double sind12 = std::sin(d12);
                        const double sind24 = std::sin(d24);
                        const double sind14 = std::sin(d14);
                        const double cosd12 = std::cos(d12);
                        const double cosd24 = std::cos(d24);
                        const double cosd14 = std::cos(d14);

                        //The angles
                        const double a1 = SM::chacos((cosd24-cosd12*cosd14)/(sind12*sind14));
                        const double a2 = SM::chacos((cosd14-cosd12*cosd24)/(sind12*sind24));
                        const double a4 = SM::chacos((cosd12-cosd14*cosd24)/(sind14*sind24));

                        phi[pind] = (a1+a2+a4) - utl::kPi;

                     } else {
                        //Complete polygon
                        const double d13 = SM::chacos(sinlat1*sinlat3+coslat1*coslat3*std::cos((world[i1][wcs->lng]-world[i3][wcs->lng])*utl::kConvertDegreesToRadians));
                        const double d24 = SM::chacos(sinlat2*sinlat4+coslat2*coslat4*std::cos((world[i2][wcs->lng]-world[i4][wcs->lng])*utl::kConvertDegreesToRadians));

                        //Cache for sines and cosines
                        const double sind12 = std::sin(d12);
                        const double sind14 = std::sin(d14);
                        const double sind32 = std::sin(d32);
                        const double sind34 = std::sin(d34);
                        const double cosd12 = std::cos(d12);
                        const double cosd13 = std::cos(d13);
                        const double cosd14 = std::cos(d14);
                        const double cosd24 = std::cos(d24);
                        const double cosd32 = std::cos(d32);
                        const double cosd34 = std::cos(d34);

                        //The angles
                        const double a1 = SM::chacos((cosd13-cosd34*cosd14)/(sind34*sind14));
                        const double a2 = SM::chacos((cosd24-cosd32*cosd34)/(sind32*sind34));
                        const double a3 = SM::chacos((cosd13-cosd12*cosd32)/(sind12*sind32));
                        const double a4 = SM::chacos((cosd24-cosd14*cosd12)/(sind14*sind12));

                        phi[pind] = (a1+a2+a3+a4) - utl::kTwoPi;
                        //std::cout<<il<<", "<<ib<<", "<<phi[pind]<<", "<<a1<<", "<<a2<<", "<<a3<<", "<<a4<<std::endl;
                        //std::cout<<d12<<", "<<d13<<", "<<d14<<", "<<d32<<", "<<d24<<", "<<d34<<std::endl;
                        //std::cout<<world[i1][wcs->lng]<<", "<<world[i2][wcs->lng]<<", "<<world[i3][wcs->lng]<<", "<<world[i4][wcs->lng]<<std::endl;
                        //std::cout<<world[i1][wcs->lat]<<", "<<world[i2][wcs->lat]<<", "<<world[i3][wcs->lat]<<", "<<world[i4][wcs->lat]<<std::endl;
                     }

                     if (phi[pind] > 0) {
                        if (spec < 0) 
                           image[il*ilmul+ib*ibmul] /= phi[pind];
                        else {
                           const size_t imind = il*ilmul+ib*ibmul;
                           for (int ispec = 0; ispec < mapCube.axis(spec); ++ispec)
                              image[imind+ispec*ispecmul] /= phi[pind];
                        }
                     }

                     ++pind;
                  }
               }
               

            }

            //Divide each HEALPix pixel into 256 smaller ones and fill those with the
            //value from the map pixel directly beneeth its center.  The HEALPix
            //pixel value is the average of those 256 smaller ones.
            const int dOrder = std::min(5,13-order_); //Order can not exceed 13
            const int nPix = (1<<dOrder)*(1<<dOrder);

            //Create a nested healpix base for us to work with
            const Healpix_Base finer(order_+dOrder, NEST);

#pragma omp parallel for default(shared) schedule(static)
            for (int p = 0; p<npix_; ++p){

               //Use pixel numbering from the nested scheme
               const int pp = ring2nest(p);

               //Set up arrays for all the pixels, use only the first index in the spectral coordinate
               double world[nPix][3];
               double phi[nPix];
               double theta[nPix];
               double imgcrd[nPix][3];
               double pixcrd[nPix][3];
               int stat[nPix];

               //Loop over all of the subpixels
               size_t count(0);
               for (int sp = pp*nPix; sp<(pp+1)*nPix; ++sp){
                  const pointing pnt = finer.pix2ang(sp);

                  //Do coordinate transformation
                  const double fl = pnt.phi;
                  const double fb = utl::kPi/2. - pnt.theta;


		  double tl, tb;
                  //Assume J2000
                  (*tfunc)(J2000,fb,fl,&tb,&tl);

                  //Convert to degrees
                  tb *= 180./utl::kPi;
                  tl *= 180./utl::kPi;
                  if (tb > bmin && tb < bmax) {

                     if ( tl > lmax )
                        tl -= 360;

                     if ( tl < lmin )
                        tl += 360;

                     if ( tl > lmin && tl < lmax ) {
                        //Reset the values, in case spec is not defined
                        world[count][0] = world[count][1]= world[count][2] = 1;
                        world[count][wcs->lng] = tl;
                        world[count][wcs->lat] = tb;
                        ++count;
                     } else {
                        //std::cout<<tl<<", "<<lmin<<", "<<lmax<<std::endl;
                     }
                  } else {
                     //std::cout<<tb<<", "<<bmin<<", "<<bmax<<std::endl;
                  }

               }

               //No valid pixels
               if (count == 0)
                  continue;

               //Do the conversion
               wcss2p(wcs, count, 3, world[0], phi, theta, imgcrd[0], pixcrd[0], stat);
               if (status != 0) {
                  ERROR("wcss2p failed, cannot convert map to healpix");
                  throw(std::runtime_error("WCS failure"));
               }

               //To store the total value for all the pixels
               std::valarray<T> totSpectra(T(0), fdimension.Spectra().size());

               //Now we must loop over the pixels and take the average
               //What happens next depends on spec
               size_t numValid(0);
               if ( spec < 0 ) {
                  for (int i(0); i < count; ++i) {
                     //No spectral information, only a single value
                     int i0 = int(round(pixcrd[i][0]-1));
                     int i1 = int(round(pixcrd[i][1]-1));

                     //They must both be within bounds to do something useful
                     if (i0 >= 0 && i0 < mapCube.axis(0) && i1 >= 0 && i1 < mapCube.axis(1)) {
                        ++numValid;
                        const int ind1 = i0 + i1*mapCube.axis(0);
                        //There is only one value in the spectra
                        totSpectra[0] += T(image[ind1]);
                     } else {
                        //std::cout<<"No spec: "<<i0<<", "<<i1<<", "<<mapCube.axis(0)<<", "<<mapCube.axis(1)<<std::endl;
                     }
                  }
               } else if (spec == 0) {
                  for (int i(0); i < count; ++i) {
                     //Spectral information in first axis
                     int i1 = int(round(pixcrd[i][1]-1));
                     int i2 = int(round(pixcrd[i][2]-1));

                     //They must both be within bounds to do something useful
                     if (i1 >= 0 && i1 < mapCube.axis(1) && i2 >= 0 && i2 < mapCube.axis(2)) {
                        ++numValid;
                        const int ind1 = i1*mapCube.axis(0) + i2*mapCube.axis(0)*mapCube.axis(1);

#if (_OPENMP >= 201307)
#pragma omp simd
#endif
                        for (size_t is=0; is<fdimension.Spectra().size(); ++is){
                           totSpectra[is] += T(image[ind1+is]);
                        }
                     } else {
                        //std::cout<<"Spec 0: "<<i1<<", "<<i2<<", "<<mapCube.axis(1)<<", "<<mapCube.axis(2)<<std::endl;
                     }
                  }
               } else if (spec == 1) {
                  for (int i(0); i < count; ++i) {
                     //Spectral information in second axis
                     int i0 = int(round(pixcrd[i][0]-1));
                     int i2 = int(round(pixcrd[i][2]-1));

                     //They must both be within bounds to do something useful
                     if (i0 >= 0 && i0 < mapCube.axis(0) && i2 >= 0 && i2 < mapCube.axis(2)) {
                        ++numValid;
                        const int ind1 = i0 + i2*mapCube.axis(0)*mapCube.axis(1);

#if (_OPENMP >= 201307)
#pragma omp simd
#endif
                        for (size_t is=0; is<fdimension.Spectra().size(); ++is){
                           totSpectra[is] += T(image[ind1+is*mapCube.axis(0)]);
                        }
                     } else {
                        //std::cout<<"Spec 1: "<<i0<<", "<<i2<<", "<<mapCube.axis(0)<<", "<<mapCube.axis(2)<<std::endl;
                     }
                  }
               } else if (spec == 2) {
                  for (int i(0); i < count; ++i) {
                     //Spectral information in third axis
                     int i0 = int(round(pixcrd[i][0]-1));
                     int i1 = int(round(pixcrd[i][1]-1));

                     //They must both be within bounds to do something useful
                     if (i0 >= 0 && i0 < mapCube.axis(0) && i1 >= 0 && i1 < mapCube.axis(1)) {
                        ++numValid;
                        const int ind1 = i0 + i1*mapCube.axis(0);

#if (_OPENMP >= 201307)
#pragma omp simd
#endif
                        for (size_t is=0; is<fdimension.Spectra().size(); ++is){
                           totSpectra[is] += T(image[ind1+is*mapCube.axis(0)*mapCube.axis(1)]);
                        }
                     } else {
                        //std::cout<<"Spec 2: "<<i0<<", "<<i1<<", "<<mapCube.axis(0)<<", "<<mapCube.axis(1)<<std::endl;
                        //std::cout<<pixcrd[i][0]<<" : "<<world[i][0]<<std::endl;
                        //std::cout<<pixcrd[i][1]<<" : "<<world[i][1]<<std::endl;
                     }
                  }
               }

               if (numValid > 0) {
                  size_t first = p*fdimension.Spectra().size();
                  for (int l = 0; l < nSpectra(); ++l){
                     fMap[first+l] = T(totSpectra[l]/double(numValid));
                  }
                  if (correctSA) {
                     const double SA = solidAngle();
                     for (int l = 0; l < nSpectra(); ++l)
                        fMap[first+l] *= SA;
                  }
               }
            }



#else
            INFO("Assuming we have a CAR fits image");
            //Assume we have a CAR fits image
            //Read the number of axes
            long axes = mapCube.axes();

            //Throw an error if the number of axes is less than 2
            if (axes < 2) {
               ERROR("Number of axis less than 2, cannot continue");
               throw(std::runtime_error("Number of axis less than 2, cannot continue"));
            }

            //We take at most 3 axes into account
            axes = std::min(long(3),axes);

            //Create a vector of keywords for the CR values
            std::stringstream ss;  //To write the numbers to
            std::valarray<double> crval(0.0, 3), crpix(1.0, 3), cdelt(1.0, 3);
            std::valarray<long> axis(3);

            //stringstreams should be initialized before str("") is called
            ss << 'a';
            for (int i = 0; i < axes; ++i) {
               axis[i] = mapCube.axis(i);
               //Seek to the beginning of the stringstream to overwrite old values
               ss.str("");
               ss << "CRPIX" << i+1;
               try {
                  mapCube.readKey(ss.str(), crpix[i]);
               } catch (CCfits::HDU::NoSuchKeyword) {} //Assume the value is 1 if undefined
               //Seek to the beginning of the stringstream to overwrite old values
               ss.str("");
               ss << "CDELT" << i+1;
               try {
                  mapCube.readKey(ss.str(), cdelt[i]);
               } catch (CCfits::HDU::NoSuchKeyword) {
                  //Assume whole sky maps and 1 for all
                  //other axis
                  if (i == 0) {
                     cdelt[i] = 360./axis[i];
                  } else if (i == 1) {
                     cdelt[i] = 180./axis[i];
                  }
               } 
               //Seek to the beginning of the stringstream to overwrite old values
               ss.str("");
               ss << "CRVAL" << i+1;
               try {
                  mapCube.readKey(ss.str(), crval[i]);
               } catch (CCfits::HDU::NoSuchKeyword) {
                  //Assume full sky maps and 0 for everything else
                  if (i == 0) {
                     crval[i] = 0 + cdelt[i]/2.;
                  } else if (i == 1) {
                     crval[i] = -90 + cdelt[i]/2.;
                  }
               } 
            }

            //Read the data and resize the skymap.  Let the skymap be of RING
            //structure, most favorable for convolution.
            //The resolution determines the order for the map, which is always
            //greater or equal to the resolution of the map
            double res = std::min(fabs(cdelt[0]),fabs(cdelt[1]));
            if (order == -1)
               order = int(log(sqrt(3./utl::kPi)*60/res)/log(2.0)) + 1; //log(nside)/log(2) + 1      log(2)->log(2.0) AWS20080519

            //If this is a proper mapcube, read in the energies value
            std::valarray<double> energies, emin, emax;
            try {
               //Try energies extension first
               CCfits::ExtHDU & energyTable = fits.extension("ENERGIES");
               int nSpectra;
               energyTable.readKey("NAXIS2", nSpectra);
               energyTable.column(1).read(energies, 1, nSpectra);
            } catch (CCfits::FITS::NoSuchHDU) {
               try{
                  //Then ebounds
                  CCfits::ExtHDU & energyTable = fits.extension("EBOUNDS");
                  int nSpectra;
                  energyTable.readKey("NAXIS2", nSpectra);
                  energyTable.column(2).read(emin, 1, nSpectra);
                  energyTable.column(3).read(emax, 1, nSpectra);
                  energies.resize(emin.size());
               } catch (CCfits::FITS::NoSuchHDU) {}
            }

            //If no energies extension is found, create the values from CRVAL and
            //CDELT values
            if (energies.size() == 0) {
               if (axes < 3) {
                  energies.resize(1);
                  energies[0] = 1;
                  axis[2] = 1;
               }else{
                  energies.resize(axis[2]);
                  //Assume CDELT represents logarithmic values
                  for (int i = 0; i < axis[2]; ++i) {
                     energies[i] = pow(10,crval[2] + i*cdelt[2]);
                  }
               }
            }

            //Now we can set up the map
            if (emin.size() == 0){
	      if (Order() != order || energies.size() != fdimension.Spectra().size() || Scheme() != RING ) {
                  Resize(order, energies, RING);
               }else{
                  setSpectra(energies);
               }
            } else {
	      if (Order() != order || emin.size() != fdimension.Spectra().size() || Scheme() != RING ) {
                  Resize(order, emin, emax, RING);
               }else{
                  setSpectra(emin, emax);
               }
            }

            //Read the image data and fill the skymap
            std::valarray<T> image;
            long npix(1);
            for (int i = 0; i < axes; ++i){
               npix *= axis[i];
            }
            mapCube.read(image,1,npix);
            fillMapcube(image, axis, crval, cdelt, crpix);
#endif

            //Read in all the keywords
            mapCube.readAllKeys();
            std::map<std::string, CCfits::Keyword*>::iterator kit;
            for (kit = mapCube.keyWord().begin(); kit != mapCube.keyWord().end(); ++kit) {
               //Find the value type
               CCfits::ValueType vtype = kit->second->keytype();
               std::string value;
               std::ostringstream os;
               switch(vtype) {
                  case CCfits::Tstring:
                     kit->second->value(value);
                     break;

                  case CCfits::Tlogical:
                     bool btmp;
                     kit->second->value(btmp);
                     os << btmp;
                     value = os.str();
                     break;

                  case CCfits::Tbyte:
                     char ctmp;
                     kit->second->value(ctmp);
                     os << ctmp;
                     value = os.str();
                     break;

                  case CCfits::Tshort:
                     short stmp;
                     kit->second->value(stmp);
                     os << stmp;
                     value = os.str();
                     break;

                  case CCfits::Tushort:
                     unsigned short ustmp;
                     kit->second->value(ustmp);
                     os << ustmp;
                     value = os.str();
                     break;

                  case CCfits::Tint:
                     int itmp;
                     kit->second->value(itmp);
                     os << itmp;
                     value = os.str();
                     break;

                  case CCfits::Tuint:
                     unsigned int uitmp;
                     kit->second->value(uitmp);
                     os << uitmp;
                     value = os.str();
                     break;

                  case CCfits::Tlong:
                     long ltmp;
                     kit->second->value(ltmp);
                     os << ltmp;
                     value = os.str();
                     break;

                  case CCfits::Tulong:
                     unsigned long ultmp;
                     kit->second->value(ultmp);
                     os << ultmp;
                     value = os.str();
                     break;

                  case CCfits::Tlonglong:
                     long long lltmp;
                     kit->second->value(lltmp);
                     os << lltmp;
                     value = os.str();
                     break;

                  case CCfits::Tfloat:
                     float ftmp;
                     kit->second->value(ftmp);
                     os << ftmp;
                     value = os.str();
                     break;

                  case CCfits::Tdouble:
                     double dtmp;
                     kit->second->value(dtmp);
                     os << dtmp;
                     value = os.str();
                     break;

                  default:
                     value = "";
               }

               fKeywords[kit->first] = value;
            }

         }
      }

      /**\brief Write the skymap to a file
       *
       * \param fileName is the name of the file to write to
       *
       * The file is overwritten without warning.
       */
      void write(const std::string & fileName, const std::string type="energy", const std::string unit="MeV", const std::map<std::string, std::string> &additionalFitsKeywords = (std::map<std::string, std::string>())) const {
         //Do nothing if there is no data
	if (fdimension.Spectra().size() == 0) return;

         //Prepend ! so the file gets overwritten
         std::string fname = fileName[0] != '!' ? "!" + fileName : fileName;
         //Append .gz to allow compression
#ifdef ENABLE_COMPRESSION
         if ( fname.substr(fname.size()-3).compare(".gz") )
            fname += ".gz";
#endif

         //Create a CCFITS object with the filename
         CCfits::FITS fits(fname, CCfits::Write );

         //Create a map of typeid's to cfitsio format characters
         std::map<const char*, std::string> formatMap;
         formatMap[typeid(char).name()] = "I";
         formatMap[typeid(short).name()] = "I";
         formatMap[typeid(int).name()] = "J";
         formatMap[typeid(long).name()] = "K";
         formatMap[typeid(float).name()] = "E";
         formatMap[typeid(double).name()] = "D";
         formatMap[typeid(unsigned char).name()] = "B";
         formatMap[typeid(unsigned short).name()] = "U";
         formatMap[typeid(unsigned int).name()] = "V";
         formatMap[typeid(unsigned long).name()] = "K";
         //Create the data table to store the actual image.  The extension name is
         //SKYMAP
         std::ostringstream strForm;
         strForm << 1 << formatMap[typeid(T).name()];
         std::vector<std::string> colNames(fdimension.Spectra().size(),"Spectra")
            , colForm(fdimension.Spectra().size(),strForm.str())
            , colUnit(fdimension.Spectra().size(),"");

         for (size_t ibin(0); ibin < fdimension.Spectra().size(); ++ibin) {
            std::ostringstream bname;
            bname << "Bin "<< ibin;
            colNames[ibin] = bname.str();
         }
         CCfits::Table *imgTable = fits.addTable("SKYMAP2", Npix(), colNames, colForm, colUnit);

         /* This is not useful for gzipped files, which is our default.
         //Write the data, split it into efficient row numbers
         long nrows(10);
         int status(0);
         fits_get_rowsize(fits.fitsPointer(),&nrows,&status);
         //default to old method if this doesn't work for some reason
         if (nrows < 1)
            nrows = Npix();
         for (size_t ipix(0); ipix < size_t(Npix()); ipix += nrows) {
            const size_t start = ipix*fdimension.Spectra().size();
            const size_t size = std::min(size_t(nrows),size_t(Npix())-ipix);
            for (size_t ibin(0); ibin < fdimension.Spectra().size(); ++ibin) {
               std::ostringstream bname;
               bname << "Bin "<< ibin;
               std::valarray<T> binMap = fMap[std::slice(start+ibin,size,fdimension.Spectra().size())];
               imgTable->column(bname.str()).write(binMap,ipix+1);
            }
         }
         */
         for (size_t ibin(0); ibin < fdimension.Spectra().size(); ++ibin) {
            std::ostringstream bname;
            bname << "Bin "<< ibin;
            std::valarray<T> binMap = fMap[std::slice(ibin,Npix(),fdimension.Spectra().size())];
            imgTable->column(bname.str()).write(binMap,1);
         }


         //Write keywords
         imgTable->addKey("PIXTYPE", "HEALPIX", "Healpix pixel scheme");
         std::string keyvalue = "RING";
         if (NEST == Order()) keyvalue = "NESTED";
         std::ostringstream osst;
         imgTable->addKey("ORDERING", keyvalue, "Ring or nested ordering of pixels");
         imgTable->addKey("NSIDE", Nside(), "Number of sides in a base pixel");
         imgTable->addKey("FIRSTPIX", 0, "Number of first pixel");
         imgTable->addKey("LASTPIX", Npix()-1, "Number of last pixel");
         osst << "Number of "<<type<<" bins in spectra";
         imgTable->addKey("NBRBINS", nSpectra(), osst.str());
         imgTable->addKey("COORDTYPE", "GAL", "");

         //Write additional keywords
         std::map<std::string, std::string>::const_iterator it;
         for (it = additionalFitsKeywords.begin(); it != additionalFitsKeywords.end(); ++it) 
            imgTable->addKey(it->first, it->second, "");

         //Create another table to store the energy
         if (fdimension.IsBinned()){
            colNames.resize(3); colUnit.resize(3); colForm.resize(3);
            colNames[0] = "CHANNEL"; colNames[1] = type + "_MIN"; colNames[2] = type + "_MAX";
            colUnit[0] = ""; colUnit[1] = unit; colUnit[2] = unit;
            colForm[0] = "I"; colForm[1] = "D"; colForm[2] = "D";
            CCfits::Table *eTable = fits.addTable("EBOUNDS", nSpectra(), colNames, colForm, colUnit);
            //Create the channel array
            std::valarray<int> channel(fdimension.Spectra().size());
            for (int i = 0; i < int(channel.size()); ++i){
               channel[i] = i+1;
            }
            eTable->column("CHANNEL").write(channel,1);
            eTable->column(type + "_MIN").write(fdimension.SpectraMin(),1);
            eTable->column(type + "_MAX").write(fdimension.SpectraMax(),1);
         }else{
            colNames.resize(1); colUnit.resize(1); colForm.resize(1);
            colNames[0] = type;
            colUnit[0] = unit;
            colForm[0] = "D";
            CCfits::Table *eTable = fits.addTable("ENERGIES", nSpectra(), colNames, colForm, colUnit);

            eTable->column(type).write(fdimension.Spectra(),1);
         }
      }

      //! Return the keywords read from the FITS file.  They are converted to strings.
      const std::map<std::string,std::string> getKeywords() { return fKeywords; }

      /**\brief Convert all the pixels to a different type and returns the new copy of the map
       *
       * \param dummy controls the output type
       */
      template <typename C>
         Skymap<C> convert(C dummy) const{
            //Create an output, the same size and same spectra
            Skymap<C> output;
            if (fdimension.IsBinned()) {
	      output.Resize(order_, fdimension.SpectraMin(), fdimension.SpectraMax(), scheme_);
            } else {
	      output.Resize(order_, fdimension.Spectra(), scheme_);
            }
            for (int i = 0; i < Npix(); ++i){
               for (int j = 0; j < nSpectra(); ++j){
                  output[i][j] = C((*this)[i][j]);
               }
            }
            return output;
         }

      /** \brief Interpolate the map to a different order.
       *
       * \param order specifies the new order.  It makes little sense
       * to interpolate to lower order, but it is not forbidden.
       *
       * Internally we use the get_interpol method of Healpix_Base to
       * calculate the interpolation value at the centers of the new
       * pixels
       */
      Skymap<T> interpolate(const int order, const bool includeZeros=true) const{
         //If the order is the same, just return
         if (order == order_) return *this;

         //Create a new skymap with the new order, keeping the spectra and the
         //scheme_
         Skymap<T> newMap;
         if (fdimension.IsBinned()) {
	   newMap.Resize(order, fdimension.SpectraMin(), fdimension.SpectraMax(), scheme_);
         }else{
	   newMap.Resize(order, fdimension.Spectra(), scheme_);
         }

         //Loop the new map and calculate the interpolation
         fix_arr<int,4> pixels;
         fix_arr<double,4> weight;
         for (int i = 0; i < newMap.Npix(); ++i) {
            get_interpol(newMap.pix2ang(i),pixels,weight);
            if ( includeZeros ) {
               for (int j = 0; j < 4; ++j){
		 for (size_t k = 0; k < fdimension.Spectra().size(); ++k){
                     newMap[i][k] += weight[j]*(*this)[pixels[j]][k];
                  }
               }
            } else {
               fix_arr<bool,4> include;
               unsigned int included(0);
               for (int j = 0; j < 4; ++j){
                  include[j] = false;
                  for (size_t k = 0; k < fdimension.Spectra().size(); ++k){
                     if ( (*this)[pixels[j]][k] != 0 ) {
                        include[j] = true;
                        ++included;
                        break;
                     }
                  }
               }
               if (included) {
                  //Fix the weights before adding the pixels
                  for (int j = 0; j < 4; ++j){
                     if (! include[j]) {
                        for (int jj = 0; jj < 4; ++jj) {
                           if (include[jj]) {
                              weight[jj] += weight[j]/included;
                           }
                        }
                        weight[j] = 0;
                     }
                  }
                  for (int j = 0; j < 4; ++j){
                     if (include[j]) {
		       for (size_t k = 0; k < fdimension.Spectra().size(); ++k){
                           newMap[i][k] += weight[j]*(*this)[pixels[j]][k];
                        }
                     }
                  }
               }
            }
         }
         return newMap;
      }

      /** \brief Calculate an md5 checksum, using the data and binning information
       *
       * Does not check auxillary information such as keywords.
       */
      std::string md5() const {
         MD5 cs;
         cs.update(reinterpret_cast<const unsigned char*>(&fMap[0]), fMap.size()*sizeof(T)/sizeof(char));
         cs.update(reinterpret_cast<const unsigned char*>(&fdimension.Spectra()[0]), fdimension.Spectra().size()*sizeof(double)/sizeof(char));
         cs.update(reinterpret_cast<const unsigned char*>(&fdimension.SpectraMin()[0]), fdimension.SpectraMin().size()*sizeof(double)/sizeof(char));
         cs.update(reinterpret_cast<const unsigned char*>(&fdimension.SpectraMax()[0]), fdimension.SpectraMax().size()*sizeof(double)/sizeof(char));
         cs.finalize();
         return cs.hexdigest();
      }


      /** \brief Rebin the map to a different order.
       *
       * \param order specifies the new order.  If it is greater than the old
       * order, the map size is increased and the new map contains multiple
       * pixels with the same values.  If it is less than the old order, the new
       * pixel values will be the average of the old pixels contained in the new
       * set.
       *
       * It is not wise to rebin the map to a higher order, rather use the
       * coordinate selector to get the pixel values.  That is what this
       * rebinning does anyway.
       */
      Skymap<T> rebin(int order, bool SAcorrect = true) const{
         //If the order is the same, just return
         if (order == order_) return *this;

         //Create a new skymap with the new order, keeping the spectra and the
         //scheme_
         Skymap<T> newMap;
         if (fdimension.IsBinned()) {
	   newMap.Resize(order, fdimension.SpectraMin(), fdimension.SpectraMax(), scheme_);
         }else{
	   newMap.Resize(order, fdimension.Spectra(), scheme_);
         }

         //What we do now depends if the order, if it is greater, loop the new map
         //and set the pixel from their coordinate
         if (order > order_){
            int dOrder = order - order_;
            const int nPix = (1<<dOrder)*(1<<dOrder);
            for (Iterator it = newMap.begin(); it != newMap.end(); ++it){
               SM::Coordinate co = it.coord();
               *it = (*this)[co];
            }
            if (!SAcorrect)
               newMap /= double(nPix);
         }else{
            //If the order is less, we have to average over several pixels
            //The difference in the number of pixels
            int dOrder = order_ - order;
            const int nPix = (1<<dOrder)*(1<<dOrder);
            //Create a memory to work with and an
            //ArraySlice
#pragma omp parallel default(shared)
            {
	      std::valarray<T> avStore(0.0, fdimension.Spectra().size());
	      ArraySlice<T> average(&avStore[0], fdimension.Spectra().size());
#pragma omp for schedule(static)
               for (long p = 0; p < newMap.Npix(); ++p){
                  average = 0;
                  //Use pixel numbering from the nested scheme
                  const int pp = (newMap.Scheme() == NEST) ? p : newMap.ring2nest(p);
                  for (long po = pp*nPix; po < (pp+1)*nPix; ++po){
                     //Convert pixel back to correct scheme if needed
                     const int pop = (Scheme() == NEST) ? po : nest2ring(po);
                     average += (*this)[pop];
                  }
                  //Divide by the number of pixels to get the average
                  if (SAcorrect) average /= double(nPix);
                  newMap[p] = average;
               }
            }//End parallel
         }
         //Assign the newMap to us
         return newMap;
      }

      /**\brief The size of the spectra */
      int nSpectra() const{ return fdimension.Spectra().size();}

      /**\brief Return the dimension of the map */
      const SkymapDim & getDimension() const { return fdimension;}

      /**\brief Set the values at which the spectra is evaluated 
       *
       * \param spectra are the new values, stored in a valarray
       *
       * \return if the size of the input values do not confirm with the spectral size, a false value is returned.
       *
       * This method changes the skymap to unbinned mode
       */
      bool setSpectra(const std::valarray<double> &spectra){
	if (spectra.size() == fdimension.Spectra().size()){
	  fdimension.Spectra() = spectra;
	  fdimension.SetBinned(false);
	  return true;
	}
	return false;
      }

      /**\brief Set the values at which the spectra is evaluated in
       * binned mode
       *
       * \param specMin are the new lower boundaries
       * \param specMax are the new upper boundaries
       *
       * \return if the size of the input values do not confirm with the spectral size, a false value is returned.
       *
       * This method changes the skymap to binned mode
       */
      bool setSpectra(const std::valarray<double> &specMin, const std::valarray<double> &specMax){
	if (specMin.size() == fdimension.Spectra().size() && specMax.size() == fdimension.Spectra().size()){
	  fdimension.SpectraMin().resize(fdimension.Spectra().size());
	  fdimension.SpectraMax().resize(fdimension.Spectra().size());
	  fdimension.Spectra() = 0.5*(specMin+specMax);
	  fdimension.SpectraMin() = specMin;
	  fdimension.SpectraMax() = specMax;
	  fdimension.SetBinned(true);
	  return true;
	}
	return false;
      }

      /**\brief Set the values at which the spectra is evaluated
       *
       * \param spectra are the new values in a C array
       * \param size is the size of the new array.
       *
       * \return if the size of the input values do not confirm with the spectral size, a false value is returned.
       */
      bool setSpectra(const double spectra[], const size_t size){
	if (size == fdimension.Spectra().size()){
	  for (size_t i=0; i<fdimension.Spectra().size(); ++i){
	    fdimension.Spectra()[i] = spectra[i];
            }
            return true;
         }
         return false;
      }

      /**\brief Return the values at which the spectra is evaluated*/
      const std::valarray<double> & getSpectra() const{return fdimension.Spectra();}

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
      
      /** \brief Fill skymap with data from CAR map
       * 
       * \param map is an one dimensional array formatted in mapcube format (order, l,b,spectra)
       * \param l are the longitude values in the grid, assumed to be uniformly spaced, must have at least 2 values
       * \param b are the latitude values in the grid, assumed to be uniformly spaced, must have at least 2 values
       * \param spectra are the spectral values, throws an error if the size of the spectra is not identical to that of the map
       * \param correctSA is a boolean which controls wether we correct for solid angle or not.
       */
      void filllbCARmap( const std::valarray<double> &map, const std::valarray<double> &l, const std::valarray<double> &b, const std::valarray<double> &spectra, bool correctSA=false) {
	//Make sure the spectra are the same size
	if (nSpectra() != int(spectra.size())) {
	  std::cerr<<"Spectral size do not match when filling skymap with mapcube data"<<std::endl;
	  throw(1);
	}
	
	//Set up crval, cdelt, and crpix and pass it on
	std::valarray<double> crval(0.0, 3), cdelt(1.0, 3), crpix(1.0, 3);
	std::valarray<long> axis(3);
	
	cdelt[0] = l[1]-l[0];
	cdelt[1] = b[1]-b[0];
	crval[0] = l[0];
	crval[1] = b[0];
	axis[0] = l.size();
	axis[1] = b.size();
	axis[2] = spectra.size();
	
	fillMapcube(map, axis, crval, cdelt, crpix, correctSA);
      }

      /**\brief Fill skymap with data from an l and b map
       *
       * \param map is an one dimensional array formatted in mapcube format (order l,b,spectra)
       * \param nl the number of l bins
       * \param nb the number of b bins
       * \param correctSA is a boolean which controls wether we correct for solid angle or not.
       * 
       * The spectral size is assumed to be the same as the skymap and no checking is performed.
       * This method splits the healpix map into finer bins and selects the undelying l and b pixel.
       * The input map must be full sky and the first pixel in the input map must have the edge at 
       * l = 0 and b = -90.  It is assumed that the pixel values are given at the center of the pixels,
       * so the pixel (ib, il) spans the area  l = [il*360/nl, (il+1)*360/nl] and b = [ib*180/nl-90, (ib+1)*180/nl-90].
       * Note that ib and il are numbered from 0.
       */
      void filllbCARarray( const T map[], const int nl, const int nb, const bool correctSA=false){
         //Set up the crval, cdelt, crpix and pass it on
         std::valarray<double> crval(0.0, 3), cdelt(1.0,3), crpix(1.0,3);
         std::valarray<long> axis(3);

         //Full skymaps
         cdelt[0] = 360./nl;
         cdelt[1] = 180./nb;
         crval[0] = cdelt[0]/2.;
         crval[1] = -90+cdelt[1]/2.;
         axis[0] = nl;
         axis[1] = nb;
         axis[2] = nSpectra();

         //The valarray for the image data
         std::valarray<T> image(map,axis[0]*axis[1]*axis[2]);

         fillMapcube(image,axis,crval,cdelt,crpix,correctSA);
      }

      /**\bried Convert to a mapcube in CAR projection in galactic
       * coordinates and write it to a file
       *
       * Spatial coordinates are given as in the fits standard, with
       * vectors of length 2: NAXIS, CRVAL, CDELT, CRPIX.  The first
       * element corresponds to longitude but the second to latitude.
       *
       * If NAXIS is 0, full sky is assumed, trying to take CDELT
       * and CRVAL into account.  If CDELT is also 0, the resolution
       * of the skymap is used.  We put a pixel at 0,0 by default.
       * 
       * \param fileName is the name of the output file, with full
       * path.
       * \param NAXIS is the number of pixels per axis (length 2)
       * \param CRVAL is the value of pixel CRPIX (length 2)
       * \param CDELT is the difference in value between two pixels
       * (length 2)
       * \param CRPIX is the pixel number that is at CRVAL (note that
       * pixels are numbered from 1 like when you count, not from 0
       * like in C/C++)
       * \param correctSA should be true if you are working with
       * counts, false if you are working with flux.
       *
       * \return a 1D valarray that contains all the pixels, where l
       * loops fastest, then b and last the energy.
       *
       * Note that if the length of NAXIS et al. is greater than 2,
       * the values are simply ignored.  There is no way to rebin in
       * energy.  On output, the arrays will be resized to length 3
       * and the corresponding value for the energy dimension added
       * to it.  It will be in log energy.
       *
       * This routine works by splitting the output pixels up into
       * 1024 smaller pixels and assigning it the value that is
       * directly under the center.  Then the final bigger pixels is
       * achieved by taking the average.
       */
      void writeMapcube(const std::string &fileName, std::vector<long> &NAXIS, std::vector<double> &CRVAL, std::vector<double> &CDELT, std::vector<double> &CRPIX, bool correctSA = false) const{
         //Check for consistency in the input parameters
         NAXIS.resize(2,0);
         CDELT.resize(2,0.0);
         CRVAL.resize(2,0.0);
         CRPIX.resize(2,1.0);
         double res = resolution();
         if (NAXIS[0] == 0) {
            if (CDELT[0] == 0) {
               NAXIS[0] = long(360/res);
               CDELT[0] = -360./double(NAXIS[0]);
            } else {
               NAXIS[0] = long(360/fabs(CDELT[0]));
            }
            CRVAL[0] = 0;
            CRPIX[0] = (NAXIS[0]+1)/2.;
         }
         if (NAXIS[1] == 0) {
            if (CDELT[1] == 0) {
               NAXIS[1] = long(180/res);
	       if (NAXIS[1] % 2 == 0) {
                 CDELT[1] = 180./double(NAXIS[1]);
                 ++NAXIS[1];
	       } else {
                 CDELT[1] = 180./double(NAXIS[1]-1);
	       }
	    } else {
               NAXIS[1] = long(180/fabs(CDELT[1]))+1;
            }
            CRVAL[1] = 0;
            CRPIX[1] = (NAXIS[1]+1)/2.;
         }
         //Fix CDELT, if it is still 0
         if (CDELT[0] == 0) {
            CDELT[0] = -360./double(NAXIS[0]);
         }
         if (CDELT[1] == 0) {
            CDELT[1] = 180./double(NAXIS[1]);
         }

         //Do the conversion
         std::valarray<T> mapCube = toMapcube(NAXIS,CRVAL,CDELT,CRPIX,correctSA);

         //Create the fits file
         CCfits::FITS fits(fileName, DOUBLE_IMG, 3, &NAXIS[0]);

         fits.pHDU().write(1, NAXIS[0]*NAXIS[1]*NAXIS[2], mapCube);

         //Write the keywords
         fits.pHDU().addKey("CRVAL1", CRVAL[0], "Value of longitude in pixel CRPIX1");
         fits.pHDU().addKey("CDELT1", CDELT[0], "Step size in longitude");
         fits.pHDU().addKey("CRPIX1", CRPIX[0], "Pixel that has value CRVAL1");
         fits.pHDU().addKey("CTYPE1", "GLON", "The type of parameter 1 (Galactic longitude in CAR projection)");
         fits.pHDU().addKey("CUNIT1", "deg", "The unit of parameter 1");
         fits.pHDU().addKey("CRVAL2", CRVAL[1], "Value of latitude in pixel CRPIX2");
         fits.pHDU().addKey("CDELT2", CDELT[1], "Step size in latitude");
         fits.pHDU().addKey("CRPIX2", CRPIX[1], "Pixel that has value CRVAL2");
         fits.pHDU().addKey("CTYPE2", "GLAT", "The type of parameter 2 (Galactic latitude in CAR projection)");
         fits.pHDU().addKey("CUNIT2", "deg", "The unit of parameter 2");
         fits.pHDU().addKey("CRVAL3", CRVAL[2], "Energy of pixel CRPIX3");
         fits.pHDU().addKey("CDELT3", CDELT[2], "log10 of step size in energy (if it is logarithmically distributed, see ENERGIES/EBOUNDS extension)");
         fits.pHDU().addKey("CRPIX3", CRPIX[2], "Pixel that has value CRVAL3");
         fits.pHDU().addKey("CTYPE3", "Energy", "Axis 3 is the spectra");
         fits.pHDU().addKey("CUNIT3", "MeV", "The unit of axis 3");

         //Write the energy extension
         if (fdimension.IsBinned()){
            std::vector<std::string> colNames(3), colForm(3), colUnit(3);
            colNames[0] = "CHANNEL"; colNames[1] = "E_MIN"; colNames[2] = "E_MAX";
            colUnit[0] = ""; colUnit[1] = ""; colUnit[2] = "";
            colForm[0] = "I"; colForm[1] = "D"; colForm[2] = "D";
            std::valarray<int> channel(nSpectra());
            for (int i = 0; i < channel.size(); ++i){
               channel[i] = i+1;
            }
            CCfits::Table *eTable = fits.addTable("EBOUNDS", nSpectra(), colNames, colForm, colUnit);
            eTable->column("CHANNEL").write(channel,1);
            eTable->column("E_MIN").write(fdimension.SpectraMin(),1);
            eTable->column("E_MAX").write(fdimension.SpectraMax(),1);
         }else{
            std::vector<std::string> colNames(1), colForm(1), colUnit(1);
            colNames[0] = "Energy";
            colUnit[0] = "";
            colForm[0] = "D";
            CCfits::Table *eTable = fits.addTable("ENERGIES", nSpectra(), colNames, colForm, colUnit);

            eTable->column("Energy").write(fdimension.Spectra(),1);
         }

      }

      /**\bried Convert to a mapcube in CAR projection in galactic
       * coordinates
       *
       * Spatial coordinates are given as in the fits standard, with
       * vectors of length 2: NAXIS, CRVAL, CDELT, CRPIX.  The first
       * element corresponds to longitude but the second to latitude.
       * No effort is put into correcting for input errors.
       *
       * \param NAXIS is the number of pixels per axis (length 2)
       * \param CRVAL is the value of pixel CRPIX (length 2)
       * \param CDELT is the difference in value between two pixels
       * (length 2)
       * \param CRPIX is the pixel number that is at CRVAL (note that
       * pixels are numbered from 1 like when you count, not from 0
       * like in C/C++)
       * \param correctSA should be true if you are working with
       * counts, false if you are working with flux.
       *
       * \return a 1D valarray that contains all the pixels, where l
       * loops fastest, then b and last the energy.
       *
       * Note that if the length of NAXIS et al. is greater than 2,
       * the values are simply ignored.  There is no way to rebin in
       * energy.  On output, the arrays will be resized to length 3
       * and the corresponding value for the energy dimension added
       * to it.  It will be in log energy.
       *
       * This routine works by splitting the output pixels up into
       * 1024 smaller pixels and assigning it the value that is
       * directly under the center.  Then the final bigger pixels is
       * achieved by taking the average.
       */
      std::valarray<T> toMapcube(std::vector<long> &NAXIS, std::vector<double> &CRVAL, std::vector<double> &CDELT, std::vector<double> &CRPIX, bool correctSA = false) const{
         //Create the output array
	std::valarray<T> output(NAXIS[0]*NAXIS[1]*fdimension.Spectra().size());

	//Calculating the finer resolution, taking the resolution of
	//the current HEALPix map into account, not making it too
	//small
	const double res = resolution();
	const int nfiner = 32; //The nominal number of subdivision, unless the pixel size of the HEALPix map is small enough
	const int nlfiner = fabs(CDELT[0]) < res ? int(fabs(CDELT[0])/res*nfiner)+1 : nfiner;
	const int nbfiner = fabs(CDELT[1]) < res ? int(fabs(CDELT[1])/res*nfiner)+1 : nfiner;
	int npfiner = nlfiner*nbfiner;
	
	//The finer resolution
	const double cdeltl = CDELT[0]/nlfiner;
	const double cdeltb = CDELT[1]/nbfiner;
	
	//Add the spectral information to the vectors
	NAXIS.resize(3);
	CRVAL.resize(3);
	CDELT.resize(3);
	CRPIX.resize(3);
	NAXIS[2] = fdimension.Spectra().size();
	CRVAL[2] = fdimension.Spectra()[0];
	CRPIX[2] = 1;
	CDELT[2] = fdimension.Spectra().size() > 1 ? log10(fdimension.Spectra()[1]/fdimension.Spectra()[0]) : 0;

         //Loop over the output pixels, dividing them up into the
         //finer grid
#pragma omp parallel for default(shared) schedule(static) private(npfiner)
         for (long ib = 0; ib < NAXIS[1]; ++ib){
            //Create storage to calculate the average
	   std::valarray<T> totStore(T(0),fdimension.Spectra().size());
	   ArraySlice<T> totSpectra(&totStore[0], fdimension.Spectra().size());
	   const size_t indb = ib*NAXIS[0];
	   //The value of b at the edge of the pixel
	   const double b = CRVAL[1] + CDELT[1]*(ib+0.5-CRPIX[1]);
	   //Correct for solid angle if necessary
	   double lbSolidAngle = 1;
	   if (correctSA){
	     lbSolidAngle = fabs(CDELT[0])*utl::kPi/180. * (sin(utl::kPi*(b+fabs(CDELT[1]))/180.) - sin(utl::kPi*b/180.));
	   }
	   for (long il = 0; il < NAXIS[0]; ++il){
	     totSpectra = 0;
	     const size_t indl = indb + il;
	     //The value of l at the edge of the pixel
	     const double l = CRVAL[0] + CDELT[0]*(il+0.5-CRPIX[0]);
	     //Loop over the finer grid and calculate the sum
	     npfiner = 0;
	     for (int ibf = 0; ibf < nbfiner; ++ibf){
	       const double bf = b + (ibf+0.5)*cdeltb;
	       if (fabs(bf) >= 90) continue;
	       for (int ilf = 0; ilf < nlfiner; ++ilf){
		 double lf = l + (ilf+0.5)*cdeltl;
		 while (lf < 0) lf += 360;
		 while (lf >= 360) lf -= 360;
		 totSpectra += (*this)[SM::Coordinate(lf,bf)];
		 ++npfiner;
	       }
	     }
	     //Now assign the average to the pixels
	     for (int j = 0; j < totSpectra.size(); ++j){
	       output[j*NAXIS[0]*NAXIS[1]+indl] = totSpectra[j]/double(npfiner)*lbSolidAngle;
	     }
	   }
         }
	 
         //Correct for solid Angle if needed
         if (correctSA)
	   output/=solidAngle();
	 
         return output;
      }
      
      /**\brief Convert to a l and b map in CAR projection
       *
       * \param nl is the number of l values in the output array
       * \param nb is the number of b values in the output array
       * \param correctSA is a boolean which tells the routine weather to correct for solid angle or not.
       * (When converting counts map, correctSA should be true.)
       *
       * \return a pointer to a dynamically allocated one dimensional c array compatible with mapcubes.  The 
       * l values are looped fastest, then b values and spectra last.  Please deallocate the array once done using
       * it to prevent memory leaks.
       */
      T * lbCARarray( const int nl, const int nb, double CRVAL[3], int CRPIX[3], double CDELT[3], const bool correctSA = false) const{
         T *array;
         array = new T[nl*nb*fdimension.Spectra().size()];
         const int nfiner = 32;
         const double lres = 360./nl;
         const double bres = 180./nb;
         const double mapres = resolution();
         //The resolution of the conversion has to be 1/10th of the resolution of
         //the resulting map.  But it only has to have 1/10th of the resolution of
         //the healpix map.
         const int nlfiner = lres < mapres ? int(lres/mapres*nfiner)+1 : nfiner;
         const int nbfiner = bres < mapres ? int(bres/mapres*nfiner)+1 : nfiner;
         const int nboxes = nbfiner*nlfiner;
         //Populate the fits header arrays
         CRVAL[0] = 0; CRVAL[1] = -90+bres/2.; CRVAL[2] = fdimension.Spectra()[0];
         CRPIX[0] = CRPIX[1] = CRPIX[2] = 1; //fits starts counting at 1
         CDELT[0] = lres; CDELT[1] = bres; CDELT[2] = fdimension.Spectra().size() > 1 ? fdimension.Spectra()[1] - fdimension.Spectra()[0] : 0;
         //Divide by the solid angle if needed to get the correct results
         double hpSolidAngle = 1;
         if (correctSA){
            hpSolidAngle = solidAngle();
         }
         //Create storage
         std::valarray<T> totStore(T(0),fdimension.Spectra().size());
         ArraySlice<T> totSpectra(&totStore[0], fdimension.Spectra().size());
#pragma omp parallel for default(shared) schedule(static)
         for (int j = 0; j<nb; ++j){
            const int ind1 = j*nl;
            double bb = bres * j - 90;
            for (int i = 0; i<nl; ++i){
               totSpectra = 0;
               int ind2 = ind1 + i;
               //First pixel's center is at 0
               double ll = lres*i - 0.5*lres;
               //Correct for solid angle if necessary
               double lbSolidAngle = 1;
               if (correctSA){
                  //Note that we use cos, since sin(x-pi/2) = -cos(x)
                  lbSolidAngle = lres * utl::kPi/180 * (cos(utl::kPi/double(nb)*j)-cos(utl::kPi/double(nb)*(j+1))); //AWS20081020
               }
               for (int jj = 0; jj<nbfiner; ++jj){
                  const double b = bb + (jj+0.5)*bres/nbfiner;
                  for (int ii = 0; ii<nlfiner; ++ii){
                     double l = ll + (ii+0.5)*lres/nlfiner;
                     if (l < 0) l += 360;
                     totSpectra += (*this)[SM::Coordinate(l,b)];
                  }
               }
               for (size_t ll = 0; ll<fdimension.Spectra().size(); ++ll){
                  int ind = ind2 +ll*nl*nb;
                  array[ind] = totSpectra[ll]/double(nboxes)/ hpSolidAngle * lbSolidAngle;
               }
            }
         }
         return array;
      }

      /**\brief Fill the skymap with data from a Healpix_Map
       *
       * \param iSpectra the point in spectra at which the map should be inserted
       * \param inmap is the Healpix_Map to insert
       */
      void fromHealpixMap(const int iSpectra, const Healpix_Map<T> & inmap){
         //Check for iSpectra bounds and the order of the map
	if (iSpectra < int(fdimension.Spectra().size()) && order_ == inmap.Order() && scheme_ == inmap.Scheme()){
#if (_OPENMP >= 201307)
#pragma omp simd
#endif
            for(int i=0; i<npix_; ++i){
	      fMap[i*fdimension.Spectra().size()+iSpectra] = inmap[i];
            }
         } else {
            std::cerr<<"Failed to insert Healpix_Map"<<std::endl;
         }
      }

      /**\brief Convert the data to a healpix map
       *
       * \param iSpectra the point in spectra for the returned map
       * \return the Healpix map at the corresponding point in spectra
       */
      Healpix_Map<T> toHealpixMap(const int iSpectra) const{
         //Create the output map and fill it with zeros
         Healpix_Map<T> outmap(order_, scheme_);

         //iSpectra must be within bounds, otherwise return a zero map
         if (iSpectra < int(fdimension.Spectra().size())){
            for(int i=0; i<npix_; ++i){
	      outmap[i] = fMap[i*fdimension.Spectra().size()+iSpectra];
            }
         } else {
            outmap.fill(0.0);
         }

         return outmap;
      }

      /**\brief Reference to a spectra given a coordinate
       *
       * \param coordinate is an instance of the coordinate class \see SM::Coordinate
       * \return a valarray of the spectra in the pixel corresponding to coordinate.
       */
      ArraySlice<T> operator [] (const SM::Coordinate & coordinate){
         return (*this)[ang2pix(coordinate.healpixAng())];
      }
      /** \overload */
      const ArraySlice<T> operator [] (const SM::Coordinate & coordinate) const{
         return (*this)[ang2pix(coordinate.healpixAng())];
      }
      /**\brief Reference to a spectra given a pixel
       *
       * \param pixel is a pixel number for the skymap
       * \return a std::valarray of the spectra in the pixel
       */
      ArraySlice<T> operator [] (const int pixel){
	return ArraySlice<T>(&fMap[pixel*fdimension.Spectra().size()], fdimension.Spectra().size());
      }
      /** \overload */
      const ArraySlice<T> operator [] (const int pixel) const{
	return ArraySlice<T>(&fMap[pixel*fdimension.Spectra().size()], fdimension.Spectra().size());
      }


      /** \brief Bidirectional iterator for the skymap */
      class Iterator : public std::iterator<std::bidirectional_iterator_tag, std::vector<T> > {
         private:
            int m_index; //!< The index for the array
            Skymap<T> & m_map; //!< A reference to a skymap
            Iterator(){} //!< The default iterator constructor should not be allowed
         public:
            /**\brief Constructor */
            Iterator(const int index, Skymap<T> & map) : m_index(index), m_map(map) {};
            /**\brief Pre-increment operator */
            Iterator & operator ++ () { ++m_index; return *this; }
            /**\brief Post-increment operator */
            Iterator operator ++ (int i) { Iterator tmp(*this); ++m_index; return tmp; }
            /**\brief Pre-decrement operator */
            Iterator & operator -- () { --m_index; return *this; }
            /**\brief Post-decrement operator */
            Iterator operator -- (int i) { Iterator tmp(*this); --m_index; return tmp; }
            /**\brief Dereference operator */
            ArraySlice<T> operator * () { return m_map[m_index];}
            /**\brief Unequal operator */
            bool operator != (const Iterator compare) { return m_index != compare.m_index; }
            /**\brief Return the corresponding coordinate */
            SM::Coordinate coord () { SM::Coordinate co(m_map.pix2ang(m_index)); return co; }
      };
      /** \brief Constant bidirectional iterator for the skymap */
      class constIterator : public std::iterator<std::bidirectional_iterator_tag, std::vector<T> > {
         private:
            int m_index; //!< The index for the array
            const Skymap<T> & m_map; //!< A reference to a skymap
            constIterator(){} //!< The default iterator constructor should not be allowed
         public:
            /**\brief Constructor*/
            constIterator(const int index, const Skymap<T> & map) : m_index(index), m_map(map) {};
            /**\brief Pre-increment operator */
            constIterator & operator ++ () { ++m_index; return *this; }
            /**\brief Post-increment operator */
            constIterator operator ++ (int i) { Iterator tmp(*this); ++m_index; return tmp; }
            /**\brief Pre-decrement operator */
            constIterator & operator -- () { --m_index; return *this; }
            /**\brief Post-decrement operator */
            constIterator operator -- (int i) { Iterator tmp(*this); --m_index; return tmp; }
            /**\brief Dereference operator */
            const ArraySlice<T> operator * () { return m_map[m_index];}
            /**\brief Unequal operator */
            bool operator != (const constIterator compare) { return m_index != compare.m_index; }
            /**\brief Return the corresponding coordinate */
            SM::Coordinate coord () { SM::Coordinate co(m_map.pix2ang(m_index)); return co; }
      };


      /**\brief Returns an Iterator for the first pixel */
      Iterator begin() { return Iterator(0, *this);}
      /** \overload */
      constIterator begin() const { return constIterator(0, *this);}
      /**\brief Returns an Iterator to the last+1 pixel */
      Iterator end() { return Iterator(npix_, *this); }
      /** \overload */
      constIterator end() const { return constIterator(npix_, *this); }

      /**\brief Convert the skymap from one coordinate system to the
       * other
       *
       * \parameter from is the coordinate system the map is
       * currently in
       * \parameter to is the coordinate system to transfer to
       * \parmeter mj is the modified Julian date
       *
       * \return Skymap converted to the new coordinate system
       *
       * Uses similar methods as transformation to a flat CAR
       * projection, creating a finer grid for the transformation and
       * rebinning to the same grid size again.
       */
      Skymap<T> CoordTransform(SM::CoordSys from, SM::CoordSys to, double mj){
         //If from and to are the same, just return this
         if (from == to)
            return *this;

         //Create the transformation function.  We do it reversed,
         //since we take the to coordinate and transfer to from
         TransformFunction f(to,from);

         //Call it once to initalize the static variables in astro
         {
            pointing fp(0,0), tp;
            f(mj, fp, tp);
         }

         //Create the output map and a finer healpix grid to do the
         //transformation on
         Skymap<T> out(*this);
         out = 0;
         const int dOrder = std::min(5,13-order_); //Order can not exceed 13
         const int nPix = (1<<dOrder)*(1<<dOrder);

         //Create a nested healpix base for us to work with
         const Healpix_Base finer(order_+dOrder, NEST);
         //Loop over the map, then the finer base and fill the map
#pragma omp parallel for default(shared) schedule(static)
         for (int p = 0; p<npix_; ++p){
            //Use pixel numbering from the nested scheme
            const int pp = (scheme_ == NEST) ? p : ring2nest(p);
            //Loop over all of the subpixels
            for (int sp = pp*nPix; sp<(pp+1)*nPix; ++sp){
               const pointing pnt = finer.pix2ang(sp);
               pointing old;
               f(mj, pnt, old);
               out[p] += (*this)[old];
            }
            out[p] /= double(nPix);
         }
         return out;
      }

      /**\brief Print the skymap
       *
       * \param os is the output stream to which to print the skymap
       *
       * First outputs the coordinate and then the spectra.  This routine is for debugging of
       * small skymaps.  The output gets very large very quickly.  The SkymapFitsio routines are
       * recommended.
       */
      void print (std::ostream & os){
         for (int i = 0; i < npix_; ++i){
            SM::Coordinate co(pix2ang(i));
            os << co;
            os << std::endl << "  "; //!< First output the coordinates
            for(size_t j = 0; j < fdimension.Spectra().size(); ++j){
	      os << fMap[i*fdimension.Spectra().size() + j] << " ";
            }
            os<<std::endl;
         }
      }

      /** \brief Returns the minimum value of the map */
      T min () const {
         return fMap.min();
      }
      /** \brief Returns the minimum value of the map at a given point in the spectra 
       *
       * \param iSpectra is the point in spectra to evaluate the minimum value.
       */
      T min(int iSpectra) const {
         T m(0);
         if (iSpectra >= 0 && iSpectra < int(fdimension.Spectra().size())){
            m = fMap[iSpectra];
            for (int i = 1; i < npix_; ++i){
	      m = std::min(fMap[i*fdimension.Spectra().size()+iSpectra], m);
            }
         }
         return m;
      }

      /** \brief Returns the maximum value of the map */
      T max () const {
         return fMap.max();
      }
      /** \brief Returns the maximum value of the map at a given point in the spectra 
       *
       * \param iSpectra is the point in spectra to evaluate the maximum value.
       */
      T max(int iSpectra) const {
         T m(0);
         if (iSpectra >= 0 && iSpectra < int(fdimension.Spectra().size())){
            m = fMap[iSpectra];
            for (int i = 1; i < npix_; ++i){
	      m = std::max(fMap[i*fdimension.Spectra().size()+iSpectra], m);
            }
         }
         return m;
      }

      /**\brief Sum of all the pixels for a point in spectra
       *
       * \param iSpectra is the point in spectra to evaluate the sum
       */
      T sum (int iSpectra) const {
         T s(0);
         //Check if iSpectra is within bounds
         if (iSpectra >= 0 && iSpectra < int(fdimension.Spectra().size())){
            for (int i = 0; i < npix_; ++i){
	      s += fMap[i*fdimension.Spectra().size()+iSpectra];
            }
         }
         return s;
      }

      /**\brief Sum all of the pixels in the map */
      T sum () const {
         return fMap.sum();
      }

      /**\brief Assigns a single value to the map */ 
      template <typename C>
         Skymap<T> & operator = (const C number){
            fMap = T(number);
            return(*this);
         }

      /**\brief Assignment for maps, basically return the same map */
      Skymap<T> & operator = (const Skymap<T> & oldMap){
         //Avoid self-assignment
         if (this != &oldMap){
            Resize(oldMap.fdimension);
            fMap = oldMap.fMap;
         }
         return(*this);
      }

      //! Add a number to the current skymap
      template <typename C>
         Skymap<T> & operator += (C number) {
            fMap += T(number);
            return (*this);
         }
      //! Add a number to a skymap
      template <typename C>
         Skymap<T> operator + (C number) const {
            Skymap<T> returnMap(*this);
            returnMap += number;
            return returnMap;
         }

      //! Withdraw a number from the current skymap
      template <typename C>
         Skymap<T> & operator -= (C number) {
            fMap -= T(number);
            return (*this);
         }
      //! Withdraw a number from a skymap
      template <typename C>
         Skymap<T> operator - (C number) const {
            Skymap<T> returnMap(*this);
            returnMap -= number;
            return returnMap;
         }

      //! Multiply the current skymap with a number
      template <typename C>
         Skymap<T> & operator *= (C number) {
            fMap *= T(number);
            return (*this);
         }
      //! Multiply a skymap with a number
      template <typename C>
         Skymap<T> operator * (C number) const {
            Skymap<T> returnMap(*this);
            returnMap *= number;
            return returnMap;
         }

      //! Divide the current skymap with a number
      template <typename C>
         Skymap<T> & operator /= (C number) {
            fMap /= T(number);
            return (*this);
         }
      //! Divide a skymap with a number
      template <typename C>
         Skymap<T> operator / (C number) const {
            Skymap<T> returnMap(*this);
            returnMap /= number;
            return returnMap;
         }

      /**\brief Sets the mask for mathematical operations with skymaps
       *
       * The ORDERBIT cannot be unset.  Skymaps have to be rebinned prior to calling the functions.
       *
       * If SPBINNEDBIT is not set we unset SPVALUEBIT
       * If SPSIZEBIT is not set we unset both SPVALUEBIT and SPBINNEDBIT.
       * Note that we can only multiply maps with same spectral size or maps with size of 1.
       *
       * We handle the Scheme automatically so it makes little sense to set the SCHEMEBIT.
       */
      static void SetMathMask(unsigned int mask) {
         mathMask = mask;

         mathMask |= SkymapDim::ORDERBIT;

         if ( (mathMask & SkymapDim::SPSIZEBIT) == 0 ) {
            mathMask &= ~SkymapDim::SPVALUEBIT;
            mathMask &= ~SkymapDim::SPBINNEDBIT;
         } else if ( (mathMask & SkymapDim::SPBINNEDBIT) == 0 ) {
            mathMask &= ~SkymapDim::SPVALUEBIT;
         }
      }

      //! Return the current mask
      static unsigned int MathMask() {
         return mathMask;
      }
            

      /**\brief Add a skymap to the current one 
       *
       * Only adds the skymap if the values at which the spectra is evaluated are the same and the skymaps are the
       * same size.
       */
      Skymap<T> & operator += (const Skymap<T> & otherMap){
         unsigned int mask = fdimension.compare(otherMap.fdimension);
         if ( mathMask != (mask&mathMask) ) {
            WARNING("Skymaps are not compatible, no action performed");
            return (*this);
         }

         //We can only add if spectral sizes are the same or otherMap has size of 1
         if ( (mask & SkymapDim::SPSIZEBIT) == 0 && otherMap.nSpectra() > 1 ) {
            WARNING("Skymaps are not compatible, no action performed");
            return (*this);
         }

         if ( (mask & SkymapDim::SCHEMEBIT) != 0 ) {

            if (otherMap.nSpectra() > 1) {
            
               fMap += otherMap.fMap;
            
            } else {
            
               for (long i = 0; i < npix_; ++i) {
                  const size_t ibase = i*nSpectra();
                  for (size_t j = ibase; j < ibase+fdimension.Spectra().size(); ++j) {
                     fMap[j] += otherMap.fMap[i];
                  }
               }
            
            }
         
         } else {

            //Need to map pixels from one to the other
            int (Skymap<T>::* pixelConvert)(int)const;
            if ( fdimension.GetScheme() == NEST ) {
               pixelConvert = &Skymap<T>::nest2ring;
            } else {
               pixelConvert = &Skymap<T>::ring2nest;
            }

            if (otherMap.nSpectra() > 1) {

#pragma omp parallel for default(shared) schedule(static)
               for (long i = 0; i < npix_; ++i) {
                  const long io = (this->*pixelConvert)(i);
                  const size_t ibase = i*nSpectra();
                  const size_t iobase = io*nSpectra();
                  for (size_t j = 0; j < fdimension.Spectra().size(); ++j) {
                     fMap[ibase+j] += otherMap.fMap[iobase+j];
                  }
               }

            } else {

#pragma omp parallel for default(shared) schedule(static)
               for (long i = 0; i < npix_; ++i) {
                  const long io = (this->*pixelConvert)(i);
                  const size_t ibase = i*nSpectra();
                  for (size_t j = ibase; j < ibase+fdimension.Spectra().size(); ++j) {
                     fMap[j] += otherMap.fMap[io];
                  }
               }

            }

         }

         return (*this);
      }

      /**\brief Add a skymap to another skymap
       *
       * Only adds the skymaps if the values at which the spectra is evaluated are the same and the skymaps are the
       * same size.
       */
      Skymap<T> operator + (const Skymap<T> & otherMap) const{
         Skymap<T> returnMap(*this);
         returnMap += otherMap;
         return (returnMap);
      }

      /**\brief Subtract a skymap from the current one
       *
       * Only subtracts the skymap if the values at which the spectra is evaluated are the same and the skymaps are the
       * same size.
       */
      Skymap<T> & operator -= (const Skymap<T> & otherMap){
         unsigned int mask = fdimension.compare(otherMap.fdimension);
         if ( mathMask != (mask&mathMask) ) {
            WARNING("Skymaps are not compatible, no action performed");
            return (*this);
         }

         //We can only add if spectral sizes are the same or otherMap has size of 1
         if ( (mask & SkymapDim::SPSIZEBIT) == 0 && otherMap.nSpectra() > 1 ) {
            WARNING("Skymaps are not compatible, no action performed");
            return (*this);
         }

         if ( (mask & SkymapDim::SCHEMEBIT) != 0 ) {

            if (otherMap.nSpectra() > 1) {
            
               fMap -= otherMap.fMap;
            
            } else {
            
               for (long i = 0; i < npix_; ++i) {
                  const size_t ibase = i*nSpectra();
                  for (size_t j = ibase; j < ibase+fdimension.Spectra().size(); ++j) {
                     fMap[j] -= otherMap.fMap[i];
                  }
               }
            
            }
         
         } else {

            //Need to map pixels from one to the other
            int (Skymap<T>::* pixelConvert)(int)const;
            if ( fdimension.GetScheme() == NEST ) {
               pixelConvert = &Skymap<T>::nest2ring;
            } else {
               pixelConvert = &Skymap<T>::ring2nest;
            }

            if (otherMap.nSpectra() > 1) {

#pragma omp parallel for default(shared) schedule(static)
               for (long i = 0; i < npix_; ++i) {
                  const long io = (this->*pixelConvert)(i);
                  const size_t ibase = i*nSpectra();
                  const size_t iobase = io*nSpectra();
                  for (size_t j = 0; j < fdimension.Spectra().size(); ++j) {
                     fMap[ibase+j] -= otherMap.fMap[iobase+j];
                  }
               }

            } else {

#pragma omp parallel for default(shared) schedule(static)
               for (long i = 0; i < npix_; ++i) {
                  const long io = (this->*pixelConvert)(i);
                  const size_t ibase = i*nSpectra();
                  for (size_t j = ibase; j < ibase+fdimension.Spectra().size(); ++j) {
                     fMap[j] -= otherMap.fMap[io];
                  }
               }

            }

         }

         return (*this);
      }

      /**\brief Subtract a skymap from another skymap
       *
       * Only subtracts the skymaps if the values at which the spectra is evaluated are the same and the skymaps are the
       * same size.
       */
      Skymap<T> operator - (const Skymap<T> & otherMap) const{
         Skymap<T> returnMap(*this);
         returnMap -= otherMap;
         return (returnMap);
      }

      /**\brief Multiply a skymap with the current one
       *
       * Only multiplies the skymap if the values at which the spectra is evaluated are the same and the skymaps are the
       * same size.
       */
      Skymap<T> & operator *= (const Skymap<T> & otherMap){
         unsigned int mask = fdimension.compare(otherMap.fdimension);
         if ( mathMask != (mask&mathMask) ) {
            WARNING("Skymaps are not compatible, no action performed");
            return (*this);
         }

         //We can only add if spectral sizes are the same or otherMap has size of 1
         if ( (mask & SkymapDim::SPSIZEBIT) == 0 && otherMap.nSpectra() > 1 ) {
            WARNING("Skymaps are not compatible, no action performed");
            return (*this);
         }

         if ( (mask & SkymapDim::SCHEMEBIT) != 0 ) {

            if (otherMap.nSpectra() > 1) {
            
               fMap *= otherMap.fMap;
            
            } else {
            
               for (long i = 0; i < npix_; ++i) {
                  const size_t ibase = i*nSpectra();
                  for (size_t j = ibase; j < ibase+fdimension.Spectra().size(); ++j) {
                     fMap[j] *= otherMap.fMap[i];
                  }
               }
            
            }
         
         } else {

            //Need to map pixels from one to the other
            int (Skymap<T>::* pixelConvert)(int)const;
            if ( fdimension.GetScheme() == NEST ) {
               pixelConvert = &Skymap<T>::nest2ring;
            } else {
               pixelConvert = &Skymap<T>::ring2nest;
            }

            if (otherMap.nSpectra() > 1) {

#pragma omp parallel for default(shared) schedule(static)
               for (long i = 0; i < npix_; ++i) {
                  const long io = (this->*pixelConvert)(i);
                  const size_t ibase = i*nSpectra();
                  const size_t iobase = io*nSpectra();
                  for (size_t j = 0; j < fdimension.Spectra().size(); ++j) {
                     fMap[ibase+j] *= otherMap.fMap[iobase+j];
                  }
               }

            } else {

#pragma omp parallel for default(shared) schedule(static)
               for (long i = 0; i < npix_; ++i) {
                  const long io = (this->*pixelConvert)(i);
                  const size_t ibase = i*nSpectra();
                  for (size_t j = ibase; j < ibase+fdimension.Spectra().size(); ++j) {
                     fMap[j] *= otherMap.fMap[io];
                  }
               }

            }

         }

         return (*this);
      }

      /**\brief Multiply a skymap with another skymap
       *
       * Only multiplies the skymaps if the values at which the spectra is evaluated are the same and the skymaps are the
       * same size.
       */
      Skymap<T> operator * (const Skymap<T> & otherMap) const{
         Skymap<T> returnMap(*this);
         returnMap *= otherMap;
         return (returnMap);
      }

      /**\brief Divide a skymap into the current one
       *
       * Only divides the skymap if the values at which the spectra is evaluated are the same and the skymaps are the
       * same size.
       */
      Skymap<T> & operator /= (const Skymap<T> & otherMap){
         unsigned int mask = fdimension.compare(otherMap.fdimension);
         if ( mathMask != (mask&mathMask) ) {
            WARNING("Skymaps are not compatible, no action performed");
            return (*this);
         }

         //We can only add if spectral sizes are the same or otherMap has size of 1
         if ( (mask & SkymapDim::SPSIZEBIT) == 0 && otherMap.nSpectra() > 1 ) {
            WARNING("Skymaps are not compatible, no action performed");
            return (*this);
         }

         if ( (mask & SkymapDim::SCHEMEBIT) != 0 ) {

            if (otherMap.nSpectra() > 1) {
            
               fMap /= otherMap.fMap;
            
            } else {
            
               for (long i = 0; i < npix_; ++i) {
                  const size_t ibase = i*nSpectra();
                  for (size_t j = ibase; j < ibase+fdimension.Spectra().size(); ++j) {
                     fMap[j] /= otherMap.fMap[i];
                  }
               }
            
            }
         
         } else {

            //Need to map pixels from one to the other
            int (Skymap<T>::* pixelConvert)(int)const;
            if ( fdimension.GetScheme() == NEST ) {
               pixelConvert = &Skymap<T>::nest2ring;
            } else {
               pixelConvert = &Skymap<T>::ring2nest;
            }

            if (otherMap.nSpectra() > 1) {

#pragma omp parallel for default(shared) schedule(static)
               for (long i = 0; i < npix_; ++i) {
                  const long io = (this->*pixelConvert)(i);
                  const size_t ibase = i*nSpectra();
                  const size_t iobase = io*nSpectra();
                  for (size_t j = 0; j < fdimension.Spectra().size(); ++j) {
                     fMap[ibase+j] /= otherMap.fMap[iobase+j];
                  }
               }

            } else {

#pragma omp parallel for default(shared) schedule(static)
               for (long i = 0; i < npix_; ++i) {
                  const long io = (this->*pixelConvert)(i);
                  const size_t ibase = i*nSpectra();
                  for (size_t j = ibase; j < ibase+fdimension.Spectra().size(); ++j) {
                     fMap[j] /= otherMap.fMap[io];
                  }
               }

            }

         }

         return (*this);
      }

      /**\brief Divide a skymap with another skymap
       *
       * Only divides the skymaps if the values at which the spectra is evaluated are the same and the skymaps are the
       * same size.
       */
      Skymap<T> operator / (const Skymap<T> & otherMap) const{
         Skymap<T> returnMap(*this);
         returnMap /= otherMap;
         return (returnMap);
      }


      /**\brief Equality operator, returns true if every pixel has the same value and the values at which the spectra
       * is evaluated is the same as well.
       */
      bool operator == (const Skymap<T> & otherMap) const{
         //Check for equivalence of the map
         bool map = fdimension == otherMap.fdimension;
         if (! map) return map;
         //And the map, if needed
         for (size_t i = 0; i < fMap.size(); ++i){
            map &= (fMap[i] == otherMap.fMap[i]);
            if (! map) {
               return map;
            }
         }
         return map;
      }

      /** \brief Non-equality operator, the inverse of the equality operator. */
      bool operator != (const Skymap<T> & otherMap) const{
         return ( ! ((*this) == otherMap));
      }

};

//Ignore the scheme bit and handle it by conversion
template <typename T>
unsigned int Skymap<T>::mathMask = SkymapDim::ORDERBIT | SkymapDim::SPBINNEDBIT | SkymapDim::SPSIZEBIT | SkymapDim::SPVALUEBIT;

#endif
