#include "Skymap.h"

void empty(double mj, double x, double y, double *z, double *w) {
   *z = x;
   *w = y;
}

SkymapDim::SkymapDim (const SkymapDim & other) {
  /*   binned = other.binned;
   healpixOrder = other.healpixOrder;
   healpixScheme = other.healpixScheme;

   spectra.resize(other.spectra.size());
   spectra = other.spectra;

   spectraMin.resize(other.spectraMin.size());
   spectraMin = other.spectraMin;

   spectraMax.resize(other.spectraMax.size());
   spectraMax = other.spectraMax;*/

  *this = other;

}

SkymapDim & SkymapDim::operator = (const SkymapDim & other) {
   if (this != &other) {
      binned = other.binned;
      healpixOrder = other.healpixOrder;
      healpixScheme = other.healpixScheme;

      spectra.resize(other.spectra.size());
      spectra = other.spectra;

      spectraMin.resize(other.spectraMin.size());
      spectraMin = other.spectraMin;

      spectraMax.resize(other.spectraMax.size());
      spectraMax = other.spectraMax;
   }
   return *this;
}

bool SkymapDim::operator == (const SkymapDim & other) const {
   //Make sure all the bits are set
   unsigned int mask = ORDERBIT | SCHEMEBIT | SPSIZEBIT | SPBINNEDBIT | SPVALUEBIT;
   return ( mask == (compare(other)&mask) );
   /*
   bool eq = (binned == other.binned);
   eq &= (spectra.size() == other.spectra.size());
   eq &= healpixOrder == other.healpixOrder;
   eq &= healpixScheme == other.healpixScheme;
   if (eq) {
      //Check for equivalence of spectra
      if (binned) {
         for (size_t i = 0; i < spectra.size(); ++i){
            eq &= (fabs((spectraMin[i] - other.spectraMin[i])/spectraMin[i]) < 1e-6);
            eq &= (fabs((spectraMax[i] - other.spectraMax[i])/spectraMax[i]) < 1e-6);
         }
      }else{
         for (size_t i = 0; i < spectra.size(); ++i){
            eq &= (fabs((spectra[i] - other.spectra[i])/spectra[i]) < 1e-6);
         }
      }
   }
   return (eq);
   */
}

bool SkymapDim::operator != (const SkymapDim & other) const{
   return (! ( (*this) == other ) );
}

unsigned int SkymapDim::compare (const SkymapDim & other) const {

   unsigned int mask(0);

   if (healpixOrder == other.healpixOrder)
      mask |= ORDERBIT;

   if (healpixScheme == other.healpixScheme)
      mask |= SCHEMEBIT;

   if (binned == other.binned)
      mask |= SPBINNEDBIT;

   if (spectra.size() == other.spectra.size())
      mask |= SPSIZEBIT;

   if ( (mask & SPBINNEDBIT) != 0 && (mask & SPSIZEBIT) != 0 ) {
      bool eq(true);
      if (binned) {
         for (size_t i = 0; i < spectra.size(); ++i){
            eq &= (fabs((spectraMin[i] - other.spectraMin[i])/spectraMin[i]) < 1e-6);
            eq &= (fabs((spectraMax[i] - other.spectraMax[i])/spectraMax[i]) < 1e-6);
         }
      }else{
         for (size_t i = 0; i < spectra.size(); ++i){
            eq &= (fabs((spectra[i] - other.spectra[i])/spectra[i]) < 1e-6);
         }
      }
      if (eq) 
         mask |= SPVALUEBIT;
   }

   return mask;

}

