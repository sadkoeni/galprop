#include "BaseSky.h"
#include "FullSky.h"
#include "SparseSky.h"

using namespace SM;

SpectralBinning::SpectralBinning(std::vector<double>&& sp) :
   spectra(sp)
{
   AssertSpectra();
}

SpectralBinning::SpectralBinning(const std::vector<double>& sp) :
   spectra(sp)
{
   AssertSpectra();
}

SpectralBinning::SpectralBinning(const std::vector<double>& spMin, const std::vector<double>& spMax) :
   spectraMin(spMin),
   spectraMax(spMax)
{
   AssertSpectra();
}

SpectralBinning::SpectralBinning(std::vector<double>&& spMin, std::vector<double>&& spMax) :
   spectraMin(spMin),
   spectraMax(spMax)
{
   AssertSpectra();
}

void SpectralBinning::SetBinning(const std::vector<double>& sp)
{
   if (sp.size() != spectra.size())
      throw (std::invalid_argument("Cannot change size in a call to SetBinning in SpectralBinning"));

   spectra = sp;

   spectraMin.resize(0);
   spectraMax.resize(0);

   AssertSpectra();
}

void SpectralBinning::SetBinning(const std::vector<double>& spMin, const std::vector<double>& spMax)
{
   if (spMin.size() != spectraMin.size() || spMax.size() != spectraMax.size())
      throw (std::invalid_argument("Cannot change size in a call to SetBinning in SpectralBinning"));

   spectraMin = spMin;
   spectraMax = spMax;

   spectra.resize(0);

   AssertSpectra();
}

size_t SpectralBinning::GetIndex(double value) const noexcept
{


   if ( IsIntegrated() ) {

      auto it = std::lower_bound(spectraMax.begin(), spectraMax.end()-1, value);

      return static_cast<size_t>(it - spectraMax.begin());

   } else {

      std::vector<double> diff;
      diff.reserve(spectra.size());

      std::transform(spectra.begin(), spectra.end(), std::back_inserter(diff), [&value](double v){return fabs(v - value);});

      auto it = std::min_element(diff.begin(), diff.end());

      return static_cast<size_t>(it - diff.begin());

   }

}

bool SpectralBinning::IsIdentical( const SpectralBinning &other) const noexcept
{
   //Use flexible relative value comparison
   static const double error(1e-6);
   bool tmp = spectra.size() == other.spectra.size();
   tmp &= spectraMin.size() == other.spectraMin.size();
   if (tmp) {
      if ( spectra.size() ) {
         for (size_t i=0; i < spectra.size(); ++i) {
            tmp &= ( fabs(spectra[i] - other.spectra[i]) < fabs(spectra[i]) * error );
         }
      } else {
         for (size_t i=0; i < spectraMin.size(); ++i) {
            tmp &= ( fabs(spectraMin[i] - other.spectraMin[i]) < fabs(spectraMin[i]) * error );
            tmp &= ( fabs(spectraMax[i] - other.spectraMax[i]) < fabs(spectraMax[i]) * error );
         }
      }
   }
   return tmp;
}

void SpectralBinning::AssertSpectra() 
{

   if (! std::is_sorted(spectra.begin(), spectra.end()) )
      throw (std::invalid_argument("Spectrum needs to be sorted in SpectralBinning"));

   if (spectraMin.size() != spectraMax.size())
      throw (std::invalid_argument("Minimum and maximum of spectral bins not same size in SpectralBinning"));

   if (! std::is_sorted(spectraMin.begin(), spectraMin.end()) || ! std::is_sorted(spectraMax.begin(), spectraMax.end()))
      throw (std::invalid_argument("Spectrum needs to be sorted in SpectralBinning"));

   for (size_t i = 0; i+1 < spectraMin.size(); ++i) 
      if (fabs(spectraMax[i]-spectraMin[i+1]) > 1e-8*spectraMax[i])
         throw (std::invalid_argument("Spectral bins must be adjacent in SpectralBinning"));

   if (spectra.size() == 0 && spectraMin.size() == 0)
      throw (std::invalid_argument("Need at least one value for the spectrum"));

}

REGISTER_SKY_CLASS(FullSky, FullSky);
REGISTER_SKY_CLASS(SparseSky, SparseSky);

//Instantiate all the create functions for the registered templates
template std::unique_ptr< BaseSky<char> > BaseSky<char>::create(const std::string &, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, const char& );
template std::unique_ptr< BaseSky<short> > BaseSky<short>::create(const std::string &, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, const short& );
template std::unique_ptr< BaseSky<int> > BaseSky<int>::create(const std::string &, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, const int& );
template std::unique_ptr< BaseSky<long> > BaseSky<long>::create(const std::string &, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, const long& );
template std::unique_ptr< BaseSky<float> > BaseSky<float>::create(const std::string &, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, const float& );
template std::unique_ptr< BaseSky<double> > BaseSky<double>::create(const std::string &, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, const double& );
template std::unique_ptr< BaseSky<unsigned char> > BaseSky<unsigned char>::create(const std::string &, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, const unsigned char& );
template std::unique_ptr< BaseSky<unsigned short> > BaseSky<unsigned short>::create(const std::string &, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, const unsigned short& );
template std::unique_ptr< BaseSky<unsigned int> > BaseSky<unsigned int>::create(const std::string &, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, const unsigned int& );
template std::unique_ptr< BaseSky<unsigned long> > BaseSky<unsigned long>::create(const std::string &, const SpectralBinning &, int, Healpix_Ordering_Scheme, CoordSys, const unsigned long& );

