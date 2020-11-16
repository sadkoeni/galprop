#include <BaseSkyFitsIO.h>
#include <FullSky.h>

using namespace SM;

int main() {

   //Create a map of order 6 and fill it with f(E,b,l) = cos(b)*cos(0.5*l) * E**-2.4
   FullSky<double> map(SpectralBinning(std::vector<double>({0.1, 0.3, 1.0, 3.0, 10.0})), 4, RING, CoordSys::GAL, 0.0);

   for ( auto it = map.begin(); it != map.end(); ++it) {

      //Get spatial coordinates
      Coordinate co(map.pix2ang(it.healpixIndex()), map.GetCoordinateSystem());
      double l, b;
      co.getCoordinates(l, b, map.GetCoordinateSystem());

      //Get energy
      const double E = map.GetBinning().GetBins()[it.energyIndex()];

      *it = cos(b)*fabs(cos(0.5*l)) * pow(E,-2.4);

   }

   //Write out the original map
   writeToFits(map, "testIAC_original.fits");

   //Do coordinate conversion to ecliptic
   auto eclMapPtr = map.CoordinateConversion(CoordSys::ECL);

   writeToFits(*eclMapPtr, "testIAC_ECLconv.fits");

   delete eclMapPtr.release();

   //Now do similar, but with interpolation.  Also do spectral interpolation and extrapolation
   FullSky<double> eclMap(SpectralBinning(std::vector<double>({0.01, 0.2, 2.0, 8.0, 30.0})), 4, RING, CoordSys::ECL, 0.0);

   map.Interpolate(eclMap, true);

   writeToFits(eclMap, "testIAC_ECLextrap.fits");

   //Now without extrapolation
   map.Interpolate(eclMap, false);
   writeToFits(eclMap, "testIAC_ECLinterp.fits");

   //Make sure we can do this with finer binning
   FullSky<double> eclMap2(SpectralBinning(std::vector<double>({0.01, 0.2, 2.0, 8.0, 30.0})), 8, RING, CoordSys::ECL, 0.0);
   map.Interpolate(eclMap2, false);
   writeToFits(eclMap2, "testIAC_ECLinterp_finer.fits");

   //Make sure we can do this with coarser binning
   FullSky<double> eclMap3(SpectralBinning(std::vector<double>({0.01, 0.2, 2.0, 8.0, 30.0})), 1, RING, CoordSys::ECL, 0.0);
   map.Interpolate(eclMap3, false);
   writeToFits(eclMap3, "testIAC_ECLinterp_coarser.fits");

   //Now lets rebin the finer map and see how that goes
   auto eclMap2ptr = eclMap2.ResetOrder(4, false);

   writeToFits(*eclMap2ptr, "testIAC_ECLinterp_finerRebin.fits");

   //And finally rebin the coarser map just for fun
   auto eclMap3ptr = eclMap3.ResetOrder(4, false);
   writeToFits(*eclMap3ptr, "testIAC_ECLinterp_coarserRebin.fits");

   return 0;
}
