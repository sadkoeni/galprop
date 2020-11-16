#include <FullSky.h>
#include <BaseSkyFitsIO.h>

int main(int argc, char* argv[]) {

   if (argc != 3) {
      std::cerr<<"Usage: "<<argv[0]<<" <inFile> <outFile>"<<std::endl;
      return 1;
   }

   //Create a map to read into
   auto map = SM::BaseSky<float>::create("FullSky", SM::SpectralBinning(std::vector<double>({0.0})), 9, RING, SM::CoordSys::GAL, 0.0);
   //std::unique_ptr< SM::BaseSky<float> > map;

   //Read in using wcs
   SM::readFromFitsWcs(map, argv[1], false);

   //Write out
   SM::writeToFits(*map, argv[2]); 

   return 0;
}
