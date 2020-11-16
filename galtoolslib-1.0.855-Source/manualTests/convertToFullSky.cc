#include <FullSky.h>
#include <BaseSkyFitsIO.h>

int main(int argc, char* argv[]) {

   if (argc != 2) {
      std::cerr<<"Usage: "<<argv[0]<<" <File>"<<std::endl;
      return 1;
   }

   //Create a map to read into
   std::unique_ptr< SM::BaseSky<double> > map(new SM::FullSky<double>(SM::SpectralBinning(std::vector<double>({0.0})), 7, RING, SM::CoordSys::GAL, 0.0));

   //Read in using wcs
   std::map<std::string, std::string> keywords;
   SM::readFromFits(map, argv[1], keywords);

   //Write out
   SM::writeToFits(*map, argv[1]); 

   return 0;
}
