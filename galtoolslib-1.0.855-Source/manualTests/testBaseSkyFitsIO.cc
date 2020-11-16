#include <FullSky.h>
#include <BaseSkyFitsIO.h>

int main() {

   //Create a small map and fill it in with values
   SM::FullSky<int> map(SM::SpectralBinning(std::vector<double>({1.0, 2.0, 6.0})), 3, RING, SM::CoordSys::GAL, 0.0); 

   for (size_t i = 0; i < size_t(map.Npix()); ++i)
      for (size_t j = 0; j < map.GetBinning().GetSize(); ++j)
         map.SetValue(i, j, i*(j+1));


   SM::writeToFits(map, "thisIsIntegerFullSky.fits");

   //Now read it in to a FullSky skymap.
   std::unique_ptr< SM::BaseSky<int> > ptr2(new SM::FullSky<int> (SM::SpectralBinning(std::vector<double>({1.0, 2.0})), 0, RING, SM::CoordSys::GAL, 0.0));

   std::map<std::string, std::string> keywords;
   SM::readFromFits(ptr2, "thisIsIntegerFullSky.fits", keywords);

   if (! map.Equivalent(*ptr2)) {
      std::cerr<<"Maps are not equivalent after reading map"<<std::endl;
      return 1;
   }

   int different = SM::ApplyFunctionAccumulate<int,int>(map, *ptr2, [](int& a, const int& b) -> int { return a != b; });

   if (different > 0) {
      std::cerr<<"Something is wrong when reading map "<<different<<std::endl;
      return 1;
   }
   
   //Now try with an empty pointer
   std::unique_ptr< SM::BaseSky<int> > ptr3;
   SM::readFromFits(ptr3, "thisIsIntegerFullSky.fits", keywords);

   if (ptr3->name() != "FullSky") {
      std::cerr<<"Incorrect type returned "<<ptr3->name()<<std::endl;
      return 1;
   }

   if (! map.Equivalent(*ptr3)) {
      std::cerr<<"Maps are not equivalent after reading map"<<std::endl;
      return 1;
   }

   different = SM::ApplyFunctionAccumulate<int,int>(map, *ptr3, [](int& a, const int& b) -> int { return a != b; });

   if (different > 0) {
      std::cerr<<"Something is wrong when reading map "<<different<<std::endl;
      return 1;
   }
   
   return 0;
}
