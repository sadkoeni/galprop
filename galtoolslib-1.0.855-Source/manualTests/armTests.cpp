#include <spiralarms.h>
#include <Registry.h>

namespace gs = GalacticStructure;

int main() {

   //Dumps a whole bunch of output comparing the brute force method and approximate method for various values of a
   const double as[] = {3, 4.5, 5.5, 7.8, 10, 14};
   const size_t nas = sizeof(as)/sizeof(double);

   //Some radius values to test
   const double rs[] = {0.1, 1.4, 2.5, 6.3, 8.5, 10.4, 14.3, 29};
   const size_t nrs = sizeof(rs)/sizeof(double);

   //Evenly distributed theta values
   const size_t nth = 100;

   for (size_t ia(0); ia < nas; ++ia) {

      gs::ArmFunction arm( as[ia], 0.1, 0.23, 100, 1.0, utl::Registry0<gs::ScaleFunction>::create("Gaussian") );

      std::cout<<" # Calculating for a="<<as[ia]<<std::endl;

      for ( size_t ir(0); ir < nrs; ++ir) {

         for (size_t it(0); it < nth; ++it) {

            const double th = 2.*it*M_PI/(nth-1);

            double acc, app;

            arm.CompareArmDist( rs[ir], th, acc, app );

            std::cout<<rs[ir]<<" "<<th<<" "<<acc<<" "<<app<<" "<<(app-acc)/acc<<" "<<app/acc<<std::endl;

         }

         std::cout<<std::endl<<std::endl;
      }

   }

   return 0;
}

