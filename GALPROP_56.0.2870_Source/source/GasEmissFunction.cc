#include "Galprop.h"
#include "DistributionFunction.h"
#include <valarray>
#include <stdexcept>

Galprop::GasEmissFunction::GasEmissFunction( const std::string &type, const std::string &gas_type, const Galprop&gp ) :
   fgp(gp),
   fgf(gas_type, 0.0, gp),
   fdf(0)
{
   //Set up arrays, depending on dimension of galprop
   std::vector<double> z(&fgp.galaxy.z[0], &fgp.galaxy.z[0]+fgp.galaxy.n_zgrid);

   if ( type == "BREMSS" || type == "Bremss" || type == "bremss" ) {
      ftype = BREMSS;

      if ( fgp.galaxy.n_spatial_dimensions == 2 ) {
	 std::vector<double> r(&fgp.galaxy.r[0], &fgp.galaxy.r[0]+fgp.galaxy.n_rgrid);

         if ( gas_type == "HII" ) {
            fdf = new DistributionFunction(fgp.galaxy.bremss_ionized_emiss, z, r);
         } else {
            fdf = new DistributionFunction(fgp.galaxy.bremss_emiss, z, r);
         }

      } else if ( fgp.galaxy.n_spatial_dimensions == 3 ) {
	 std::vector<double> x(&fgp.galaxy.x[0], &fgp.galaxy.x[0]+fgp.galaxy.n_xgrid);
	 std::vector<double> y(&fgp.galaxy.y[0], &fgp.galaxy.y[0]+fgp.galaxy.n_ygrid);

         if ( gas_type == "HII" ) {
            fdf = new DistributionFunction(fgp.galaxy.bremss_ionized_emiss, z, x, y);
         } else {
            fdf = new DistributionFunction(fgp.galaxy.bremss_emiss, z, x, y);
         }
      }

   } else if ( type == "PION" || type == "Pion" || type == "pion" || type == "PI0" || type == "Pi0" || type == "pi0" ) {
      ftype = PION;
      if ( fgp.galaxy.n_spatial_dimensions == 2 ) {
	 std::vector<double> r(&fgp.galaxy.r[0], &fgp.galaxy.r[0]+fgp.galaxy.n_rgrid);
	 fdf = new DistributionFunction(fgp.galaxy.pi0_decay_emiss, z, r);
      } else if ( fgp.galaxy.n_spatial_dimensions == 3 ) {
	 std::vector<double> x(&fgp.galaxy.x[0], &fgp.galaxy.x[0]+fgp.galaxy.n_xgrid);
	 std::vector<double> y(&fgp.galaxy.y[0], &fgp.galaxy.y[0]+fgp.galaxy.n_ygrid);
	 fdf = new DistributionFunction(fgp.galaxy.pi0_decay_emiss, z, x, y);
      }
   } else if ( type == "TEST" || type == "Test" || type == "test" ) {
      ftype = TEST;
      test = fgp.galaxy.pi0_decay_emiss;
      test = 1;
      if ( fgp.galaxy.n_spatial_dimensions == 2 ) {
	 std::vector<double> r(&fgp.galaxy.r[0], &fgp.galaxy.r[0]+fgp.galaxy.n_rgrid);
	 fdf = new DistributionFunction(test, z, r);
      } else if ( fgp.galaxy.n_spatial_dimensions == 3 ) {
	 std::vector<double> x(&fgp.galaxy.x[0], &fgp.galaxy.x[0]+fgp.galaxy.n_xgrid);
	 std::vector<double> y(&fgp.galaxy.y[0], &fgp.galaxy.y[0]+fgp.galaxy.n_ygrid);
	 fdf = new DistributionFunction(test, z, x, y);
      }
   } else {
      throw(std::invalid_argument("Type can only be \"PION\", \"BREMSS\", or \"TEST\" in GasEmissFunction"));
   }
}

Galprop::GasEmissFunction::~GasEmissFunction() {
   delete fdf;
}

std::valarray<double> Galprop::GasEmissFunction::operator () ( const double x, const double y, const double z, const vec3 &dir ) const {
   //Just multiply the gas function with the distribution function, easy peasy
   return (*fdf)(x,y,z,dir)*fgf(x,y,z,dir);
}
