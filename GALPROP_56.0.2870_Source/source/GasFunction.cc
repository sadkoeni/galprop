#include "Galprop.h"
#include "galprop_internal.h"
#include <stdexcept>
#include <cmath>

Galprop::GasFunction::GasFunction( const std::string &type, double rPLindex, const Galprop& gp) :
   fgp(gp),
   frInd(rPLindex)
{
   if ( type == "HI" )
      ftype = HI;
   else if ( type == "H2" )
      ftype = H2;
   else if ( type == "CO" )
      ftype = CO;
   else if ( type == "HII" )
      ftype = HII;
   else
      throw(std::invalid_argument("Type can only be \"HI\", \"H2\", \"CO\", or \"HII\" in GasFunction"));
}

double Galprop::GasFunction::operator () ( const double x, const double y, const double z, const vec3 &dir ) const{
   const double r = sqrt(x*x + y*y);
   const double rScale = pow(r, frInd);
   switch (ftype) {
      case HI:
	 return nHI3D( x, y, z )*rScale;
      case CO:
	 return 2*nH23D( x, y, z, 1.0 )*rScale;
      case H2:
	 return 2*nH23D( x, y, z, fgp.fX_CO(r) )*rScale;
      case HII:
	 return nHII3D( x, y, z )*rScale;
      default:
	 return 0;
   }
}
