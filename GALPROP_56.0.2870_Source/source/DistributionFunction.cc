#include "DistributionFunction.h"
#include "Distribution.h"

#include <valarray>
#include <cmath>
#include <stdexcept>

DistributionFunction::DistributionFunction(const Distribution &dist, const std::vector<double> &z, const std::vector<double> &r) : 
   fdist(dist),
   fz(z),
   fr(r)
{
   //Check for dimensions in distribution
   if ( dist.n_spatial_dimensions != 2 )
      throw(std::invalid_argument ( "Distribution must have 2 spatial dimensions when DistributionFunction is initialized with radial boundaries" ));
   if ( dist.n_rgrid != r.size() )
      throw(std::invalid_argument ( "DistributionFunction: Size of r array must equal size of r grid in Distribution" ));
   if ( dist.n_zgrid != z.size() )
      throw(std::invalid_argument ( "DistributionFunction: Size of z array must equal size of z grid in Distribution" ));
}

DistributionFunction::DistributionFunction(const Distribution &dist, const std::vector<double> &z, const std::vector<double> &x, const std::vector<double> &y) : 
   fdist(dist),
   fz(z),
   fx(x),
   fy(y)
{
   //Check for dimensions in distribution
   if ( dist.n_spatial_dimensions != 3 )
      throw(std::invalid_argument ( "Distribution must have 3 spatial dimensions when DistributionFunction is initialized with Cartesian boundaries" ));
   if ( dist.n_xgrid != x.size() )
      throw(std::invalid_argument ( "DistributionFunction: Size of x array must equal size of x grid in Distribution" ));
   if ( dist.n_ygrid != y.size() )
      throw(std::invalid_argument ( "DistributionFunction: Size of y array must equal size of y grid in Distribution" ));
   if ( dist.n_zgrid != z.size() )
      throw(std::invalid_argument ( "DistributionFunction: Size of z array must equal size of z grid in Distribution" ));
}


std::valarray<double> DistributionFunction::operator () ( const double x, const double y, const double z, const vec3 &dir ) const {

   std::valarray<double> output(fdist.n_pgrid);

/*
   //We always interpolate in z, 
   size_t izu(fz.size()-1);
   size_t izl(0);
   while (izu-izl > 1) {
      size_t iz = (izu+izl)/2;
      if ( z < fz[iz] ) {
	 izu = iz;
      } else {
	 izl = iz;
      }
   }
*/
   //Assume the grid is split into linearly even bins
   const double dz = (z - fz[0])/(fz[fz.size()-1]-fz[0]);
   const size_t izl = dz <= 0. ? 0 : dz < 1. ? dz*(fz.size()-1) : fz.size()-2;
   const size_t izu = izl + 1;

   const double rz = (fz[izu]-z)/(fz[izu]-fz[izl]);
   const double lzf = rz < 0. ? 0. : rz > 1. ? 1. : rz;
   const double uzf = 1-lzf;

   //Check for the dimensions of the distribution, to find how to interpolate
   if ( fdist.n_spatial_dimensions == 2 ) {
      const double r = sqrt(x*x+y*y);
/*
      size_t iru(fr.size()-1);
      size_t irl(0);
      while (iru-irl > 1){
	 size_t ir = (iru+irl)/2;
	 if ( r < fr[ir] ) {
	    iru = ir;
	 } else {
	    irl = ir;
	 }
      }
*/
      const double dr = (r - fr[0]) / (fr[fr.size()-1] - fr[0]);
      const size_t irl = dr <= 0. ? 0 : dr < 1. ? dr*(fr.size()-1) : fr.size()-2;
      const size_t iru = irl + 1;

      const double rr = (fr[iru]-r)/(fr[iru]-fr[irl]);
      const double lrf = rr < 0. ? 0. : rr > 1. ? 1. : rr;
      const double urf = 1-lrf;

      for ( size_t i = 0; i < fdist.n_pgrid; ++i ) 
	 output[i] = lzf*lrf*fdist.d2[irl][izl].s[i]
	           + lzf*urf*fdist.d2[iru][izl].s[i]
		   + uzf*lrf*fdist.d2[irl][izu].s[i]
		   + uzf*urf*fdist.d2[iru][izu].s[i];

   } else if ( fdist.n_spatial_dimensions == 3 ) {
/*
      size_t ixu(fx.size()-1);
      size_t ixl(0);
      while (ixu-ixl > 1){
	 size_t ix = (ixu+ixl)/2;
	 if ( x < fx[ix] ) {
	    ixu = ix;
	 } else {
	    ixl = ix;
	 }
      }
*/
      const double dx = (x - fx[0]) / (fx[fx.size()-1] - fx[0]);
      const size_t ixl = dx <= 0. ? 0 : dx < 1. ? dx*(fx.size()-1) : fx.size()-2;
      const size_t ixu = ixl + 1;

      const double rx = (fx[ixu]-x)/(fx[ixu]-fx[ixl]);
      const double lxf = rx < 0. ? 0. : rx > 1. ? 1. : rx;
      const double uxf = 1-lxf;

/*
      size_t iyu(fy.size()-1);
      size_t iyl(0);
      while (iyu-iyl > 1){
	 size_t iy = (iyu+iyl)/2;
	 if ( y < fy[iy] ) {
	    iyu = iy;
	 } else {
	    iyl = iy;
	 }
      }
*/
      const double dy = (y - fy[0]) / (fy[fy.size()-1] - fy[0]);
      const size_t iyl = dy <= 0. ? 0 : dy < 1. ? dy*(fy.size()-1) : fy.size()-2;
      const size_t iyu = iyl + 1;

      const double ry = (fy[iyu]-y)/(fy[iyu]-fy[iyl]);
      const double lyf = ry < 0. ? 0. : ry > 1. ? 1. : ry;
      const double uyf = 1-lyf;
      for ( size_t i = 0; i < fdist.n_pgrid; ++i ) 
	 output[i] = lxf*lyf*lzf*fdist.d3[ixl][iyl][izl].s[i]
	           + lxf*lyf*uzf*fdist.d3[ixl][iyl][izu].s[i]
		   + lxf*uyf*lzf*fdist.d3[ixl][iyu][izl].s[i]
		   + lxf*uyf*uzf*fdist.d3[ixl][iyu][izu].s[i]
	           + uxf*lyf*lzf*fdist.d3[ixu][iyl][izl].s[i]
	           + uxf*lyf*uzf*fdist.d3[ixu][iyl][izu].s[i]
		   + uxf*uyf*lzf*fdist.d3[ixu][iyu][izl].s[i]
		   + uxf*uyf*uzf*fdist.d3[ixu][iyu][izu].s[i];
   }

   return output;
}
