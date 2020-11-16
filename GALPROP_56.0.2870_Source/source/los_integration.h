/*
 * Generic line of sight integration routine.
 * Provided a vector of functions (x,y,z), it returns an array (R) where the
 * output is binned into Galacto-centric rings.
 */

#ifndef _los_integration_h_
#define _los_integration_h_

#include <vec3.h>
#include "integ.h"
#include "constants.h"
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <memory>
#include "los_intersection.h"
#include "integ.h"
using namespace std;

namespace SM {
//! The function to be integrated over the LOS
template <typename T>
class LOSfunction {
public:
   virtual ~LOSfunction() {}
   //! x,y,z are the spatial coordinates for evaluation and dir is the direction the camera is pointing
   virtual T operator () ( const double x, const double y, const double z, const vec3 &dir ) const = 0;
};

//! Out of bounds error
class LOSboundError : public out_of_range {
public:
    explicit LOSboundError ( const string& arg ) :
        out_of_range ( "Out of bounds error in LOS integration: " + arg )
   {}
};

//! Invalid argument error
class LOSargumentError : public invalid_argument {
public:
   explicit LOSargumentError ( const string& arg ) :
      invalid_argument ( "Invalid argument in LOS integration: " + arg )
   {}
};

/** \brief The LOS integrator
 *
 * Assumes boundaries are in units of kpc but integrates in units of cm
 * The plane b = 0 is parallel to z = 0, b = 90 points in positive z, l = 0 point in negative x, and l = 90 points in negative y from
 * camera location. 
 *
 * If size of Rbins is < 2 then no binning in R is done.
 * If Rbins are given and the bins don't encompass the entire integration region, the rest is added as a single additional bin.
 * 
 */
template <typename T>
class LOSintegrator {
    double fxmin, fxmax, fymin, fymax, fzmin, fzmax, fRmax;
    std::vector<double> fcamera;
    std::vector<double> fRbins;
    const double fstepSize;
    const double fminStepSize;
    const double frelTolerance;
    const double fabsTolerance;
    const unsigned int fdimensions;
    size_t fnRbins; //Number of bins, different depending on binsEncompass
    size_t faddBoundaries; //0 if rings do not encompass 0 but do encompass the entire region.

    size_t findAnnulus(double R) const;
    void checkBins();
    void output(T v) const;

    static bool lessLimit(const std::pair< std::pair<double,double>, int> &a, const std::pair< std::pair<double,double>, int> &b);

    // Class for integration with the simpsons method in integ.h
    class simFun : public integ::vecFun<T> {
          //Camera location to calculate coordinates
          const std::vector<double> &fcamera;
          // The functions to integrate
          const std::vector< std::unique_ptr< LOSfunction<T> > > &ffuncs;
          const std::unique_ptr< LOSfunction<T> > &fabsorption;
          // Need temporary storage for tau
          mutable std::vector<T> tau;
          mutable std::vector<double> stau;
          mutable unsigned int itau;
          // The los direction
          const double fcosl, fcosb, fsinl, fsinb;
          const vec3 fdir;
          // Given a distance s along the line of sight, returns a x,y,z
          void coordinates(const double s, double &x, double &y, double &z) const;

       public:
          simFun(const double cosl, const double cosb, const double sinl, const double sinb, 
                const std::vector<double> & camera, 
                const std::vector< std::unique_ptr< LOSfunction<T> > > &funcs, 
                const std::unique_ptr< LOSfunction<T> > &absorption) :
             fcamera(camera),
             ffuncs(funcs),
             fabsorption(absorption), 
             fcosl(cosl),
             fcosb(cosb),
             fsinl(sinl),
             fsinb(sinb),
             fdir(-cosb*cosl, -cosb*sinl, sinb) 
       {
             
             //Initialize the cache for tau, only need to store 5 values for tau
             itau = 0;
             if (absorption.get() != nullptr ) {
               tau.resize(5, (*fabsorption)(0,0,0,vec3(0,0,0)));
               stau.resize(5, 0.0);
               itau = 1;
             }
          }
          std::vector<T> operator () ( double s ) const;
    };

public:
    LOSintegrator ( const double Rmax, const double zmin, const double zmax, const std::vector<double> &Rbins, const std::vector<double> &camera, const double stepSize, const double relTolerance , const double minStepSize) :
            fzmin ( zmin ),
            fzmax ( zmax ),
            fRmax ( Rmax ),
            fcamera ( camera ),
            fRbins ( Rbins ),
            fstepSize ( stepSize ),
            fminStepSize ( minStepSize ),
            frelTolerance ( relTolerance ),
            fabsTolerance ( 0 ),
            fdimensions ( 2 ) {
        if ( fzmin >= fzmax )
            throw ( LOSboundError ( "zmin must be > zmax" ) );
        if ( fRmax <= 0 )
            throw ( LOSboundError ( "Rmax must be > 0" ) );
	checkBins();
    }

    LOSintegrator ( const double xmin, const double xmax, const double ymin, const double ymax, const double zmin, const double zmax, const std::vector<double> &Rbins, const std::vector<double> &camera, const double stepSize, const double relTolerance, const double minStepSize ) :
            fxmin ( xmin ),
            fxmax ( xmax ),
            fymin ( ymin ),
            fymax ( ymax ),
            fzmin ( zmin ),
            fzmax ( zmax ),
            fcamera ( camera ),
            fRbins ( Rbins ),
            fstepSize ( stepSize ),
            fminStepSize ( minStepSize ),
            frelTolerance ( relTolerance ),
            fabsTolerance ( 0 ),
            fdimensions ( 3 ) {
        if ( fzmin >= fzmax )
            throw ( LOSboundError ( "zmin must be < zmax" ) );
        if ( fxmin >= fxmax )
            throw ( LOSboundError ( "xmin must be < xmax" ) );
        if ( fymin >= fymax )
            throw ( LOSboundError ( "ymin must be < ymax" ) );
	checkBins();
    }

    /** \brief Integrate all functions funcs along the line of sight in the directon l,b
     *
     * Currently we use the trapesoid method with a fixed stepsize.
     * 
     * absorption should provide the absorption coefficient ( in units cm^-1 ).
     *
     * The return double vector array has the dimension (funcs.size(), Rbins.size()-1) if Rbins encompasses the entire integration region.
     * If Rbins.size() < 2, the second dimension is 1.  If Rbins does not encompass the entire region, the size is Rbins.size().
     * The region not in one of the defined Rbins is always returned last in the array
     */
    std::vector<std::vector<T> > integrate ( const double l, const double b, 
          const std::vector< std::unique_ptr< LOSfunction<T> > > &funcs, 
          const std::unique_ptr< LOSfunction<T> > &absorption) const;
    std::vector<std::vector<T> > integrate ( const double l, const double b, 
          const std::vector< std::unique_ptr< LOSfunction<T> > > &funcs) const { return integrate(l, b, funcs, nullptr); }
};

template <typename T>
void LOSintegrator<T>::output(T t) const {
   std::cout<<t<<std::endl;
}

template <> void SM::LOSintegrator<std::valarray<double> >::output (std::valarray<double> t) const;

/*
template <typename T>
inline void LOSintegrator<T>::coordinates( const double d, const double cosl, const double sinl, const double cosb, const double sinb, double &x, double &y, double &z ) const{
   x = fcamera[0] - d*cosb*cosl;
   y = fcamera[1] - d*cosb*sinl;
   z = fcamera[2] + d*sinb;
}
*/

template <typename T>
inline size_t LOSintegrator<T>::findAnnulus(double R) const {
   size_t il(0);
   size_t iu(fRbins.size()-1);
   while (iu-il>1){
      size_t i = (il+iu)/2;
      if ( R < fRbins[i] ){
	 iu = i;
      } else {
	 il = i;
      }
   }
   return il;
}

template <typename T>
void LOSintegrator<T>::checkBins() {
   //Check for empty rbins
   if (fRbins.size() > 1) {
      //Make sure they are monotonically increasing
      if (fRbins[0] < 0) {
         throw(LOSargumentError("Rbins must be >= 0"));
      }
      for ( size_t i = 0; i < fRbins.size()-1; ++i ) {
         if ( fRbins[i] >= fRbins[i+1] )
            throw(LOSargumentError("Rbins must be monotonically increasing"));
      }
      //Find the number of Rbins we need in the integration
      double rmax;
      if (fdimensions == 2) {
         rmax = fRmax;
      } else {
         const double xmax = std::max(fabs(fxmax), fabs(fxmin));
         const double ymax = std::max(fabs(fymax), fabs(fymin));
         rmax = sqrt(xmax*xmax+ymax*ymax);
      }

      if (fRbins[fRbins.size()-1] < rmax) {
         fnRbins = fRbins.size();
         faddBoundaries = 1;
      } else if (fRbins[0] > 0) {
         fnRbins = fRbins.size();
         faddBoundaries = 0;
      } else {
         fnRbins = fRbins.size() - 1;
         faddBoundaries = 1;
      }
   } else {
      fnRbins = 1;
      fRbins.resize(1);
      fRbins[0] = 1; //Just has to be greater than 0
      faddBoundaries = 0;
   }

}

template <typename T>
bool LOSintegrator<T>::lessLimit(const std::pair< std::pair<double,double>, int> &a, const std::pair< std::pair<double,double>, int> &b) {
    return a.first.first < b.first.first;
}

template <typename T>
void LOSintegrator<T>::simFun::coordinates( const double s, double &x, double &y, double &z) const {
   x = fcamera[0] - s*fcosb*fcosl;
   y = fcamera[1] - s*fcosb*fsinl;
   z = fcamera[2] + s*fsinb;
}

template <typename T>
std::vector<T> LOSintegrator<T>::simFun::operator () (double s) const {

   //Calculate the correct coordinates
   double x(0), y(0), z(0);
   coordinates(s, x, y, z);

   //Keep track of absorption, does not work when simvec uses cache
   if (fabsorption.get() != nullptr) {
      if (itau == 5) {
         if (s > stau[4]) {
            tau[0] = tau[4];
            stau[0] = stau[4];
         }
         itau = 1;
      }
      tau[itau] = (*fabsorption)(x,y,z,fdir)*(s-stau[itau-1]) * kpc2cm + tau[itau-1];
      stau[itau] = s;
   }

   //Calculate the function values
   std::vector<T> output;
   output.reserve(ffuncs.size());
   for (size_t i = 0; i < ffuncs.size(); ++i) {
      output.push_back((*ffuncs[i])(x,y,z,fdir));
   }

   //Correct for absorption
   if (fabsorption.get() != nullptr) {
      for (size_t i = 0; i < ffuncs.size(); ++i) {
         output[i] *= exp(-tau[itau]);
      }
      ++itau;
   }

   return output;
}

template <typename T>
std::vector< std::vector<T> > LOSintegrator<T>::integrate ( const double l, const double b, 
      const std::vector< std::unique_ptr< LOSfunction<T> > > &funcs, 
      const std::unique_ptr< LOSfunction<T> > &absorption) const {

   // Calculate the unit vector for the camera direction
   const double dtr= M_PI/180.;                           // conversion degrees to radians
   const double sinb=sin ( b*dtr );
   const double cosb=cos ( b*dtr );
   const double sinl=sin ( l*dtr );
   const double cosl=cos ( l*dtr );

   const vec3 dir(-cosb*cosl, -cosb*sinl, sinb);

   // Initialize the output vector.  As we
   // expect valarrays we have to initialize on copy constructor rather than
   // assignment.
   std::vector< std::vector<T> > integ(funcs.size());
   for ( size_t i = 0; i < funcs.size(); ++i ) {
      T tmp((*funcs[i])(0.,0.,0.,dir));
      tmp = 0;
      integ[i].resize(fnRbins, tmp);
   }

   //Find the boundaries for this line of sight
   std::vector< std::pair<double,double> > boundaries(fnRbins+faddBoundaries);
   if (fdimensions == 2) {
      boundaries[boundaries.size()-1] = IntersectionCylinder(l,b,fRmax, fzmin, fzmax, fcamera);
   } else {
      boundaries[boundaries.size()-1] = IntersectionBox(l,b,fxmin,fxmax,fymin,fymax,fzmin,fzmax,fcamera);
   }

   const std::pair<double,double> &outerBound = boundaries[boundaries.size()-1];

   //If both are negative, return the empty array
   if (outerBound.first < 0 && outerBound.second < 0)
      return integ;

   //Find the boundaries for the rings, if any
   for (size_t i = 0; i < boundaries.size() - 1; ++i) {
      boundaries[i] = IntersectionCylinder(l,b,fRbins[i], fzmin, fzmax, fcamera);
   }

   //Now we will merge the boundaries together to add everything into integration limits for the rings
   std::vector< std::pair< std::pair<double,double>, int > > limits;

   //If the innermost ring doesn't cover 0, add that part to the last bin if needed
   if ( fRbins[0] > 0 && ! ( boundaries[0].first < 0 && boundaries[0].second < 0) ) {
      //Make sure they are within outer boundary
      if (boundaries[0].first < outerBound.second && boundaries[0].second > outerBound.first) {
         std::pair<double,double> limit;
         limit.first = std::max(boundaries[0].first,outerBound.first);
         limit.second = std::min(boundaries[0].second,outerBound.second);
         limits.push_back(std::pair< std::pair<double,double>, int> (limit, fnRbins-1));
      }
   }

   //Loop over the rings and add them to the limits as needed
   //Note that this also takes care of outermost boundary if Rbins does not encompass the integration region
   for (size_t i = 0; i < boundaries.size()-1; ++i) {
      //If los does not pass outer boundary of ring, nothing should be done
      if (boundaries[i+1].first < 0 && boundaries[i+1].second < 0)
         continue;

      // Check that at least part of the ring is within the boundaries of the box
      if (boundaries[i+1].second < outerBound.first)
         continue;
      if (boundaries[i+1].first > outerBound.second)
         continue;

      //If it does not pass inner boundary, just use outer boundary, making sure to be within outer boundaries of the integration zone
      if (boundaries[i].first < 0 && boundaries[i].second < 0) {

         std::pair<double,double> limit;
         limit.first = std::max(boundaries[i+1].first, outerBound.first);
         limit.second = std::min(boundaries[i+1].second, outerBound.second);
         limits.push_back(std::pair< std::pair<double,double>, int> (limit, i));
      
      } else {

         //Need to add first leg if first intersection of inner boundary is larger than first intersecion of outer boundary (i.e. not both 0)
         if (boundaries[i+1].first < boundaries[i].first) {
            //Make sure we are within bounds of outer box
            if ( boundaries[i].first > outerBound.first || boundaries[i+1].first < outerBound.second ) {
               std::pair<double,double> limit;
               limit.first = std::max(outerBound.first, boundaries[i+1].first);
               limit.second = std::min(outerBound.second, boundaries[i].first);
               limits.push_back(std::pair< std::pair<double,double>, int> (limit, i));
            }
         }

         //Similarily, if the second intersection of inner boundary is smaller than second interseciont of outer boundary, add it
         if (boundaries[i].second < boundaries[i+1].second) {
            //Make sure we are within bounds of outer box
            if ( boundaries[i+1].second > outerBound.first && boundaries[i].second < outerBound.second  ) {
               std::pair<double,double> limit;
               limit.first = std::max(outerBound.first, boundaries[i].second);
               limit.second = std::min(outerBound.second, boundaries[i+1].second);
               limits.push_back(std::pair< std::pair<double,double>, int> (limit, i));
            }
         }
      }
   }

   //Sort the limits to get extinction calculation correct
   std::sort(limits.begin(), limits.end(), LOSintegrator<T>::lessLimit);

/*
   for (size_t i = 0; i < boundaries.size(); ++i)
      std::cout<<i<<", "<<boundaries[i].first<<", "<<boundaries[i].second<<std::endl;

   for (size_t j = 0; j < limits.size(); ++j) 
      std::cout<<limits[j].second<<", "<<j<<", "<<limits[j].first.first<<", "<<limits[j].first.second<<std::endl;
*/

   //Use the Simpsons variable step integrator to integrate.

   simFun fun(cosl, cosb, sinl, sinb, fcamera, funcs, absorption);
   for (size_t j = 0; j < limits.size(); ++j) {
      //Fix the step size so we won't have many fine steps on the end of the LOS.  Make it such that the step size is always reduced for better accuracy
      const double step = (limits[j].first.second - limits[j].first.first)/(int((limits[j].first.second - limits[j].first.first)/fstepSize) + 1);
      //Only use cache when not using absorption
      std::vector<T> out = simvec( limits[j].first.first, limits[j].first.second, step, fminStepSize, frelTolerance, fabsTolerance, fun, absorption.get() == nullptr);
      for (size_t i = 0; i < funcs.size(); ++i) {
         integ[i][limits[j].second] += out[i];
      }
   }

   //Scale it with kpc2cm to get the units correct
   for ( size_t i = 0; i < funcs.size(); ++i )
      for ( size_t j = 0; j < integ[i].size(); ++j )
         integ[i][j] *= kpc2cm;
   
   return integ;
}

/*
 * Store the boundary method for future reference.  It might be needed in the future
 *
       // Perform check for a ring boundary if dl > 0.
       // Note that this splits the integration value between the annuli as if the pixels
       // where square shaped
       if ( dl > 0 && ( l < 90 || l > 270 ) ) {
	  double lmin = l-dl*0.49;
	  double lmax = l+dl*0.49;
	  double xmin, ymin, zmin, xmax, ymax, zmax;
	  coordinates(d, cos(lmin*dtr), sin(lmin*dtr), cosb, sinb, xmin, ymin, zmin);
	  coordinates(d, cos(lmax*dtr), sin(lmax*dtr), cosb, sinb, xmax, ymax, zmax);
	  double Rmin = sqrt(xmin*xmin+ymin*ymin);
	  double Rmax = sqrt(xmax*xmax+ymax*ymax);

      	  size_t minAnn, maxAnn;
	  minAnn = findAnnulus(Rmin);
	  maxAnn = findAnnulus(Rmax);

	  // Find the edge of the ring, to make sure we only apply this
	  // correction for tangential split of pixels
	  double ledge(0);
	  if ( l > 270. )
	     ledge=360 - asin ( fRbins[maxAnn+1]/abs(fxorig) ) / dtr;
	  else
    	     ledge=asin ( fRbins[minAnn+1]/abs(fxorig) ) / dtr;

	  //std::cout<<"Edge test: "<<ledge<<", "<<lmin<<", "<<lmax<<std::endl;
	  //std::cout<<fRbins[maxAnn+1]<<", "<<fxorig<<", "<<dtr<<std::endl;

      	  if ( lmin < ledge && ledge < lmax && minAnn != maxAnn ) {
	     if ( currAnn != minAnn ) {
		//We apply the weight so that when all annuli are summed
		//up, it would give very similar results as if no splitting
		//between annuli was performed.  This is of course not
		//entirely true as we change the integration path in the
		//split pixel.
		double weight = (ledge-lmin)/(lmax-lmin);
		if ( noSplit ) { // Last step we had no split
		   for ( size_t i = 0; i < funcs.size(); ++i ) {
		      T tmp((*funcs[i])(xx,yy,zz));
		      integ[i][minAnn] += (*funcs[i])(xmin,ymin,zmin)*weight;
		      integ[i][currAnn] += tmp*(1-weight);
		      integ[i][prevAnn] += tmp;
		   }
		} else {
		   for ( size_t i = 0; i < funcs.size(); ++i ) {
		      T tmp((*funcs[i])(xx,yy,zz));
		      integ[i][minAnn] += (*funcs[i])(xmin,ymin,zmin)*weight*2;
		      integ[i][currAnn] += tmp*(1-weight);
		      integ[i][prevAnn] += tmp*(1-weight);
		   }
		}
	     } else {
		double weight = (lmax-ledge)/(lmax-lmin);
		if ( noSplit ) { // Last step we had no split
		   for ( size_t i = 0; i < funcs.size(); ++i ) {
		      T tmp((*funcs[i])(xx,yy,zz));
		      integ[i][maxAnn] += (*funcs[i])(xmax,ymax,zmax)*weight;
		      integ[i][currAnn] += tmp*(1-weight);
		      integ[i][prevAnn] += tmp;
		   }
		} else {
		   for ( size_t i = 0; i < funcs.size(); ++i ) {
		      T tmp((*funcs[i])(xx,yy,zz));
		      integ[i][maxAnn] += (*funcs[i])(xmax,ymax,zmax)*weight*2;
		      integ[i][currAnn] += tmp*(1-weight);
		      integ[i][prevAnn] += tmp*(1-weight);
		   }
		}
	     }
	     noSplit = false;
	  } else {
	     //There should be a final addition to the split ring, but it
	     //is only a minor error.
	     noSplit = true;
	  }

       } //End split pixel test

       if ( noSplit ) {
	  // Perform the unsplit integration.  Since we are using the trapezoid
	  // rule, we add the value to both previous and current annuli.
	  // If the are the same, we get twice the value as we are supposed
	  // to.  Otherwise, we both end and start the integral.
	  for ( size_t i = 0; i < funcs.size(); ++i ) {
	     T tmp((*funcs[i])(xx,yy,zz));
	     integ[i][currAnn] += tmp;
	     integ[i][prevAnn] += tmp;
	  }
       }

       prevAnn = currAnn;
    }
    for ( size_t i = 0; i < funcs.size(); ++i )
       for ( size_t j = 0; j < integ[i].size(); ++j )
   	  integ[i][j] *= fstepSize;

    //Now we must add the last step, assuming it is not split, which is
    //safe in most cases
    coordinates(dmax,cosl,sinl,cosb,sinb,xx,yy,zz);
    double RR = sqrt(xx*xx+yy*yy);
    //std::cout<<dmax<<", "<<xx<<", "<<yy<<", "<<zz<<", "<<RR<<std::endl;
    currAnn = findAnnulus(RR);
    for ( size_t i = 0; i < funcs.size(); ++i ) 
       integ[i][currAnn] +=(*funcs[i])(xx,yy,zz)*(dmax-d);

    //We must also multiply with kpc2cm and divide by 2 to get the correct
    //answer
    for ( size_t i = 0; i < funcs.size(); ++i )
       for ( size_t j = 0; j < integ[i].size(); ++j )
   	  integ[i][j] *= kpc2cm/2;

    return integ;
    */

} //End namespace
#endif
