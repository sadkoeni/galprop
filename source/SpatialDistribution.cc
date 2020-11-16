#include "SpatialDistribution.h"

#include <Registry.h>

#include <cmath>

double SpatialDistribution::operator() (double r, double z) const
{
   //Use as a starting point at least 4 points or a maximum distance of 500 pc
   size_t nPoints = size_t(ceil(M_PI*2*r/500.));
   if (nPoints < 4)
      nPoints = 4;

   //Calculate the average over the points
   double newAverage(0);
   for (size_t i(0); i < nPoints; ++i) {
      const double phi = 2*M_PI*double(i)/double(nPoints);
      const double x = r*cos(phi);
      const double y = r*sin(phi);
      newAverage += (*this)(x,y,z);
   }
   newAverage /= nPoints;

   //Double the points and stop once the relative difference is less than 10^-4
   double oldAverage;
   do {
      oldAverage = newAverage;
      newAverage=0;
      for (size_t i(0); i < nPoints; ++i) {
         const double phi = 2*M_PI*(double(i)+0.5)/double(nPoints);
         const double x = r*cos(phi);
         const double y = r*sin(phi);
         newAverage += (*this)(x,y,z);
      }
      newAverage /= nPoints;
      newAverage = (newAverage + oldAverage)*0.5;
      nPoints *= 2;
   } while ( fabs((oldAverage - newAverage)/newAverage) > 1e-4 );

   return newAverage;
}
