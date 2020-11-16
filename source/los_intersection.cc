#include "los_intersection.h"
#include "constants.h"
#include "Units.h"
#include "PhysicalConstants.h"

#include <valarray>
#include <vector>
#include <utility>

using namespace std;

pair<double,double> SM::IntersectionCylinder(const double l, const double b, 
                                                    const double rMax, 
                                                    const double zMin, const double zMax,
                                                    const vector<double>& cameraLocation) 
{ 

   pair<double,double> s(-1.,-1.);

   const double cameraX = cameraLocation[0];
   const double cameraY = cameraLocation[1];
   const double cameraZ = cameraLocation[2];

   const double lRad = l*utl::kConvertDegreesToRadians, bRad = b*utl::kConvertDegreesToRadians;

   const double sinb = sin(bRad), cosb = cos(bRad), sinl = sin(lRad), cosl = cos(lRad);  

   // Intersection with cylindrical volume 
   // Method: solve quadratic equation for cylinder with ray intersection
   // Let point X = (xo + tpx, y0 + tpy, z0 + tpz) intersect with 
   // cylinder with equation x^2 + y^2 = rMax^2, z = +/- zMax.

   //It is already normalized so no need to create a pointing or a vector
   //const pointing co(vec3(-cosb*cosl, -cosb*sinl, sinb));

   //vec3 dir(-cosb*cosl, -cosb*sinl, sinb);

   const double pX = -cosb*cosl, pY = -cosb*sinl, pZ = sinb;

   const double cameraR2 = cameraX*cameraX + cameraY*cameraY;

   const double C = (cameraR2 - rMax*rMax);

   const double A = pX*pX + pY*pY;

   const double B = 2.*(pX*cameraX + pY*cameraY);

   const double discriminant = B*B - 4.*A*C;

   // No real roots, no intersection with cylinder

   if (discriminant < 0.) 
      return s;

   // If A == 0 we are looking straight up and have no intersection with the cylinder

   if ( A > 0. ) {  // TODO We should have a range to prevent buffer overflow but I am lazy

      // t2 <= t1 because discriminant >= 0 and A > 0
      const double t1 = 0.5/A*(-B + sqrt(discriminant));

      const double t2 = 0.5/A*(-B - sqrt(discriminant));

      // Both negative means there is no intersection in the right direction
      if (t1 < 0 && t2 < 0)
         return s;

      const double z1 = cameraZ + t1*pZ, z2 = cameraZ + t2*pZ;

      // If both points above/below the cylinder there is no intersection

      if (z1 > zMax && z2 > zMax)
         return s;

      if (z1 < zMin && z2 < zMin)
         return s;

      // If both points are within the z boundaries, t1 and t2 are the answer

      s.first = t2;
      s.second = t1;


      // Find the intersection with the z boundaries that should be used

      if (z1 > zMax) {

         s.second = (zMax - cameraZ)/pZ;

      } else if (z1 < zMin) {

         s.second = (zMin - cameraZ)/pZ;

      }

      if (z2 > zMax) {

         s.first = (zMax - cameraZ)/pZ;

      } else if (z2 < zMin) {

         s.first = (zMin - cameraZ)/pZ;

      }

   } else {// Looking up/down

      // If camera is not within rmax, no intersection
      if (C > 0) 
         return s;

      const double t1 = zMin - cameraZ;
      const double t2 = zMax - cameraZ;

      // Both are negative so we have no intersection
      if (t1 < 0. && t2 < 0.) 
         return s;

      if ( sinb > 0 ) {//Up, so zmin should give smaller distance

         s.first = t1;
         s.second = t2;

      } else {

         s.first = t2;
         s.second = t1;

      }
   }

   // If first coordinate less than 0 then we put it to 0 like it should be

   if (s.first < 0)
      s.first = 0;

   return s;
}


pair<double,double> SM::IntersectionBox(const double l, const double b, 
                                               const double xMin, const double xMax,
                                               const double yMin, const double yMax,
                                               const double zMin, const double zMax,
                                               const vector<double>& cameraLocation) 
{ 

   pair<double,double> s(-1.,-1.);

   const double cameraX = cameraLocation[0];
   const double cameraY = cameraLocation[1];
   const double cameraZ = cameraLocation[2];

   const double lRad = l*utl::kConvertDegreesToRadians, bRad = b*utl::kConvertDegreesToRadians;

   const double sinb = sin(bRad), cosb = cos(bRad), sinl = sin(lRad), cosl = cos(lRad);  

   // Intersection with box volume
   // Method: equation for plane with ray intersection for all faces of box
   // Let point X = (xo + tpx, y0 + tpy, z0 + tpz) intersect with 
   // plane with equation Ax + By + Cz + D for all faces of box.
   // Use algorithm from `An efficient and robust ray-box intersection
   // algorithm' -- Williams et al. http://www.cs.utah.edu/~awilliam/box/

   //Since it is already normalized, just use constants pX, pY, and pZ
   //const pointing co(vec3(-cosb*cosl, -cosb*sinl, sinb));

   //vec3 dir = co.to_vec3();

   const double pX = -cosb*cosl, pY = -cosb*sinl, pZ = sinb;

   vector<double> invDir;
   invDir.reserve(3);

   vector<size_t> index;
   index.reserve(3);

   // Check for directions that are parallel to the planes, possibly causing error
   // TODO Make this test do it within numerical accuracy to avoid overflow
   if (pX != 0) {

      invDir.push_back(1./pX);
      index.push_back(0);

   } else if (cameraX > xMax || cameraX < xMin)
      return s;

   if (pY != 0){

      invDir.push_back(1./pY);
      index.push_back(1);

   } else if (cameraY > yMax || cameraY < yMin)
      return s;

   if (pZ != 0) {

      invDir.push_back(1./pZ);
      index.push_back(2);

   } else if (cameraZ > zMax || cameraZ < zMin)
      return s;

   valarray<int> sign(0, index.size());

   for (size_t i = 0; i < index.size(); ++i) 
      sign[i] = invDir[i] < 0;

   valarray< valarray<double> > boxSize;

   boxSize.resize(2);
   boxSize[0].resize(3);
   boxSize[1].resize(3);

   boxSize[0][0] = xMin, boxSize[0][1] = yMin, boxSize[0][2] = zMin;

   boxSize[1][0] = xMax, boxSize[1][1] = yMax, boxSize[1][2] = zMax;

   valarray<double> tMin(0., 2), tMax(0., 2);

   tMin[0] = (boxSize[sign[0]][index[0]] - cameraLocation[index[0]])*invDir[0];
   tMax[0] = (boxSize[1 - sign[0]][index[0]] - cameraLocation[index[0]])*invDir[0];

   for (size_t i = 1; i < index.size(); ++i) {

      tMin[1] = (boxSize[sign[i]][index[i]] - cameraLocation[index[i]])*invDir[i];
      tMax[1] = (boxSize[1 - sign[i]][index[i]] - cameraLocation[index[i]])*invDir[i];

      if ((tMin[0] > tMax[1]) || (tMin[1] > tMax[0]))
         return s; // No intersection

      if (tMin[1] > tMin[0])
         tMin[0] = tMin[1];

      if (tMax[1] < tMax[0])
         tMax[0] = tMax[1];

   }

   //Since tMin always <= tMax, only need to test tMax
   if (tMax[0] < 0.)
      return s;

   s.first = tMin[0] < 0 ? 0 : tMin[0];
   s.second = tMax[0];

   return s;

}


