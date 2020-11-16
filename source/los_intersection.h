#ifndef los_intersection_h
#define los_intersection_h

#include <vector>
#include <utility>

using namespace std;

namespace SM {
/** \bried Helper functions to calculate intersections with a cylinder and a box
 * 
 * Returns the distance along the line of sight for the pair of points intersecting the shape
 * l=0 points in the negative x direction, b=pi/2 points in the z direction 
 * So shifted Galactic coordinates assuming the sun is at positive x
 * First part of pair is always smaller of the intersection points.
 * If camera is within shape, return 0 for the first coordinate
 * Returns negative values in the pair if no intersection can be found (-1,-1)
 */
pair<double,double> IntersectionCylinder(const double l, const double b,
                                                const double rMax,
                                                const double zMin, const double zMax,
                                                const vector<double>& cameraLocation);

pair<double,double> IntersectionBox(const double l, const double b,
                                           const double xMin, const double xMax,
                                           const double yMin, const double yMax,
                                           const double zMin, const double zMax,
                                           const vector<double>& cameraLocation);
}
#endif
