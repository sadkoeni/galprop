#ifndef SPATIAL_DISTRIBUTION_H
#define SPATIAL_DISTRIBUTION_H

#include <Parameters.h>

class SpatialDistribution
{
   public:
      /** The only constructor
       *
       * All tunable parameters should be read from the pars object.
       * See the derived class documentation for details.
       */
      SpatialDistribution(const utl::Parameters &pars) {}
      virtual ~SpatialDistribution() {}

      /** Access the distribution in an azimuthally averaged fashion.
       *
       * Default operator takes the azimuthal average of the 3D distribution.
       * Override for 2D distributions and distributions with analytic form
       * for the azimuthal average.
       *
       * \param r is the radial part of cylindrical coordinates centered on the GC
       * \param z is the distance above the galactic plane, positive z in the northern hemisphere
       *
       * The units for the distances are in kpc.
       */
      virtual double operator() (double r, double z) const;

      /** Access the distribution in 3D.
       *
       * Derived classes must define this.
       *
       * \param x, \param y, and \param z are part of a right handed Cartesian coordinate system
       * centered on the GC.  The sun is at positive x and postive z are in the northern hemisphere.
       *
       * Note that in this system, the x-axis points towards Galactic longitudes of 180 degrees and the
       * y-axis to longitudes of 270 degrees.
       *
       * The units for the distances are in kpc.
       */
      virtual double operator() (double x, double y, double z) const = 0;

};

#endif
