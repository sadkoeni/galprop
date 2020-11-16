#ifndef Coordinate_h
#define Coordinate_h
#include "PhysicalConstants.h"
#include "pointing.h"

//Use Trafo from healpix to do coordinate conversion.
#define J2000 2000.

namespace SM {
   //Keep these functions for backwards compatibility
   void eq_gal2(double mj, double DEC, double RA, double *b, double *l);
   void eq_ecl2(double mj, double DEC, double RA, double *lt, double *lg);
   void gal_eq2(double mj, double b, double l, double *DEC, double *RA);
   void ecl_eq2(double mj, double lt, double lg, double *DEC, double *RA);
   //Functions to convert from gal to ecl and wise versa
   void gal_ecl2(double mj, double b, double l, double *lt, double *lg);
   void ecl_gal2(double mj, double lt, double lg, double *b, double *l);

   //! acos function that checks input boundaries
   double chacos(double x);

   /** \brief Specifies the coordinate system */
   enum class CoordSys {EQ, GAL, ECL};

   /** \brief Coordinates for easy access of skymaps
    *
    * Handles conversion between conventional l and b with the equator at b = 0 to healpix coordinates
    */
   class Coordinate {
      private:
         double m_lon, m_lat; //!< The coordinates, stored in radians and the current system
         //! Normalize the angles to 0,2pi and -pi/2, pi/2
         static void normalize(double &lon, double &lat);
         CoordSys m_cSys;
      public:
         /**\brief Constructor that takes in galactic coordinates (l,b) in
          * degrees.
          *
          * \param l the longitude in degrees, GC at 0
          * \param b the latitude in degrees, GC at 0
          *
          * This is for backwards compatibility.  Use of coordinator
          */
         Coordinate(const double l, const double b);
         /**\brief Default constructor initializes to 0 in Galactic coordinates*/
         Coordinate() : m_lon(0), m_lat(0), m_cSys(CoordSys::GAL) {} 
         /**\brief Constructor that takes in healpix coordinate pointing in Galactic coordinates.
          *
          * For backwards compatibility, use explicit coordinate system in future.
          */
         Coordinate(const pointing & point);

         /** \brief Create coordinate from other coordinate systems 
          *
          * lon and lat should be in radians
          */
         Coordinate(double lon, double lat, CoordSys sys);

         /** \brief Use a pointing and coordinate system */
         Coordinate(const pointing & point, CoordSys sys);

         /** \brief Return angles for another coordinate system.  
          *
          * Returns radians. 
          *
          * Throws an error if not compiled with astro and coordinate transform is requested
          */
         void getCoordinates(double &lon, double &lat, CoordSys sys) const;

         /** \brief Return pointing for another coordinate system 
          *
          * Throws an error if not compiled with astro and coordinate transform is requested
          */
         pointing healpixAng(CoordSys sys) const;


         /** \brief Return the value of l in degrees 
          *
          * This is for backwards compatibility, use getCoordinates instead.
          */
         double l() const;
         /** \brief Return the value of b in degrees
          *
          * This is for backwards compatibility, use getCoordinates instead
          */
         double b() const;

         /**\brief Get the angular coordinate for healpix
          *
          * \return a pointing object to use with healpix.  No coordinate conversion performed
          */
         pointing healpixAng(void) const;

         /**\brief Output operator
          *
          * The format is (l,b) in Galactic coordinates in units of degrees
          */
         friend std::ostream & operator << (std::ostream &os, const Coordinate & coord);

   };
}

#endif
