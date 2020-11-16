#include "Coordinate.h"
#include "pointing.h"
#include "PhysicalConstants.h"
#include <iostream>

#include <ErrorLogger.h>

#include <trafos.h>

//Create static instances to do the conversion.  Currently assume J2000 always.
static Trafo tcg(J2000, J2000, Equatorial, Galactic);
static Trafo tgc(J2000, J2000, Galactic, Equatorial);
static Trafo tce(J2000, J2000, Equatorial, Ecliptic);
static Trafo tec(J2000, J2000, Ecliptic, Equatorial);
static Trafo teg(J2000, J2000, Ecliptic, Galactic);
static Trafo tge(J2000, J2000, Galactic, Ecliptic);

void SM::eq_gal2(double mj, double DEC, double RA, double *b, double *l){
   pointing in(utl::kPiOnTwo-DEC, RA);
   auto out = tcg(in);
   *b = utl::kPiOnTwo-out.theta;
   *l = out.phi;
}
void SM::eq_ecl2(double mj, double DEC, double RA, double *lt, double *lg){
   pointing in(utl::kPiOnTwo-DEC, RA);
   auto out = tce(in);
   *lt = utl::kPiOnTwo-out.theta;
   *lg = out.phi;
}
void SM::gal_eq2(double mj, double b, double l, double *DEC, double *RA){
   pointing in(utl::kPiOnTwo-b, l);
   auto out = tgc(in);
   *DEC = utl::kPiOnTwo-out.theta;
   *RA = out.phi;
}
void SM::ecl_eq2(double mj, double lt, double lg, double *DEC, double *RA){
   pointing in(utl::kPiOnTwo-lt, lg);
   auto out = tec(in);
   *DEC = utl::kPiOnTwo-out.theta;
   *RA = out.phi;
}
//Functions to convert from gal to ecl and wise versa
void SM::gal_ecl2(double mj, double b, double l, double *lt, double *lg){
   pointing in(utl::kPiOnTwo-b, l);
   auto out = tge(in);
   *lt = utl::kPiOnTwo-out.theta;
   *lg = out.phi;
}
void SM::ecl_gal2(double mj, double lt, double lg, double *b, double *l){
   pointing in(utl::kPiOnTwo-lt, lg);
   auto out = teg(in);
   *b = utl::kPiOnTwo-out.theta;
   *l = out.phi;
}

double SM::chacos(double x) {
   if (x <= -1)
      return utl::kPi;
   if (x >= 1)
      return 0;
   return std::acos(x);
}

SM::Coordinate::Coordinate(const double l, const double b) {
   m_lon = l*utl::kConvertDegreesToRadians;
   m_lat = b*utl::kConvertDegreesToRadians;
   m_cSys = CoordSys::GAL;
   normalize(m_lon, m_lat);
}

void SM::Coordinate::normalize(double &lon, double &lat) {
   //Normalize the latitude range to -pi,pi
   while (lat < -utl::kPi)
      lat += utl::kTwoPi;
   while (lat > utl::kPi)
      lat -= utl::kTwoPi;

   //Now normalize it to -pi/2, pi/2
   if (lat > utl::kPiOnTwo) {
      lat = utl::kPi-lat;
      lon += utl::kPi;
   }
   if (lat < -utl::kPiOnTwo) {
      lat = -utl::kPi-lat;
      lon += utl::kPi;
   }

   //Normalize the longitude range to 0,2pi
   while (lon < 0)
      lon += utl::kTwoPi;
   while (lon >= utl::kTwoPi)
      lon -= utl::kTwoPi;

}


SM::Coordinate::Coordinate(const pointing &point) {
   m_lat = utl::kPiOnTwo - point.theta;
   m_lon = point.phi;
   m_cSys = CoordSys::GAL;
}

SM::Coordinate::Coordinate(double lon, double lat, CoordSys sys) : m_lon(lon), m_lat(lat), m_cSys(sys)
{
   normalize(m_lon, m_lat);
}

SM::Coordinate::Coordinate(const pointing &point, CoordSys sys) {
   m_lat = utl::kPiOnTwo - point.theta;
   m_lon = point.phi;
   m_cSys = sys;
}

void SM::Coordinate::getCoordinates(double &lon, double &lat, CoordSys sys) const
{
   switch (m_cSys) {
      case CoordSys::GAL:
         switch (sys) {
            case CoordSys::GAL:
               lon = m_lon;
               lat = m_lat;
               break;
            case CoordSys::EQ:
               gal_eq2(J2000, m_lat, m_lon, &lat, &lon);
               break;
            case CoordSys::ECL:
               gal_ecl2(J2000, m_lat, m_lon, &lat, &lon);
               break;
         }
         break;
      case CoordSys::EQ:
         switch (sys) {
            case CoordSys::GAL:
               eq_gal2(J2000, m_lat, m_lon, &lat, &lon);
               break;
            case CoordSys::EQ:
               lon = m_lon;
               lat = m_lat;
               break;
            case CoordSys::ECL:
               eq_ecl2(J2000, m_lat, m_lon, &lat, &lon);
               break;
         }
         break;
      case CoordSys::ECL:
         switch (sys) {
            case CoordSys::GAL:
               ecl_gal2(J2000, m_lat, m_lon, &lat, &lon);
               break;
            case CoordSys::EQ:
               ecl_eq2(J2000, m_lat, m_lon, &lat, &lon);
               break;
            case CoordSys::ECL:
               lon = m_lon;
               lat = m_lat;
               break;
         }
         break;
   }
   normalize(lon, lat);
}

pointing SM::Coordinate::healpixAng(CoordSys sys) const
{
   double lon(0), lat(0);
   getCoordinates(lon, lat, sys);

   auto ang = pointing(utl::kPiOnTwo-lat, lon);
   ang.normalize();

   return ang;
}

pointing SM::Coordinate::healpixAng() const{
   auto ang = pointing(utl::kPiOnTwo-m_lat, m_lon);
   ang.normalize();
   return ang;
}

double SM::Coordinate::l() const {
   double lon(0), lat(0);
   getCoordinates(lon, lat, CoordSys::GAL);

   return lon*utl::kConvertRadiansToDegrees;
}

double SM::Coordinate::b() const {
   double lon(0), lat(0);
   getCoordinates(lon, lat, CoordSys::GAL);

   return lat*utl::kConvertRadiansToDegrees;
}


namespace SM {
std::ostream & operator << (std::ostream &os, const SM::Coordinate &coord) {
   os << "(" << coord.l() <<","<<coord.b()<<")";
   return os;
}
}
