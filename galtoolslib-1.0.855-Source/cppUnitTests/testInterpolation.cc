#include "testInterpolation.h"
#include <Interpolation.h>
#include <vector>
#include <valarray>
#include <stdexcept>

CPPUNIT_TEST_SUITE_REGISTRATION( testInterpolation );

using namespace utl;

void testInterpolation::setUp() {
   
   //Create the xValues for the interpolation
   xValues.resize(100);
   for (size_t i(0); i < xValues.size(); ++i)
      xValues[i] = 0.1 + 5.*sqrt(double(i));

}

void testInterpolation::tearDown() {
}

void testInterpolation::index() {

   //Create a vector to search
   std::vector<double> values(10);
   values[0] = -3.2;
   values[1] = -1.2;
   values[2] = -0.2;
   values[3] = 3.2;
   values[4] = 4;
   values[5] = 6.2;
   values[6] = 7.2;
   values[7] = 8.2;
   values[8] = 9.2;
   values[9] = 13.2;

   //Check couple of values of index below
   size_t expected, value;
   expected = 1;
   value = utl::indexBelow(values, -0.5);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   expected = 7;
   value = utl::indexBelow(values, 8.5);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //Similar for index above
   expected = 3;
   value = utl::indexAbove(values, 3.0);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   expected = 5;
   value = utl::indexAbove(values, 5.0);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //Now check boundaries
   expected = 0;
   value = utl::indexBelow(values, -5.0);
   CPPUNIT_ASSERT_EQUAL(expected, value);
   value = utl::indexAbove(values, -5.0);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   expected = 9;
   value = utl::indexBelow(values, 50.);
   CPPUNIT_ASSERT_EQUAL(expected, value);
   value = utl::indexAbove(values, 50.);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //Now make sure that indexBelow returns values <= and index above >=
   expected = 4;
   value = utl::indexBelow(values, 4);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   expected = 4;
   value = utl::indexAbove(values, 4);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //Test that the range is correct for index above and below
   size_t lowExp, highExp, low, high;
   lowExp = 2; highExp = 5;
   utl::findIndexRange(values, 2., 5., low, high);
   CPPUNIT_ASSERT_EQUAL(lowExp, low);
   CPPUNIT_ASSERT_EQUAL(highExp, high);

   //Make sure the range is correct at the boundaries
   lowExp = 8; highExp = 9;
   utl::findIndexRange(values, 20., 50., low, high);
   CPPUNIT_ASSERT_EQUAL(lowExp, low);
   CPPUNIT_ASSERT_EQUAL(highExp, high);

   lowExp = 0; highExp = 1;
   utl::findIndexRange(values, -50., -20., low, high);
   CPPUNIT_ASSERT_EQUAL(lowExp, low);
   CPPUNIT_ASSERT_EQUAL(highExp, high);

   //Test for run_time errors
   CPPUNIT_ASSERT_THROW(utl::findIndexRange(values, 20., 10., low, high), std::runtime_error);

   //Do all the same tests for valarrays
   std::valarray<double> vvalues(&values[0], values.size());

   //Check couple of vvalues of index below
   expected = 1;
   value = utl::indexBelow(vvalues, -0.5);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   expected = 7;
   value = utl::indexBelow(vvalues, 8.5);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //Similar for index above
   expected = 3;
   value = utl::indexAbove(vvalues, 3.0);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   expected = 5;
   value = utl::indexAbove(vvalues, 5.0);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //Now check boundaries
   expected = 0;
   value = utl::indexBelow(vvalues, -5.0);
   CPPUNIT_ASSERT_EQUAL(expected, value);
   value = utl::indexAbove(vvalues, -5.0);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   expected = 9;
   value = utl::indexBelow(vvalues, 50.);
   CPPUNIT_ASSERT_EQUAL(expected, value);
   value = utl::indexAbove(vvalues, 50.);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //Now make sure that indexBelow returns vvalues <= and index above >
   expected = 4;
   value = utl::indexBelow(vvalues, 4);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   expected = 4;
   value = utl::indexAbove(vvalues, 4);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //Test that the range is correct for index above and below
   lowExp = 2; highExp = 5;
   utl::findIndexRange(vvalues, 2., 5., low, high);
   CPPUNIT_ASSERT_EQUAL(lowExp, low);
   CPPUNIT_ASSERT_EQUAL(highExp, high);

   //Make sure the range is correct at the boundaries
   lowExp = 8; highExp = 9;
   utl::findIndexRange(vvalues, 20., 50., low, high);
   CPPUNIT_ASSERT_EQUAL(lowExp, low);
   CPPUNIT_ASSERT_EQUAL(highExp, high);

   lowExp = 0; highExp = 1;
   utl::findIndexRange(vvalues, -50., -20., low, high);
   CPPUNIT_ASSERT_EQUAL(lowExp, low);
   CPPUNIT_ASSERT_EQUAL(highExp, high);

   //Test for run_time errors
   CPPUNIT_ASSERT_THROW(utl::findIndexRange(vvalues, 20., 10., low, high), std::runtime_error);

}

double testInterpolation::linearFunction(double x) {
   return 0.4*x - 3.2;
}

double testInterpolation::powerLawFunction(double x) {
   return 1.43e-3 * pow(x/3.4, -0.432);
}

void testInterpolation::linear() {

   //Create the yValues
   std::vector<double> yValues(xValues.size());

   for (size_t i(0); i < xValues.size(); ++i)
      yValues[i] = linearFunction(xValues[i]);
   
   //Calculate the range of the x values so we can test within range and outside of range
   const double r = xValues[xValues.size()-1] - xValues[0];

   utl::LinearInterpolation lint(xValues, yValues);

   //Do a test for arbitrary values in the range
   const size_t N(10);
   for (size_t i(0); i < N; ++i) {
      const double x = xValues[0] + 0.1*r + 0.8*r/(N-1)*i;
      const double expected = linearFunction(x);
      const double value = lint(x);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, value, 1e-8);
   }

   //Extrapolation to both high and low values
   double x = xValues[0] - 0.5*r;
   double expected = linearFunction(x);
   double value = lint(x);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, value, 1e-8);

   x = xValues[0] + 1.5*r;
   expected = linearFunction(x);
   value = lint(x);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, value, 1e-8);

}

void testInterpolation::powerlaw() {

   //Create the yValues
   std::vector<double> yValues(xValues.size());

   for (size_t i(0); i < xValues.size(); ++i)
      yValues[i] = powerLawFunction(xValues[i]);
   
   //Calculate the range of the x values so we can test within range and outside of range
   const double r = xValues[xValues.size()-1] - xValues[0];

   utl::PowerLawInterpolation pint(xValues,yValues);

   //Do a test for arbitrary values in the range
   const size_t N(10);
   for (size_t i(0); i < N; ++i) {
      const double x = xValues[0] + 0.1*r + 0.8*r/(N-1)*i;
      const double expected = powerLawFunction(x);
      const double value = pint(x);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, value, 1e-8);
   }

   //Extrapolation to both high and low values
   double x = xValues[0] * 0.5;
   double expected = powerLawFunction(x);
   double value = pint(x);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, value, 1e-8);

   x = xValues[xValues.size()-1] * 2.;
   expected = powerLawFunction(x);
   value = pint(x);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, value, 1e-8);

   //Make sure the exceptions are thrown
   CPPUNIT_ASSERT_THROW(pint(-1), std::domain_error);

   yValues[0] = -3.;
   utl::PowerLawInterpolation pint1(xValues,yValues);
   CPPUNIT_ASSERT_THROW(pint1(xValues[0]), std::domain_error);

   yValues[0] = powerLawFunction(xValues[0]);
   xValues[0] = -1.2;
   utl::PowerLawInterpolation pint2(xValues,yValues);
   CPPUNIT_ASSERT_THROW(pint2(xValues[1]*0.5), std::domain_error);

}

void testInterpolation::lagrange() {

}
