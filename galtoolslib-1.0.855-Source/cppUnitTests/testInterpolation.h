#ifndef testInterpolation_h
#define testInterpolation_h

#include <cppunit/extensions/HelperMacros.h>

class testInterpolation : public CppUnit::TestFixture
{
   CPPUNIT_TEST_SUITE( testInterpolation );
   CPPUNIT_TEST( index );
   CPPUNIT_TEST( linear );
   CPPUNIT_TEST( powerlaw );
   CPPUNIT_TEST( lagrange );
   CPPUNIT_TEST_SUITE_END();

   private:
      double linearFunction(double x);
      double powerLawFunction(double x);

      std::vector<double> xValues;
      
   public:
      void setUp();
      void tearDown();

      void index();
      void linear();
      void powerlaw();
      void lagrange();
};

#endif
