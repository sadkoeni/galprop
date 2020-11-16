#ifndef testFullsky_h
#define testFullsky_h

#include <cppunit/extensions/HelperMacros.h>
#include <FullSky.h>

class testFullsky : public CppUnit::TestFixture
{
   CPPUNIT_TEST_SUITE( testFullsky );
   CPPUNIT_TEST( assignment );
   CPPUNIT_TEST( interpolation );
   CPPUNIT_TEST( iterators );
   CPPUNIT_TEST( basicFunctions );
   CPPUNIT_TEST( applyFunctions );
   CPPUNIT_TEST( coordinateConversion );
   CPPUNIT_TEST_SUITE_END();

   public:
      void setUp();
      void tearDown();

      void basicFunctions();
      void assignment();
      void iterators();
      void interpolation();
      void applyFunctions();
      void coordinateConversion();
};

#endif
