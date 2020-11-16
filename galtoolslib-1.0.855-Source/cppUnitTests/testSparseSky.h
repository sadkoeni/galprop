#ifndef testSparsesky_h
#define testSparsesky_h

#include <cppunit/extensions/HelperMacros.h>
#include <SparseSky.h>

class testSparsesky : public CppUnit::TestFixture
{
   CPPUNIT_TEST_SUITE( testSparsesky );
   CPPUNIT_TEST( assignment );
   CPPUNIT_TEST( basicFunctions );
   CPPUNIT_TEST( interpolation );
   CPPUNIT_TEST( iterators );
   CPPUNIT_TEST( applyFunctions );
   CPPUNIT_TEST( coordinateConversion );
   CPPUNIT_TEST_SUITE_END();

   public:
      void setUp();
      void tearDown();

      void assignment();
      void basicFunctions();
      void iterators();
      void interpolation();
      void applyFunctions();
      void coordinateConversion();
};

#endif
