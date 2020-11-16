#ifndef testSkymap_h
#define testSkymap_h

#include <cppunit/extensions/HelperMacros.h>
#include <Skymap.h>

class testSkymap : public CppUnit::TestFixture
{
   CPPUNIT_TEST_SUITE( testSkymap );
   CPPUNIT_TEST( assignment );
   CPPUNIT_TEST( addition );
   CPPUNIT_TEST( subtraction );
   CPPUNIT_TEST( multiplication );
   CPPUNIT_TEST( division );
   CPPUNIT_TEST_SUITE_END();

   private:
      Skymap<double> NestMapBinned1,
         NestMap1,
         NestMapBinned4,
         NestMap4,
         NestMap10,
         RingMap4;
      
   public:
      void setUp();
      void tearDown();

      void assignment();
      void addition();
      void subtraction();
      void multiplication();
      void division();
};

#endif
