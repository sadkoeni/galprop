#ifndef testRadialProfile_h
#define testRadialProfile_h

#include <cppunit/extensions/HelperMacros.h>
#include <radialprofiles.h>

class testRadialProfile : public CppUnit::TestFixture
{
   CPPUNIT_TEST_SUITE( testRadialProfile );
   CPPUNIT_TEST( createFromXML );
   CPPUNIT_TEST( constantProfile );
   CPPUNIT_TEST( pulsarProfile );
   CPPUNIT_TEST( exponentialProfile );
   CPPUNIT_TEST( freudenreichWarpProfile );
   CPPUNIT_TEST( cutOffPolynomialProfile );
   CPPUNIT_TEST( expHoleProfile );
   CPPUNIT_TEST( expGaussHoleProfile );
   CPPUNIT_TEST( sech2Profile );
   CPPUNIT_TEST( gaussianProfile );
   CPPUNIT_TEST( multipleGaussianProfile );
   CPPUNIT_TEST( constantCoreProfile );
   CPPUNIT_TEST( splineProfile );
   CPPUNIT_TEST( logsplineProfile );
   CPPUNIT_TEST_SUITE_END();
   
   public:
      void setUp();
      void tearDown();

      void createFromXML();
      void constantProfile();
      void pulsarProfile();
      void exponentialProfile();
      void freudenreichWarpProfile();
      void cutOffPolynomialProfile();
      void expHoleProfile();
      void expGaussHoleProfile();
      void sech2Profile();
      void gaussianProfile();
      void multipleGaussianProfile();
      void constantCoreProfile();
      void splineProfile();
      void logsplineProfile();

   private:
      void compareProfiles(const std::unique_ptr<GalacticStructure::RadialProfile> &expected, const std::unique_ptr<GalacticStructure::RadialProfile> &test);
};

#endif

