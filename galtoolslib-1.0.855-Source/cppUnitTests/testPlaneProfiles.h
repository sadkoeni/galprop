#ifndef testPlaneProfile_h
#define testPlaneProfile_h

#include <cppunit/extensions/HelperMacros.h>
#include <planeprofiles.h>

class testPlaneProfile : public CppUnit::TestFixture
{
   CPPUNIT_TEST_SUITE( testPlaneProfile );
   CPPUNIT_TEST( createFromXML );
   CPPUNIT_TEST_SUITE_END();
   
   public:
      void setUp();
      void tearDown();

      void createFromXML();

   private:
      void compareProfiles(const std::unique_ptr<GalacticStructure::PlaneProfile> &expected, const std::unique_ptr<GalacticStructure::PlaneProfile> &test, const std::string &message);
};

#endif

