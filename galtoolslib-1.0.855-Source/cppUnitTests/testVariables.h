#ifndef testVariables_h
#define testVariables_h

#include <cppunit/extensions/HelperMacros.h>

class testVariables : public CppUnit::TestFixture
{
   CPPUNIT_TEST_SUITE( testVariables );
   CPPUNIT_TEST( addFromParameters );
   CPPUNIT_TEST( addFromXML );
   CPPUNIT_TEST( addFromVariables );
   CPPUNIT_TEST( addManually );
   CPPUNIT_TEST( modifyBounds );
   CPPUNIT_TEST( modifyPriors );
   CPPUNIT_TEST_SUITE_END();
   
   public:
      void setUp();
      void tearDown();

      void addFromParameters();
      void addFromXML();
      void addFromVariables();
      void addManually();
      void modifyBounds();
      void modifyPriors();

};

#endif
