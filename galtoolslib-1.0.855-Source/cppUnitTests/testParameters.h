#ifndef testParameters_h
#define testParameters_h

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestAssert.h>
#include <sstream>

//Template for comparing vectors
namespace CppUnit{

   template<typename T>
      struct assertion_traits<std::vector<T> > {

         static bool equal(const std::vector<T>& x, const std::vector<T>& y) {
            return x == y;
         }

         static std::string toString( const std::vector<T>& x ) {
            if (x.size() == 0)
               return "";

            std::ostringstream ost;
            ost << '"' << x[0] << '"';
            for (size_t i(1); i < x.size(); ++i)
               ost << ", " << '"' << x[i] << '"';

            return ost.str();
         }

      };

}


class testParameters : public CppUnit::TestFixture
{
   CPPUNIT_TEST_SUITE( testParameters );
   CPPUNIT_TEST( parser );
   CPPUNIT_TEST( valueConversion );
   CPPUNIT_TEST( comments );
   CPPUNIT_TEST( vectors );
   CPPUNIT_TEST( errors );
   CPPUNIT_TEST( unused );
   CPPUNIT_TEST_SUITE_END();

   private:
      
   public:
      void setUp();
      void tearDown();

      void parser();
      void valueConversion();
      void comments();
      void vectors();
      void errors();
      void unused();
};

#endif
