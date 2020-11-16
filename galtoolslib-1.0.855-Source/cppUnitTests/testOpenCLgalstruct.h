#ifndef testOpenCLgalstruct_h
#define testOpenCLgalstruct_h

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestAssert.h>
#include <sstream>

#ifdef HAVE_OPENCL
#include <CL/cl2.hpp>
#endif

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


class testOpenCLgalstruct : public CppUnit::TestFixture
{
   CPPUNIT_TEST_SUITE( testOpenCLgalstruct );
   CPPUNIT_TEST( scaleFunctions );
   CPPUNIT_TEST( radialProfiles );
   CPPUNIT_TEST( planeProfiles );
   CPPUNIT_TEST( spiralArms );
   CPPUNIT_TEST( cylindricalProfiles );
   CPPUNIT_TEST_SUITE_END();

   private:
#ifdef HAVE_OPENCL
      cl::Buffer *b_pars, *b_radius, *b_theta, *b_z, *b_results;
      std::vector<cl::Platform> platforms;
      std::vector<cl::Device> devices;
      cl::Context *context;
      cl::CommandQueue *cq;
      size_t n_buffer;
      std::vector<cl_float> radius;
      std::vector<cl_float> theta;
      std::vector<cl_float> z;

      std::vector<cl::Kernel> compileCode(const std::string &code);
#endif
      
   public:
      void setUp();
      void tearDown();

      void scaleFunctions();
      void radialProfiles();
      void planeProfiles();
      void spiralArms();
      void cylindricalProfiles();
};

#endif
