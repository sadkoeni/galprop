#include "testOpenCLgalstruct.h"
#include <radialprofiles.h>
#include <planeprofiles.h>
#include <cylindricalprofiles.h>
#include <spiralarms.h>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <Registry.h>

//No need to register it unless we have OpenCL available
#ifdef HAVE_OPENCL
CPPUNIT_TEST_SUITE_REGISTRATION( testOpenCLgalstruct );
#endif

using namespace utl;

void testOpenCLgalstruct::setUp() 
{
#ifdef HAVE_OPENCL
   cl_int error;
   cl::Platform::get(&platforms);

   //Perform the tests on the CPU
   platforms[0].getDevices(CL_DEVICE_TYPE_GPU, &devices);
   //Fall back on CPU if needed
   if (devices.size() == 0)
      platforms[0].getDevices(CL_DEVICE_TYPE_CPU, &devices);

   context = new cl::Context(devices);

   //Set up the buffers, 30000 should be enough for all relevant parameters
   n_buffer = 30000;
   b_pars = new cl::Buffer(*context, CL_MEM_READ_ONLY, 2000*sizeof(cl_float), NULL, &error);
   b_radius = new cl::Buffer(*context, CL_MEM_READ_ONLY, n_buffer*sizeof(cl_float), NULL, &error);
   b_theta = new cl::Buffer(*context, CL_MEM_READ_ONLY, n_buffer*sizeof(cl_float), NULL, &error);
   b_z = new cl::Buffer(*context, CL_MEM_READ_ONLY, n_buffer*sizeof(cl_float), NULL, &error);
   b_results = new cl::Buffer(*context, CL_MEM_WRITE_ONLY, n_buffer*sizeof(cl_float), NULL, &error);

   //Set up the points at which to test the profiles
   radius.resize(n_buffer);
   theta.resize(n_buffer);
   z.resize(n_buffer);
   for (size_t i(0); i < n_buffer; ++i){
      radius[i] = double(i)*100.47842314/double(n_buffer);
      theta[i] = double(i%1578)/double(1577) * 2.*M_PI;
      z[i] = double(i%878)/double(877) * 8.-4.;
   }

   //Create a command queue on the device
   cq = new cl::CommandQueue(*context, devices[0]);

   //Write the positions to the buffers
   cq->enqueueWriteBuffer(*b_radius, true, 0, n_buffer*sizeof(cl_float), &radius[0]);
   cq->enqueueWriteBuffer(*b_theta, true, 0, n_buffer*sizeof(cl_float), &theta[0]);
   cq->enqueueWriteBuffer(*b_z, true, 0, n_buffer*sizeof(cl_float), &z[0]);
#endif
}

void testOpenCLgalstruct::tearDown() 
{
#ifdef HAVE_OPENCL
   delete b_pars;
   delete b_radius;
   delete b_theta;
   delete b_z;
   delete b_results;
   delete context;
   delete cq;
#endif
}

#ifdef HAVE_OPENCL
std::vector<cl::Kernel> testOpenCLgalstruct::compileCode(const std::string &code)
{

   //Compile it
   cl::Program::Sources source(std::vector<std::string>(1,code));
   cl::Program program = cl::Program(*context, source);

   cl_int error;
   error = program.build(devices);
   if ( error != CL_SUCCESS) {
      std::cerr<<"Error when building the source "<<error<<std::endl;
      std::string buildLog;
      for (size_t j(0); j < devices.size(); ++j) {
         error = program.getBuildInfo(devices[j], CL_PROGRAM_BUILD_LOG, &buildLog);
         if (error != CL_SUCCESS) {
            std::cerr<<"Error when retrieving the build log "<<error<<std::endl;
            throw(std::runtime_error("OPENCL error"));
         }
         std::cerr<<"Build log for device "<<j<<" while building the test program:"<<std::endl;
         std::cerr<<buildLog<<std::endl;
      }
      throw(std::runtime_error("OPENCL error"));
   }

   std::vector<cl::Kernel> kernels;
   program.createKernels(&kernels);

   return kernels;
}
#endif

void testOpenCLgalstruct::scaleFunctions() {

#ifdef HAVE_OPENCL
   //Use a vector of unique pointers for this
   std::vector< std::unique_ptr<GalacticStructure::ScaleFunction> > functions;

   //Here we can use all registered functions
   std::vector<std::string> names = utl::Registry0<GalacticStructure::ScaleFunction>::GetRegisteredNames();
   for (auto name : names)
      functions.push_back(utl::Registry0<GalacticStructure::ScaleFunction>::create(name));

   //Loop over the profiles and perform the test
   for (size_t i=0; i < functions.size(); ++i) {
      //Get the code and compile it
      std::string fName = "tmp";
      std::string fcode = functions[i]->getOpenCLFunction(fName);

      //Add our kernel to the end of the string
      fcode += "__kernel void run(__global const float *radius, __global float *results)\n"
         "{\n"
         "results[get_global_id(0)] = "+fName+"(radius[get_global_id(0)]);\n"
         "}\n";

      std::vector<cl::Kernel> kernels = compileCode(fcode);

      //Set the arguments to the kernel
      kernels[0].setArg(0, *b_radius);
      kernels[0].setArg(1, *b_results);

      //Run the kernel
      cq->enqueueNDRangeKernel(kernels[0], cl::NullRange, cl::NDRange(radius.size()), cl::NullRange);

      //Read the results
      std::vector<cl_float> openCLresults(radius.size());
      cq->enqueueReadBuffer(*b_results, true, 0, radius.size()*sizeof(cl_float), &openCLresults[0]);

      //Loop over the results and compare to the normal operations
      double error = 1e-30;
      if (names[i] == "Sech1_2")
         error = 1e-19;

      for (size_t j=0; j < radius.size(); ++j) {
         const double expected = (*functions[i])(radius[j]);
         CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(names[i], expected, openCLresults[j], std::max(expected*1e-3,error));
      }
   }
#endif

}

void testOpenCLgalstruct::radialProfiles() {

#ifdef HAVE_OPENCL
   //Use a vector of unique pointers for this
   std::vector< std::unique_ptr<GalacticStructure::RadialProfile> > profiles;

   //Push back the profiles to test
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::ConstantRadialProfile()));
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::PulsarRadialProfile( 1.0, 0.7, 3.0, 6.0, "test" )));
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::ExponentialRadialProfile( 3.4, 8.5, "R0_test", "RS_test" )));
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::ExpHoleRadialProfile( 3.4, 8.5, 1.2, 1.3, "R0_test", "RS_test", "Rh_test", "hi_test" )));
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::ExpGaussianHoleRadialProfile( 0.4, 0.5, 2.2, "sigma0_test", "mu0_test", "Rh_test" )));
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::Sech2RadialProfile( 4.6, "R0")));
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::GaussianRadialProfile( 4.6, 1.1, "R0", "Roff")));
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::ConstantCoreRadialProfile( 
               std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::ExponentialRadialProfile( 1.2, 3.3, "R0", "RS")), 3.31234312, "Rcore")));
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::FreudenreichWarpRadialProfile( 4.5, 0.1, -0.2, 0.4, "Rw", "c1", "c2", "c3")));

   //For the cut off polynomial
   std::vector<double> ci(5);
   std::vector<std::string> ciname(5);
   ci[0] = 0.3;
   ci[1] = 0.1;
   ci[2] = 0.5;
   ci[3] = 0.2;
   ci[4] = 0.4;
   ciname[0] = "c0";
   ciname[1] = "c1";
   ciname[2] = "c2";
   ciname[3] = "c3";
   ciname[4] = "c4";
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::CutOffPolynomialRadialProfile( 6.7, ci, "Rc", ciname)));
   

   //For the MultipleGaussians
   utl::Parameters pars;
   pars.setParameter("Pfx_n_00", 0.42);
   pars.setParameter("Pfx_R0_00", 5.6);
   pars.setParameter("Pfx_Roff_00", 3.4);
   pars.setParameter("Pfx_n_01", 14.3);
   pars.setParameter("Pfx_R0_01", 0.05);
   pars.setParameter("Pfx_Roff_01", -7.8);
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::MultipleGaussianRadialProfile( pars, "Pfx", 2 )));

   //Points for the spline interpolation
   std::vector<std::pair<double, double> > points;
   points.push_back(std::make_pair(2.34, 3.4));
   points.push_back(std::make_pair(8.34, 83.4));
   points.push_back(std::make_pair(12.4, 2.1));
   points.push_back(std::make_pair(13.6, 7.4));
   points.push_back(std::make_pair(25.7, 53.4));
   points.push_back(std::make_pair(82.8, 34.1));
   std::vector<std::pair<std::string,std::string> > vnames;
   for (size_t i=0; i < points.size(); ++i) {
      std::ostringstream ost;
      ost<<i;
      vnames.push_back(std::make_pair("r"+ost.str(),"v"+ost.str()));
   }
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::SplineRadialProfile( points, vnames, gsl_interp_linear)));
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::SplineRadialProfile( points, vnames, gsl_interp_cspline)));
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::SplineRadialProfile( points, vnames, gsl_interp_akima)));
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::SplineLogRadialProfile( points, vnames, gsl_interp_linear)));
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::SplineLogRadialProfile( points, vnames, gsl_interp_cspline)));
   profiles.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::SplineLogRadialProfile( points, vnames, gsl_interp_akima)));

   //Loop over the profiles and perform the test
   for (size_t i=0; i < profiles.size(); ++i) {
      //Get the code and compile it
      std::string fName = "tmp";
      std::string fcode = profiles[i]->getOpenCLFunction(fName, 0);

      //Add our kernel to the end of the string
      fcode += "__kernel void run(__global const float *radius, __global float *results, __constant float *pars)\n"
         "{\n"
         "results[get_global_id(0)] = "+fName+"(radius[get_global_id(0)], pars);\n"
         "}\n";

      std::vector<cl::Kernel> kernels = compileCode(fcode);

      //Set the arguments to the kernel
      kernels[0].setArg(0, *b_radius);
      kernels[0].setArg(1, *b_results);
      kernels[0].setArg(2, *b_pars);

      //Modify the parameters by some small value
      for (size_t k=0; k < 3; ++k) {

         if (k > 0) {
            auto var = profiles[i]->getVariables();

            for (const auto &n : var.getNames() ) {
               var[n] = 1.04*var[n] + 0.04;
            }

            profiles[i]->updateVariableValues(var);
         }

         //Get the parameters and write them to the buffer
         std::vector<cl_float> pars = profiles[i]->getOpenCLPars();

         cq->enqueueWriteBuffer(*b_pars, true, 0, pars.size()*sizeof(cl_float), &pars[0]);

         //Run the kernel
         cq->enqueueNDRangeKernel(kernels[0], cl::NullRange, cl::NDRange(radius.size()), cl::NullRange);

         //Read the results
         std::vector<cl_float> openCLresults(radius.size());
         cq->enqueueReadBuffer(*b_results, true, 0, radius.size()*sizeof(cl_float), &openCLresults[0]);

         //Loop over the results and compare to the normal operations
         std::ostringstream buf;
         buf<<"Profile "<<i<<" with variables "<<k<<std::endl;
         for (size_t j=0; j < radius.size(); ++j) {
            const double expected = (*profiles[i])(radius[j]);
            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(buf.str(), expected, openCLresults[j], std::max(fabs(expected)*1e-3, 1e-35));
         }

      }
   }
#endif

}

void testOpenCLgalstruct::planeProfiles() {

#ifdef HAVE_OPENCL
   //Use a vector of unique pointers for this
   std::vector< std::unique_ptr<GalacticStructure::PlaneProfile> > profiles;

   //Push back the profiles to test
   utl::Parameters pars;
   pars.setParameter("norm", "1.0");
   profiles.push_back(std::unique_ptr<GalacticStructure::PlaneProfile>(
            new GalacticStructure::PlaneProfile1Mode(pars, "norm", 
               std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::PulsarRadialProfile( 1.0, 0.7, 3.0, 6.0, "test" )) )
            ));

   //Points for the spline interpolation
   std::vector<std::pair<double, double> > points;
   points.push_back(std::make_pair(2.34, 3.4));
   points.push_back(std::make_pair(8.34, 83.4));
   points.push_back(std::make_pair(12.4, 2.1));
   points.push_back(std::make_pair(13.6, 7.4));
   points.push_back(std::make_pair(25.7, 53.4));
   points.push_back(std::make_pair(82.8, 34.1));
   std::vector<std::pair<std::string,std::string> > vnames;
   for (size_t i=0; i < points.size(); ++i) {
      std::ostringstream ost;
      ost<<i;
      vnames.push_back(std::make_pair("r"+ost.str(),"v"+ost.str()));
   }
   profiles.push_back(std::unique_ptr<GalacticStructure::PlaneProfile>(
            new GalacticStructure::PlaneProfile1Mode(pars, "norm", 
               std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::SplineRadialProfile( points, vnames, gsl_interp_cspline)) )
            ));

   //Create profiles for the NMode plane profile, use 3 modes
   std::vector< std::unique_ptr<GalacticStructure::RadialProfile> > zeroProfs, normProfs;
   normProfs.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::PulsarRadialProfile( 1.0, 0.7, 3.0, 6.0, "test" )));
   normProfs.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::ExponentialRadialProfile( 3.4, 8.5, "R0_test", "RS_test" )));
   normProfs.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::ExpHoleRadialProfile( 3.4, 8.5, 1.2, 1.3, "ehR0_test", "ehRS_test", "Rh_test", "hi_test" )));
   zeroProfs.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::Sech2RadialProfile( 15.6, "sR0")));
   zeroProfs.push_back(std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::GaussianRadialProfile( 5.6, -10.1, "gR0", "gRoff")));
   //The normalization constants
   std::vector<std::string> normNames(normProfs.size());
   std::vector<std::string> zeroNames(zeroProfs.size());
   for (size_t i(0); i < normProfs.size(); ++i) {
      std::ostringstream buf;
      buf<<"norm"<<i;
      normNames[i] = buf.str();
      buf.str("");
      buf<<4.1323/(i+1);
      pars.setParameter(normNames[i], buf.str());
   }
   for (size_t i(0); i < zeroProfs.size(); ++i) {
      std::ostringstream buf;
      buf<<"zero"<<i;
      zeroNames[i] = buf.str();
      buf.str("");
      buf<<1.1323*i;
      pars.setParameter(zeroNames[i], buf.str());
   }
   profiles.push_back(std::unique_ptr<GalacticStructure::PlaneProfile>(
            new GalacticStructure::PlaneProfileNMode(pars, normNames, zeroNames, std::move(normProfs), std::move(zeroProfs) )
            ));

   //Arm profile
   pars.setParameter("arm_norm", "4.3");
   pars.setParameter("arm_a", "9.3");
   pars.setParameter("arm_rMin", "1.3");
   pars.setParameter("arm_phiMin", "3.3");
   pars.setParameter("arm_rMax", "40.3");
   pars.setParameter("arm_width", "0.3");
   profiles.push_back(std::unique_ptr<GalacticStructure::PlaneProfile>(
            new GalacticStructure::PlaneProfileArm(pars,
               std::unique_ptr<GalacticStructure::ArmFunction>(
                  new GalacticStructure::ArmFunction(5.4, 2.0, 0.3, 4*M_PI, 0.4, 
                     utl::Registry0<GalacticStructure::ScaleFunction>::create("Sech2"))
                  ),
               std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::SplineRadialProfile( points, vnames, gsl_interp_cspline)),
               "arm", true)
            ));


   //Symmetric arm profile
   std::vector< std::unique_ptr<GalacticStructure::ArmFunction> > armFuncs;
   for (size_t i(0); i < 4; ++i)
      armFuncs.push_back( std::unique_ptr<GalacticStructure::ArmFunction> (
               new GalacticStructure::ArmFunction(5.4, 2.0, 0.3, 4*M_PI, 0.4, 
                  utl::Registry0<GalacticStructure::ScaleFunction>::create("Sech2"))
               ));

   pars.setParameter("arm_norm", "4.3");
   pars.setParameter("arm_a", "9.3");
   pars.setParameter("arm_rMin", "1.3");
   pars.setParameter("arm_phiMin", "3.3");
   pars.setParameter("arm_rMax", "40.3");
   pars.setParameter("arm_width", "0.3");
   profiles.push_back(std::unique_ptr<GalacticStructure::PlaneProfile>(
            new GalacticStructure::PlaneProfileSymmetricArms(pars,
               std::move(armFuncs),
               std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::SplineRadialProfile( points, vnames, gsl_interp_cspline)),
               "arm", true)
            ));


   //Loop over the profiles and perform the test
   for (size_t i=0; i < profiles.size(); ++i) {
      //Get the code and compile it
      std::string fName = "tmp";
      std::string fcode = profiles[i]->getOpenCLFunction(fName, 0);

      //Add our kernel to the end of the string
      fcode += "__kernel void run(__global const float *radius, __global const float *theta, __global float *results, __constant float *pars)\n"
         "{\n"
         "results[get_global_id(0)] = "+fName+"(radius[get_global_id(0)], theta[get_global_id(0)], pars);\n"
         "}\n";

      std::vector<cl::Kernel> kernels = compileCode(fcode);

      //Set the arguments to the kernel
      kernels[0].setArg(0, *b_radius);
      kernels[0].setArg(1, *b_theta);
      kernels[0].setArg(2, *b_results);
      kernels[0].setArg(3, *b_pars);

      //Modify the parameters by some small value
      for (size_t k=0; k < 3; ++k) {

         if (k > 0) {
            auto var = profiles[i]->getVariables();

            for (const auto &n : var.getNames() ) {
               var[n] = 1.012432*var[n] + 0.0412341;
            }

            profiles[i]->updateVariableValues(var);
         }

         //Get the parameters and write them to the buffer
         std::vector<cl_float> pars = profiles[i]->getOpenCLPars();

         cq->enqueueWriteBuffer(*b_pars, true, 0, pars.size()*sizeof(cl_float), &pars[0]);

         //Run the kernel
         cq->enqueueNDRangeKernel(kernels[0], cl::NullRange, cl::NDRange(radius.size()), cl::NullRange);

         //Read the results
         std::vector<cl_float> openCLresults(radius.size());
         cq->enqueueReadBuffer(*b_results, true, 0, radius.size()*sizeof(cl_float), &openCLresults[0]);

         //Loop over the results and compare to the normal operations
         for (size_t j=0; j < radius.size(); ++j) {
            std::ostringstream buf;
            buf<<"Profile "<<i<<" at index "<<j<<" ("<<radius[j]<<", "<<theta[j]<<")";
            const double expected = (*profiles[i])(radius[j], theta[j]);
            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(buf.str(),expected, openCLresults[j], std::max(5e-3*fabs(expected), 1e-30));
         }
      }
   }


#endif

}


void testOpenCLgalstruct::spiralArms() {

#ifdef HAVE_OPENCL
   //Use a vector of unique pointers for this
   std::vector< std::unique_ptr<GalacticStructure::ArmFunction> > profiles;

   //Push back the profiles to test
   profiles.push_back(std::unique_ptr<GalacticStructure::ArmFunction>(
            new GalacticStructure::ArmFunction(5.4, 2.0, 0.3, 4*M_PI, 0.4, 
               utl::Registry0<GalacticStructure::ScaleFunction>::create("Sech2"))
            ));

   profiles.push_back(std::unique_ptr<GalacticStructure::ArmFunction>(
            new GalacticStructure::ArmFunction(4.4, 2.2, 4.3, 1*M_PI, 0.8, 
               utl::Registry0<GalacticStructure::ScaleFunction>::create("Gaussian"))
            ));

   profiles.push_back(std::unique_ptr<GalacticStructure::ArmFunction>(
            new GalacticStructure::ArmFunction(10.4, 7.0, 3.2, 0.2*M_PI, 1.4, 
               utl::Registry0<GalacticStructure::ScaleFunction>::create("Step"))
            ));


   //Loop over the profiles and perform the test
   for (size_t i=0; i < profiles.size(); ++i) {
      //Get the code and compile it
      std::string fName = "tmp";
      std::string fcode = profiles[i]->getOpenCLFunction(fName, 0);

      //Add our kernel to the end of the string
      fcode += "__kernel void run(__global const float *radius, __global const float *theta, __global float *results, __constant float *pars)\n"
         "{\n"
         "results[get_global_id(0)] = "+fName+"(radius[get_global_id(0)], theta[get_global_id(0)], pars);\n"
         "}\n";

      std::vector<cl::Kernel> kernels = compileCode(fcode);

      //Set the arguments to the kernel
      kernels[0].setArg(0, *b_radius);
      kernels[0].setArg(1, *b_theta);
      kernels[0].setArg(2, *b_results);
      kernels[0].setArg(3, *b_pars);

      //Get the parameters and write them to the buffer
      std::vector<cl_float> pars = profiles[i]->getOpenCLPars();

      cq->enqueueWriteBuffer(*b_pars, true, 0, pars.size()*sizeof(cl_float), &pars[0]);

      //Run the kernel
      cq->enqueueNDRangeKernel(kernels[0], cl::NullRange, cl::NDRange(radius.size()), cl::NullRange);

      //Read the results
      std::vector<cl_float> openCLresults(radius.size());
      cq->enqueueReadBuffer(*b_results, true, 0, radius.size()*sizeof(cl_float), &openCLresults[0]);

      //Loop over the results and compare to the normal operations
      for (size_t j=0; j < radius.size(); ++j) {
         std::ostringstream buf;
         buf<<"Profile "<<i<<" at index "<<j;
         const double expected = (*profiles[i])(radius[j], theta[j]);
         CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(buf.str(), expected, openCLresults[j], std::max(1e-3*fabs(expected), 1e-30));
      }
   }
#endif

}


void testOpenCLgalstruct::cylindricalProfiles() {

#ifdef HAVE_OPENCL
   //Use a vector of unique pointers for this
   std::vector< std::unique_ptr<GalacticStructure::CylindricalProfile> > profiles;

   //Push back the profiles to test
   utl::Parameters pars;
   pars.setParameter("density", "1.0");
   pars.setParameter("center", "0.0");
   pars.setParameter("northScale", "3.0");
   pars.setParameter("southScale", "1.5");
   profiles.push_back(
         std::unique_ptr<GalacticStructure::CylindricalProfile> (
            new GalacticStructure::GenericDiskProfile ( pars, 
               std::unique_ptr<GalacticStructure::PlaneProfile>(
                  new GalacticStructure::PlaneProfile1Mode(pars, "density", 
                     std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::PulsarRadialProfile( 1.0, 0.7, 3.0, 6.0, "test" )) 
                     )
                  ),
               std::unique_ptr<GalacticStructure::PlaneProfile>(
                  new GalacticStructure::PlaneProfile1Mode(pars, "center", 
                     std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::ExponentialRadialProfile( 2.5, 0.7, "R0test", "RStest" )) 
                     )
                  ),
               std::unique_ptr<GalacticStructure::PlaneProfile>(
                  new GalacticStructure::PlaneProfile1Mode(pars, "northScale", 
                     std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::PulsarRadialProfile( 2.0, 3.1, 1.4, 5.3, "test2" )) 
                     )
                  ),
               std::unique_ptr<GalacticStructure::PlaneProfile>(
                  new GalacticStructure::PlaneProfile1Mode(pars, "southScale", 
                     std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::GaussianRadialProfile( 5.4, 1.1, "g1", "g2" )) 
                     )
                  ),
               utl::Registry0<GalacticStructure::ScaleFunction>::create("Sech2"),
               utl::Registry0<GalacticStructure::ScaleFunction>::create("Exponential")
               )
            )
         );

   //Points for the spline interpolation
   std::vector<std::pair<double, double> > points;
   points.push_back(std::make_pair(2.34, 3.4));
   points.push_back(std::make_pair(8.34, 83.4));
   points.push_back(std::make_pair(12.4, 2.1));
   points.push_back(std::make_pair(13.6, 7.4));
   points.push_back(std::make_pair(25.7, 53.4));
   points.push_back(std::make_pair(82.8, 34.1));
   std::vector<std::pair<std::string,std::string> > vnames;
   for (size_t i=0; i < points.size(); ++i) {
      std::ostringstream ost;
      ost<<i;
      vnames.push_back(std::make_pair("r"+ost.str(),"v"+ost.str()));
   }

   std::vector<std::pair<double, double> > points1;
   points1.push_back(std::make_pair(2.34, 1.4));
   points1.push_back(std::make_pair(8.34, -3.4));
   points1.push_back(std::make_pair(12.4, 2.1));
   points1.push_back(std::make_pair(13.6, 1.4));
   points1.push_back(std::make_pair(25.7, -0.4));
   points1.push_back(std::make_pair(82.8, 4.1));
   std::vector<std::pair<std::string,std::string> > vnames1;
   for (size_t i=0; i < points1.size(); ++i) {
      std::ostringstream ost;
      ost<<i;
      vnames1.push_back(std::make_pair("r1"+ost.str(),"v1"+ost.str()));
   }

   std::vector<std::pair<double, double> > points2;
   points2.push_back(std::make_pair(2.34, 1.4));
   points2.push_back(std::make_pair(8.34, 3.4));
   points2.push_back(std::make_pair(12.4, 2.1));
   points2.push_back(std::make_pair(13.6, 1.4));
   points2.push_back(std::make_pair(25.7, 0.4));
   points2.push_back(std::make_pair(82.8, 4.1));
   std::vector<std::pair<std::string,std::string> > vnames2;
   for (size_t i=0; i < points2.size(); ++i) {
      std::ostringstream ost;
      ost<<i;
      vnames2.push_back(std::make_pair("r2"+ost.str(),"v2"+ost.str()));
   }
   std::vector<std::pair<std::string,std::string> > vnames3;
   for (size_t i=0; i < points2.size(); ++i) {
      std::ostringstream ost;
      ost<<i;
      vnames3.push_back(std::make_pair("r3"+ost.str(),"v3"+ost.str()));
   }

   profiles.push_back(
         std::unique_ptr<GalacticStructure::CylindricalProfile> (
            new GalacticStructure::GenericDiskProfile ( pars, 
               std::unique_ptr<GalacticStructure::PlaneProfile>(
                  new GalacticStructure::PlaneProfile1Mode(pars, "density", 
                     std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::SplineLogRadialProfile( points, vnames, gsl_interp_cspline)) 
                     )
                  ),
               std::unique_ptr<GalacticStructure::PlaneProfile>(
                  new GalacticStructure::PlaneProfile1Mode(pars, "center", 
                     std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::SplineRadialProfile( points1, vnames1, gsl_interp_linear)) 
                     )
                  ),
               std::unique_ptr<GalacticStructure::PlaneProfile>(
                  new GalacticStructure::PlaneProfile1Mode(pars, "northScale", 
                     std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::SplineLogRadialProfile( points2, vnames2, gsl_interp_linear)) 
                     )
                  ),
               std::unique_ptr<GalacticStructure::PlaneProfile>(
                  new GalacticStructure::PlaneProfile1Mode(pars, "southScale", 
                     std::unique_ptr<GalacticStructure::RadialProfile>(new GalacticStructure::SplineLogRadialProfile( points2, vnames3, gsl_interp_linear)) 
                     )
                  ),
               utl::Registry0<GalacticStructure::ScaleFunction>::create("Gaussian"),
               utl::Registry0<GalacticStructure::ScaleFunction>::create("Gaussian")
               )
            )
         );

   //Test the FreudenreichBar profile
   profiles.push_back(
         std::unique_ptr<GalacticStructure::CylindricalProfile> (
            new GalacticStructure::FreudenreichBarProfile( 10, 5, 1.0, 0.5, 0.5, 2.1, 1.5, 3.0, 1.0, 2.0, "phi", "pitch", "barX", "barY", "barZ", "barPerp", "barPara", "Rend", "Hend", "norm")
            )
         );

   //Test the FerriereCMZ profile
   profiles.push_back(
         std::unique_ptr<GalacticStructure::CylindricalProfile> (
            new GalacticStructure::FerriereCMZProfile( 0.1, 0.4, 0.5, 0.005, 0.003, 2.1, "phi", "Xmax", "Hc", "x0", "y0", "norm")
            )
         );

   profiles.push_back(
         std::unique_ptr<GalacticStructure::CylindricalProfile> (
            new GalacticStructure::FerriereDiskProfile( 0.1, 0.4, 0.5, 4.0, 2.0, 0.5, 40, "alpha", "beta", "phid", "Xmax", "Xmin", "Hd", "norm")
            )
         );

   //Test the GenericBar profile
   profiles.push_back(
         std::unique_ptr<GalacticStructure::CylindricalProfile> (
            new GalacticStructure::GenericBarProfile( 0.6, 0.01, 0.5, 2.0, 0.05, 4.0, 1.5, 0.05, 0.01, 3.1, -1.5, 40, "phi", "alpha", "ySc", "rSc", "zSc", "ri", "zi", "x0", "z0", "ei", "pi", "norm")
            )
         );
   

   //Loop over the profiles and perform the test
   for (size_t i=0; i < profiles.size(); ++i) {
      //Get the code and compile it
      std::string fName = "tmp";
      std::string fcode = profiles[i]->getOpenCLFunction(fName, 0);

      //Add our kernel to the end of the string
      fcode += "__kernel void run(__global const float *radius, __global const float *theta, __global const float *z, __global float *results, __constant float *pars)\n"
         "{\n"
         "results[get_global_id(0)] = "+fName+"(radius[get_global_id(0)], theta[get_global_id(0)], z[get_global_id(0)], pars);\n"
         "}\n";

      std::vector<cl::Kernel> kernels = compileCode(fcode);

      //Set the arguments to the kernel
      kernels[0].setArg(0, *b_radius);
      kernels[0].setArg(1, *b_theta);
      kernels[0].setArg(2, *b_z);
      kernels[0].setArg(3, *b_results);
      kernels[0].setArg(4, *b_pars);

      //Modify the parameters by some small value
      for (size_t k=0; k < 3; ++k) {

         if (k > 0) {
            auto var = profiles[i]->getVariables();

            for (const auto &n : var.getNames() ) {
               var[n] = 1.012432*var[n] + 0.0412341;
            }

            profiles[i]->updateVariableValues(var);
         }

         //Get the parameters and write them to the buffer
         std::vector<cl_float> pars = profiles[i]->getOpenCLPars();

         cq->enqueueWriteBuffer(*b_pars, false, 0, pars.size()*sizeof(cl_float), &pars[0]);

         //Run the kernel
         cq->enqueueNDRangeKernel(kernels[0], cl::NullRange, cl::NDRange(radius.size()), cl::NullRange);

         //Read the results
         std::vector<cl_float> openCLresults(radius.size());
         cq->enqueueReadBuffer(*b_results, true, 0, radius.size()*sizeof(cl_float), &openCLresults[0]);

         //Loop over the results and compare to the normal operations
         for (size_t j=0; j < radius.size(); ++j) {
            std::ostringstream buf;
            buf<<"Profile "<<i<<" at index "<<j<<" ("<<radius[j]<<", "<<theta[j]<<", "<<z[j]<<")";
            const double expected = (*profiles[i])(radius[j], theta[j], z[j]);
            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(buf.str(), expected, openCLresults[j], std::max(1e-3*fabs(expected), 1e-30));
         }
      }
   }
#endif

}
