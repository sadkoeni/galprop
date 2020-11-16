#include "testParameters.h"
#include <Parameters.h>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

CPPUNIT_TEST_SUITE_REGISTRATION( testParameters );

using namespace utl;

void testParameters::setUp() {}

void testParameters::tearDown() {}

void testParameters::parser() {

   //Create a stringstream object with our parameters
   std::string parstr("");
   parstr += "Parameter1 = Value1\n";
   parstr += "Parameter2=Value2\n";
   parstr += "Parameter with space = Value with some spaces\n";
   parstr += "Parameter3 = Value with = sign\n";
   parstr += "Parameter4 = Value with spaces at end    \n";
   parstr += "Parameter5\t = \tValue with tabs\t all\t over\t\n";
   parstr += "\t  Parameter6=  Tabs are ignored in front of parameter name\n";
   parstr += "Parameter7 = Value8\n";
   parstr += "Parameter7 = Value7\n";

   std::istringstream is(parstr);

   //parse the parameterts and check values
   utl::Parameters pars(is);

   std::string value, expected;

   //Normal parameter with spaces around the = sign
   expected = "Value1";
   pars.getParameter("Parameter1", value);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //No spaces around = sign
   expected = "Value2";
   pars.getParameter("Parameter2", value);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //Parameters and values with spaces, do not return spaces
   expected = "Value";
   pars.getParameter("Parameter with space", value);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //Parameters and values with spaces, return spaces
   expected = "Value with some spaces";
   pars.getParameter("Parameter with space", value, true);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //The split happens on the first = sign, so parameter names cannot contain = signs but values can.
   expected = "Value with = sign";
   pars.getParameter("Parameter3", value, true);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //Spaces at end of line are ignored
   expected = "Value with spaces at end";
   pars.getParameter("Parameter4", value, true);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //Tabs are ignored at = sign and end, but not in value
   expected = "Value with tabs\t all\t over";
   pars.getParameter("Parameter5", value, true);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //White space in front of parameter is ignored
   expected = "Tabs are ignored in front of parameter name";
   pars.getParameter("Parameter6", value, true);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   //Parameters are overwritten
   expected = "Value7";
   pars.getParameter("Parameter7", value);
   CPPUNIT_ASSERT_EQUAL(expected, value);

}

void testParameters::valueConversion() {

   //Create a stringstream object with our parameters
   std::string parstr("");
   parstr += "Integer = 1\n";
   parstr += "Float = 1.3425e-5\n";
   parstr += "BoolTrue = 1\n";
   parstr += "BoolFalse = 0\n";
   parstr += "VectorCommaInteger = 1,2,3,4,5,6,7,8,,,\n";
   parstr += "VectorCommaFloat = 1.32e9,2.32,3.333,4.12e-4,5.4321e-23,6.00,7.1,8.32132,,,\n";
   parstr += "VectorSpaceInteger = 1 2 3 4 5 6 7 8\n";
   parstr += "VectorStringComma = This,is,a,comma,separated,string,vector\n";
   parstr += "VectorStringSpace = This is a space separated string vector\n";

   std::istringstream is(parstr);

   //parse the parameterts and check values
   utl::Parameters pars(is);

   //Try to read values to their respective types
   int ivalue, iexpected;
   double dvalue, dexpected;
   bool bvalue, bexpected;
   std::vector<int> vivalue, viexpected;
   std::vector<double> vdvalue, vdexpected;
   std::vector<std::string> vsvalue, vsexpected;

   iexpected = 1;
   pars.getParameter("Integer", ivalue);
   CPPUNIT_ASSERT_EQUAL(iexpected, ivalue);

   dexpected = 1.3425e-5;
   pars.getParameter("Float", dvalue);
   CPPUNIT_ASSERT_EQUAL(dexpected, dvalue);

   bexpected = true;
   pars.getParameter("BoolTrue", bvalue);
   CPPUNIT_ASSERT_EQUAL(bexpected, bvalue);

   bexpected = false;
   pars.getParameter("BoolFalse", bvalue);
   CPPUNIT_ASSERT_EQUAL(bexpected, bvalue);

   viexpected.resize(8);
   for (int i(1); i < 9; ++i)
      viexpected[i-1] = i;
   pars.getParameter("VectorCommaInteger", vivalue);
   CPPUNIT_ASSERT_EQUAL(viexpected, vivalue);
   pars.getParameter("VectorSpaceInteger", vivalue);
   CPPUNIT_ASSERT_EQUAL(viexpected, vivalue);

   vdexpected.resize(8);
   vdexpected[0] = 1.32e9;
   vdexpected[1] = 2.32;
   vdexpected[2] = 3.333;
   vdexpected[3] = 4.12e-4;
   vdexpected[4] = 5.4321e-23;
   vdexpected[5] = 6.00;
   vdexpected[6] = 7.1;
   vdexpected[7] = 8.32132;
   pars.getParameter("VectorCommaFloat", vdvalue);
   CPPUNIT_ASSERT_EQUAL(vdexpected, vdvalue);

   vsexpected.resize(7);
   vsexpected[0] = "This";
   vsexpected[1] = "is";
   vsexpected[2] = "a";
   vsexpected[3] = "comma";
   vsexpected[4] = "separated";
   vsexpected[5] = "string";
   vsexpected[6] = "vector";
   pars.getParameter("VectorStringComma", vsvalue);
   CPPUNIT_ASSERT_EQUAL(vsexpected, vsvalue);

   vsexpected[3] = "space";
   pars.getParameter("VectorStringSpace", vsvalue);
   CPPUNIT_ASSERT_EQUAL(vsexpected, vsvalue);

}

void testParameters::comments() {

   //Create a stringstream object with our parameters
   std::string parstr("");
   parstr += "#Whole = comment line\n";
   parstr += "CoMmEnT: =  another comment line\n";
   parstr += "This line # = is only partially commented out\n";
   parstr += "Another line = that is # only partially commented out\n";

   std::istringstream is1(parstr);

   std::cerr<<"Expect to see errors starting with \"Failed parsing line:\""<<std::endl;

   //Test the default comment string (#) first
   utl::Parameters pars1(is1);

   std::string value, expected;

   expected = "another comment line";
   pars1.getParameter("CoMmEnT:", value, true);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   expected = "that is";
   pars1.getParameter("Another line", value, true);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   CPPUNIT_ASSERT_THROW(pars1.getParameter("This line", value), utl::Parameters::ParameterError);

   //Then with the comment string (CoMmEnT:) 
   std::istringstream is2(parstr);
   utl::Parameters pars2(is2, "CoMmEnT:");

   expected = "is only partially commented out";
   pars2.getParameter("This line #", value, true);
   CPPUNIT_ASSERT_EQUAL(expected, value);

   expected = "that is # only partially commented out";
   pars2.getParameter("Another line", value, true);
   CPPUNIT_ASSERT_EQUAL(expected, value);

}

void testParameters::vectors() {

   //Only works with recent libraries
#ifdef USE_NATIVE_REGEX
   //Create a stringstream object with our parameters
   std::string parstr("");
   parstr += "Vector1%4 = 4\n";
   parstr += "Vector2 = 1 2 3 4 5 6 7 8\n";
   parstr += "Vector2%5 = 20\n";

   std::istringstream is(parstr);

   //parse the parameterts and check values
   utl::Parameters pars(is);

   std::vector<int> vivalue, viexpected;

   viexpected.resize(5);
   viexpected[4] = 4;
   pars.getParameter("Vector1", vivalue);
   CPPUNIT_ASSERT_EQUAL(viexpected, vivalue);

   viexpected.resize(8);
   for (int i(1); i < 9; ++i)
      viexpected[i-1] = i;
   viexpected[5] = 20;
   pars.getParameter("Vector2", vivalue);
   CPPUNIT_ASSERT_EQUAL(viexpected, vivalue);

   //Make sure values are updated
   viexpected[4] = 4;
   pars.getParameter("Vector1", vivalue);
   CPPUNIT_ASSERT_EQUAL(viexpected, vivalue);
#endif
}

void testParameters::errors() {

   //Create a stringstream object with our parameters
   std::string parstr("");
   parstr += "Parameter1 = Value1\n";

   std::istringstream is(parstr);

   //parse the parameterts and check values
   utl::Parameters pars(is);

   CPPUNIT_ASSERT_THROW(pars.getParameter("Does not exist", parstr), utl::Parameters::ParameterError);

   int i;
   CPPUNIT_ASSERT_THROW(pars.getParameter("Parameter1", i), utl::Parameters::ParameterError);

}

void testParameters::unused() {
   //Set up a parameter object
   std::string parstr("");
   parstr += "Vector1%4 = 4\n";
   parstr += "Vector2 = 1 2 3 4 5 6 7 8\n";
   parstr += "Vector2%5 = 20\n";
   parstr += "Integer = 1\n";
   parstr += "Float = 1.3425e-5\n";

   std::istringstream is(parstr);
   utl::Parameters pars(is);

   //Make sure all parameters are unused in the beginning
   auto unused = pars.getUnusedParameters();

   CPPUNIT_ASSERT(unused.find("Integer") != unused.end());
   CPPUNIT_ASSERT(unused.find("Float") != unused.end());
   CPPUNIT_ASSERT(unused.find("Vector2") != unused.end());
#ifdef USE_NATIVE_REGEX
   CPPUNIT_ASSERT_EQUAL(unused.size(), size_t(4));
   CPPUNIT_ASSERT(unused.find("Vector1") != unused.end());
   CPPUNIT_ASSERT(unused.find("Vector1%4") == unused.end());
   CPPUNIT_ASSERT(unused.find("Vector2%5") == unused.end());
#else
   CPPUNIT_ASSERT_EQUAL(unused.size(), size_t(5));
   CPPUNIT_ASSERT(unused.find("Vector1") == unused.end());
   CPPUNIT_ASSERT(unused.find("Vector1%4") != unused.end());
   CPPUNIT_ASSERT(unused.find("Vector2%5") != unused.end());
#endif

   //Access one parameter and make sure it got used and nothing else
   double f;
   pars.getParameter("Float", f);

   unused = pars.getUnusedParameters();
   CPPUNIT_ASSERT(unused.find("Integer") != unused.end());
   CPPUNIT_ASSERT(unused.find("Float") == unused.end());
   CPPUNIT_ASSERT(unused.find("Vector2") != unused.end());
#ifdef USE_NATIVE_REGEX
   CPPUNIT_ASSERT(unused.find("Vector1") != unused.end());
   CPPUNIT_ASSERT(unused.find("Vector1%4") == unused.end());
   CPPUNIT_ASSERT(unused.find("Vector2%5") == unused.end());
#else
   CPPUNIT_ASSERT(unused.find("Vector1") == unused.end());
   CPPUNIT_ASSERT(unused.find("Vector1%4") != unused.end());
   CPPUNIT_ASSERT(unused.find("Vector2%5") != unused.end());
#endif

   //Try access through a full string 
   std::string str;
   pars.getParameter("Integer", str, true);

   unused = pars.getUnusedParameters();
   CPPUNIT_ASSERT(unused.find("Integer") == unused.end());

   //Access through a vector
   std::vector<double> vec;
   pars.getParameter("Vector2", vec);

   unused = pars.getUnusedParameters();
   CPPUNIT_ASSERT(unused.find("Vector2") == unused.end());

#ifdef USE_NATIVE_REGEX
   //Test access through vector
   pars.getParameter("Vector1", vec);

   unused = pars.getUnusedParameters();
   CPPUNIT_ASSERT(unused.find("Vector1") == unused.end());
#endif
   
}
