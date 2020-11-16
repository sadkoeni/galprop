#include <testRadialProfiles.h>
#include <radialprofiles.h>
#include <iostream>
#include <stdexcept>
#include <gsl/gsl_interp.h>
#include <Variables.h>
#include <cmath>

CPPUNIT_TEST_SUITE_REGISTRATION( testRadialProfile );

using namespace GalacticStructure;

void testRadialProfile::setUp()
{}

void testRadialProfile::tearDown()
{}

void testRadialProfile::constantProfile()
{

   ConstantRadialProfile pr;

   for (size_t i(0); i < 20; ++i) {
      const double R = i*1.5;
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, pr(R), 1e-8);
   }
}

void testRadialProfile::pulsarProfile()
{
   const double a(1.5), b(3.4), Roff(1.3), R0(6.4);

   PulsarRadialProfile pr(a, b, Roff, R0, "Pfx");

   for (size_t i(0); i < 20; ++i) {
      const double R = i*1.5;
      const double expected = pow(b/a,a)*exp(a-b) * pow((R+Roff)/(R0+Roff),a) * exp(-b*(R-R0)/(R0+Roff));

      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, pr(R), 1e-8);
   }

}

void testRadialProfile::exponentialProfile()
{
  const double R0(7.9), RS(8.5);
  ExponentialRadialProfile pr(R0, RS, "Pfx_R0", "Pfx_RS");
  
  for (size_t i(0); i < 20; ++i) {
    const double R = i*1.5;
    const double expected = exp(-(R-RS)/R0);

    CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, pr(R), 1e-8);
  }
}

void testRadialProfile::freudenreichWarpProfile()
{
   const double Rw(10.4), c1(0.1), c2(0.03), c3(0.003);
   FreudenreichWarpRadialProfile pr( Rw, c1, c2, c3, "Pfx_Rw", "Pfx_c1", "Pfx_c2", "Pfx_c3");

   for (size_t i(0); i < 20; ++i) {
      const double R = i*1.5;
      const double u = R - Rw;
      const double expected = u > 0 ? c1*u + c2*u*u + c3*u*u*u : 0;

      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, pr(R), 1e-8);
   }

   //Make sure setting variables work
   const double Rwn(1.4), c1n(0.3), c2n(0.06), c3n(0.033);
   pr.setRw(Rwn);
   pr.setc1(c1n);
   pr.setc2(c2n);
   pr.setc3(c3n);

   for (size_t i(0); i < 20; ++i) {
      const double R = i*1.5;
      const double u = R - Rwn;
      const double expected = u > 0 ? c1n*u + c2n*u*u + c3n*u*u*u : 0;

      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, pr(R), 1e-8);
   }
}

void testRadialProfile::cutOffPolynomialProfile()
{
   double Rc(3.4);
   std::vector<double> ci(3);
   std::vector<std::string> ciname(3);
   ci[0] = 0.3;
   ci[1] = -0.2;
   ci[2] = 0.4;
   ciname[0] = "Pfx_c00";
   ciname[1] = "Pfx_c01";
   ciname[2] = "Pfx_c02";

   CutOffPolynomialRadialProfile pr(Rc, ci, "Pfx_Rc", ciname);

   for (size_t i(0); i < 20; ++i) {
      const double R = i*1.5;
      const double u = R - Rc;
      const double expected = u > 0 ? ci[0] + ci[1]*u + ci[2]*u*u : ci[0];

      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, pr(R), 1e-8);
   }

   //Try changing variables
   ci[0] = -0.3;
   ci[1] = 0.2;
   ci[2] = -0.33;
   Rc = 7.7;
   
   pr.setRc(Rc);
   pr.setci(ci);

   for (size_t i(0); i < 20; ++i) {
      const double R = i*1.5;
      const double u = R - Rc;
      const double expected = u > 0 ? ci[0] + ci[1]*u + ci[2]*u*u : ci[0];

      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, pr(R), 1e-8);
   }

   //Make sure it throws when trying something stupid
   std::vector<double> ci0(0);
   std::vector<std::string> ciname0(0);

   CPPUNIT_ASSERT_THROW(pr.setci(ci0), std::runtime_error);
   CPPUNIT_ASSERT_THROW(new CutOffPolynomialRadialProfile(Rc, ci0, "Pfx_Rc", ciname), std::runtime_error);
   CPPUNIT_ASSERT_THROW(new CutOffPolynomialRadialProfile(Rc, ci0, "Pfx_Rc", ciname0), std::runtime_error);

}

void testRadialProfile::sech2Profile()
{
   const double R0(5.6);
   Sech2RadialProfile pr(R0, "Pfx_R0");

   for (size_t i(0); i < 20; ++i) {
      const double R = i*1.5;
      const double expected = pow(1./cosh(R/R0),2);

      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, pr(R), 1e-8);
   }
}

void testRadialProfile::gaussianProfile()
{
   const double R0(7.9);
   const double Roff(2.3);
   GaussianRadialProfile pr(R0, Roff, "Pfx_R0", "Pfx_Roff");

   for (size_t i(0); i < 20; ++i) {
      const double R = i*1.5;
      const double expected = exp(-pow(R+Roff,2)/(R0*R0));

      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, pr(R), 1e-8);
   }
}

void testRadialProfile::multipleGaussianProfile()
{
   const double n1(2.9);
   const double R01(7.9);
   const double Roff1(2.3);
   GaussianRadialProfile pr1(R01, Roff1, "Pfx_R01", "Pfx_Roff1");
   const double n2(8.9);
   const double R02(0.9);
   const double Roff2(-3.3);
   GaussianRadialProfile pr2(R02, Roff2, "Pfx_R02", "Pfx_Roff2");

   utl::Parameters pars;
   pars.setParameter("Pfx_n_00", n1);
   pars.setParameter("Pfx_R0_00", R01);
   pars.setParameter("Pfx_Roff_00", Roff1);
   pars.setParameter("Pfx_n_01", n2);
   pars.setParameter("Pfx_R0_01", R02);
   pars.setParameter("Pfx_Roff_01", Roff2);
   MultipleGaussianRadialProfile mg(pars, "Pfx", 2);

   for (size_t i(0); i < 20; ++i) {
      const double R = i*1.5;
      const double expected = n1*pr1(R)+n2*pr2(R);

      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, mg(R), 1e-8);
   }
}

void testRadialProfile::constantCoreProfile()
{
   std::unique_ptr<RadialProfile> g(new GaussianRadialProfile(3.2, 4.5, "R0","Roff"));
   const double Rcore = 4.3;
   ConstantCoreRadialProfile pr( std::move(g), Rcore, "Rcore");

   GaussianRadialProfile opr(3.2, 4.5, "R0","Roff");

   for (size_t i(0); i < 20; ++i) {
      const double R = i*1.5;
      const double expected = R < Rcore ? 1 : opr(R-Rcore);

      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, pr(R), 1e-8);
   }
}

void testRadialProfile::splineProfile()
{
   std::vector<std::pair<std::string,std::string> > varNames(3);
   varNames[0].first = "R0";
   varNames[0].second = "V0";
   varNames[1].first = "R1";
   varNames[1].second = "V1";
   varNames[2].first = "R2";
   varNames[2].second = "V2";

   std::vector<std::pair<double,double> > values(3);
   values[0].first = 0.0;
   values[0].second = 1.0;
   values[1].first = 15.0;
   values[1].second = 4.0;
   values[2].first = 45.0;
   values[2].second = 2.0;

   SplineRadialProfile pr(values, varNames, gsl_interp_linear);

   for (size_t i(0); i < 20; ++i) {
      const double R = i*1.5;
      double expected;
      if (R < values[1].first)
         expected = values[0].second + (values[1].second - values[0].second)/(values[1].first - values[0].first)*(R-values[0].first);
      else
         expected = values[1].second + (values[2].second - values[1].second)/(values[2].first - values[1].first)*(R-values[1].first);
      

      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, pr(R), 1e-8);
   }
}

void testRadialProfile::logsplineProfile()
{
   std::vector<std::pair<std::string,std::string> > varNames(3);
   varNames[0].first = "R0";
   varNames[0].second = "V0";
   varNames[1].first = "R1";
   varNames[1].second = "V1";
   varNames[2].first = "R2";
   varNames[2].second = "V2";

   std::vector<std::pair<double,double> > values(3);
   values[0].first = 0.0;
   values[0].second = 1.0;
   values[1].first = 15.0;
   values[1].second = 40.0;
   values[2].first = 45.0;
   values[2].second = 1.3e-4;

   SplineLogRadialProfile pr(values, varNames, gsl_interp_linear);

   for (size_t i(0); i < 20; ++i) {
      const double R = i*1.5;
      double expected;
      if (R < values[1].first)
         expected = exp( log(values[0].second) + (log(values[1].second) - log(values[0].second))/(values[1].first - values[0].first)*(R-values[0].first));
      else
         expected = exp( log(values[1].second) + (log(values[2].second) - log(values[1].second))/(values[2].first - values[1].first)*(R-values[1].first));

      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, pr(R), 1e-8);
   }
}

void testRadialProfile::expHoleProfile()
{
   const double Rs(8.5), R0(6.7), Rh(2.3), hi(3.4);
   ExpHoleRadialProfile pr(R0, Rs, Rh, hi, "R0", "Rs", "Rh", "hi");

   for (size_t i(0); i < 20; ++i) {
      const double R = i*1.5;
      const double expected = exp(-(R-Rs)/R0)*(1-exp(-pow(R/Rh,hi)));

      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, pr(R), 1e-8);
   }
}

void testRadialProfile::expGaussHoleProfile()
{
   const double sigma = 2.3;
   const double mu = 4.3;
   const double rh = 9.2;

   ExpGaussianHoleRadialProfile pr( sigma, mu, rh, "sigma0", "mu0", "Rh");

   const double rsmooth = sigma*sigma/rh + mu;
   const double f0 = exp(-rsmooth/rh) * exp(0.5*(rsmooth - mu)*(rsmooth - mu)/sigma/sigma);
   for (size_t i(0); i < 20; ++i) {
      const double R = i*1.5;
      const double expected = R < rsmooth ? f0*exp(-0.5*(R - mu)*(R - mu)/sigma/sigma) : exp(-R/rh);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, pr(R), 1e-8);
   }
}

void testRadialProfile::createFromXML()
{
   //Set up the XML
   std::string xmlString = 
      "<someBranch>"
      "<RadialProfile type='Constant'>"
      "<name>ConstantProfile</name>"
      "</RadialProfile>"
      "<RadialProfile type='Pulsar'>"
      "<name>PulsarProfile</name>"
      "<variable id='alpha'> <value>2.0</value> </variable>"
      "<variable id='beta'> <value>5.0</value> </variable>"
      "<variable id='Roff'> <value>1.0</value> </variable>"
      "<variable id='R0'> <value>8.5</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='Exponential'>"
      "<name>ExponentialProfile</name>"
      "<variable id='R0'> <value>2.5</value> </variable>"
      "<variable id='RS'> <value>8.5</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='Sech2'>"
      "<name>Sech2Profile</name>"
      "<variable id='R0'> <value>2.5</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='Gaussian'>"
      "<name>GaussianProfile</name>"
      "<variable id='R0'> <value>2.5</value> </variable>"
      "<variable id='Roff'> <value>3.0</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='CutOffPolynomial'>"
      "<name>CutOffPolynomialProfile</name>"
      "<polynomialDegree>2</polynomialDegree>"
      "<variable id='Rc'> <value>4.5</value> </variable>"
      "<variable id='c00'> <value>0.1</value> </variable>"
      "<variable id='c01'> <value>-0.3</value> </variable>"
      "<variable id='c02'> <value>0.2</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='Spline'>"
      "<name>SplineProfile</name>"
      "<numberOfPoints>5</numberOfPoints>"
      "<interpolationType>CSPLINE</interpolationType>"
      "<variable id='radius_00'> <value>0.0</value> </variable>"
      "<variable id='value_00'> <value>1.0</value> </variable>"
      "<variable id='radius_01'> <value>5.0</value> </variable>"
      "<variable id='value_01'> <value>3.0</value> </variable>"
      "<variable id='radius_02'> <value>10</value> </variable>"
      "<variable id='value_02'> <value>3.0</value> </variable>"
      "<variable id='radius_03'> <value>20</value> </variable>"
      "<variable id='value_03'> <value>1.0</value> </variable>"
      "<variable id='radius_04'> <value>50</value> </variable>"
      "<variable id='value_04'> <value>0.1</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='SplineLog'>"
      "<name>SplineLogProfile</name>"
      "<numberOfPoints>5</numberOfPoints>"
      "<interpolationType>CSPLINE</interpolationType>"
      "<variable id='radius_00'> <value>0.0</value> </variable>"
      "<variable id='value_00'> <value>1.0</value> </variable>"
      "<variable id='radius_01'> <value>5.0</value> </variable>"
      "<variable id='value_01'> <value>3.0</value> </variable>"
      "<variable id='radius_02'> <value>10</value> </variable>"
      "<variable id='value_02'> <value>3.0</value> </variable>"
      "<variable id='radius_03'> <value>20</value> </variable>"
      "<variable id='value_03'> <value>1e-1</value> </variable>"
      "<variable id='radius_04'> <value>50</value> </variable>"
      "<variable id='value_04'> <value>1e-4</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='Constant'>"
      "</RadialProfile>"
      "<RadialProfile type='Pulsar'>"
      "<name>PulsarMissingVariable</name>"
      "<variable id='alpha'> <value>2.0</value> </variable>"
      "<variable id='beta'> <value>5.0</value> </variable>"
      "<variable id='Roff'> <value>1.0</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='Spline'>"
      "<name>SplineWithoutNumberOfPoints</name>"
      "<interpolationType>CSPLINE</interpolationType>"
      "</RadialProfile>"
      "<RadialProfile type='SplineLog'>"
      "<name>SplineLogWihtoutType</name>"
      "<numberOfPoints>5</numberOfPoints>"
      "</RadialProfile>"
      "</someBranch>";

   //And the branch
   utl::ReaderStringInput xmlInput(xmlString);
   utl::Reader xmlReader(xmlInput);
   utl::Branch b = xmlReader.GetTopBranch().GetFirstChild();

   //Loop through the profiles and compare them to manually set up profiles
   std::unique_ptr<RadialProfile> xmlpr = RadialProfile::createProfile(b);
   std::unique_ptr<RadialProfile> cmppr(new ConstantRadialProfile());

   compareProfiles(cmppr, xmlpr);

   b = b.GetNextSibling();
   xmlpr.release();
   cmppr.release();

   xmlpr = RadialProfile::createProfile(b);
   cmppr = std::unique_ptr<RadialProfile>(new PulsarRadialProfile(2.0, 5.0, 1.0, 8.5, "Pfx"));

   compareProfiles(cmppr, xmlpr);

   b = b.GetNextSibling();
   xmlpr.release();
   cmppr.release();

   xmlpr = RadialProfile::createProfile(b);
   cmppr = std::unique_ptr<RadialProfile>(new ExponentialRadialProfile(2.5, 8.5, "Pfx_R0", "Pfx_RS"));

   compareProfiles(cmppr, xmlpr);

   //Make sure the variable name is correct
   std::vector<std::string> varNames = xmlpr->getVariables().getNames();
   const std::string expectedr0 = "ExponentialProfile_R0";
   const std::string expectedrs = "ExponentialProfile_RS";
   CPPUNIT_ASSERT_EQUAL(expectedr0, varNames[0]);
   CPPUNIT_ASSERT_EQUAL(expectedrs, varNames[1]);

   b = b.GetNextSibling();
   xmlpr.release();
   cmppr.release();

   xmlpr = RadialProfile::createProfile(b);
   cmppr = std::unique_ptr<RadialProfile>(new Sech2RadialProfile(2.5, "Pfx_R0"));

   compareProfiles(cmppr, xmlpr);

   b = b.GetNextSibling();
   xmlpr.release();
   cmppr.release();

   xmlpr = RadialProfile::createProfile(b);
   cmppr = std::unique_ptr<RadialProfile>(new GaussianRadialProfile(2.5, 3.0, "Pfx_R0", "Pfx_Roff"));

   compareProfiles(cmppr, xmlpr);

   b = b.GetNextSibling();
   xmlpr.release();
   cmppr.release();

   std::vector<std::string> cnames(3);
   std::vector<double> cis(3);
   cnames[0] = "c0";
   cnames[1] = "c1";
   cnames[2] = "c2";
   cis[0] = 0.1;
   cis[1] = -0.3;
   cis[2] = 0.2;
   xmlpr = RadialProfile::createProfile(b);
   cmppr = std::unique_ptr<RadialProfile>(new CutOffPolynomialRadialProfile(4.5, cis, "Pfx_Rc", cnames));

   compareProfiles(cmppr, xmlpr);

   b = b.GetNextSibling();
   xmlpr.release();
   cmppr.release();

   //Set up the variable and values for the spline
   std::vector<std::pair<std::string, std::string> > names(5);
   std::vector<std::pair<double, double> > values(5);

   names[0].first = "radius_00";
   names[0].second = "values_00";
   names[1].first = "radius_01";
   names[1].second = "values_01";
   names[2].first = "radius_02";
   names[2].second = "values_02";
   names[3].first = "radius_03";
   names[3].second = "values_03";
   names[4].first = "radius_04";
   names[4].second = "values_04";
   values[0].first = 0.0;
   values[0].second = 1.0;
   values[1].first = 5.0;
   values[1].second = 3.0;
   values[2].first = 10.0;
   values[2].second = 3.0;
   values[3].first = 20.0;
   values[3].second = 1.0;
   values[4].first = 50.0;
   values[4].second = 0.1;

   xmlpr = RadialProfile::createProfile(b);
   cmppr = std::unique_ptr<RadialProfile>(new SplineRadialProfile(values, names, gsl_interp_cspline));

   compareProfiles(cmppr, xmlpr);

   b = b.GetNextSibling();
   xmlpr.release();
   cmppr.release();

   //Change values for log spline
   values[3].second = 0.1;
   values[4].second = 1e-4;

   xmlpr = RadialProfile::createProfile(b);
   cmppr = std::unique_ptr<RadialProfile>(new SplineLogRadialProfile(values, names, gsl_interp_cspline));

   compareProfiles(cmppr, xmlpr);

   b = b.GetNextSibling();

   //No name, should be an error
   CPPUNIT_ASSERT_THROW(RadialProfile::createProfile(b), std::runtime_error);

   b = b.GetNextSibling();

   //Missing a variable
   CPPUNIT_ASSERT_THROW(RadialProfile::createProfile(b), utl::Variables::VariableError);

   b = b.GetNextSibling();

   //Spline missing elements
   CPPUNIT_ASSERT_THROW(RadialProfile::createProfile(b), std::runtime_error);

   b = b.GetNextSibling();

   //Spline missing elements
   CPPUNIT_ASSERT_THROW(RadialProfile::createProfile(b), std::runtime_error);
}

void testRadialProfile::compareProfiles(const std::unique_ptr<GalacticStructure::RadialProfile> &expected, const std::unique_ptr<GalacticStructure::RadialProfile> &test)
{
   for (size_t i(0); i < 10; ++i) {
      double R = i*1.5;
      CPPUNIT_ASSERT_DOUBLES_EQUAL( (*expected)(R), (*test)(R), 1e-8);
   }
}
