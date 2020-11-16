#include <testPlaneProfiles.h>
#include <planeprofiles.h>
#include <spiralarms.h>
#include <iostream>
#include <stdexcept>
#include <Variables.h>
#include <Parameters.h>

CPPUNIT_TEST_SUITE_REGISTRATION( testPlaneProfile );

using namespace GalacticStructure;

void testPlaneProfile::setUp()
{}

void testPlaneProfile::tearDown()
{}

void testPlaneProfile::createFromXML()
{
   //Set up the XML
   const std::string xmlString = 
      "<someBranch>"
      "<PlaneProfile type='1Mode'>"
      "<name>1ModeProfile</name>"
      "<RadialProfile type='Constant'>"
      "<name>ConstantProfile</name>"
      "</RadialProfile>"
      "<variable id='norm'> <value>2.0</value> </variable>"
      "</PlaneProfile>"

      "<PlaneProfile type='2Mode'>"
      "<name>2ModeProfile</name>"
      "<RadialProfile type='Pulsar' intent='norm1'>"
      "<name>PulsarProfile</name>"
      "<variable id='alpha'> <value>2.0</value> </variable>"
      "<variable id='beta'> <value>5.0</value> </variable>"
      "<variable id='Roff'> <value>1.0</value> </variable>"
      "<variable id='R0'> <value>8.5</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='Exponential' intent='norm2'>"
      "<name>ExponentialProfile</name>"
      "<variable id='R0'> <value>2.5</value> </variable>"
      "<variable id='RS'> <value>8.5</value> </variable>"     
      "</RadialProfile>"
      "<RadialProfile type='Sech2' intent='zero2'>"
      "<name>Sech2Profile</name>"
      "<variable id='R0'> <value>2.5</value> </variable>"
      "</RadialProfile>"
      "<variable id='norm1'> <value>2.0</value> </variable>"
      "<variable id='norm2'> <value>2.0</value> </variable>"
      "</PlaneProfile>"

      "<PlaneProfile type='3Mode'>"
      "<name>3ModeProfile</name>"
      "<RadialProfile type='Gaussian' intent='norm1'>"
      "<name>GaussianProfile</name>"
      "<variable id='R0'> <value>2.5</value> </variable>"
      "<variable id='Roff'> <value>3.0</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='Gaussian' intent='norm2'>"
      "<name>GaussianProfile2</name>"
      "<variable id='R0'> <value>2.5</value> </variable>"
      "<variable id='Roff'> <value>3.0</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='Gaussian' intent='norm3'>"
      "<name>GaussianProfile3</name>"
      "<variable id='R0'> <value>2.5</value> </variable>"
      "<variable id='Roff'> <value>3.0</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='Spline' intent='zero2'>"
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
      "<RadialProfile type='SplineLog' intent='zero3'>"
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
      "<variable id='norm1'> <value>2.0</value> </variable>"
      "<variable id='norm2'> <value>2.0</value> </variable>"
      "<variable id='norm3'> <value>2.0</value> </variable>"
      "</PlaneProfile>"

      "<PlaneProfile type='NMode'>"
      "<name>NModeProfile</name>"
      "<numberOfModes>3</numberOfModes>"
      "<RadialProfile type='Gaussian' intent='norm1'>"
      "<name>GaussianProfile</name>"
      "<variable id='R0'> <value>2.5</value> </variable>"
      "<variable id='Roff'> <value>3.0</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='Gaussian' intent='norm2'>"
      "<name>GaussianProfile2</name>"
      "<variable id='R0'> <value>2.5</value> </variable>"
      "<variable id='Roff'> <value>3.0</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='Gaussian' intent='norm3'>"
      "<name>GaussianProfile3</name>"
      "<variable id='R0'> <value>2.5</value> </variable>"
      "<variable id='Roff'> <value>3.0</value> </variable>"
      "</RadialProfile>"
      "<RadialProfile type='Spline' intent='zero2'>"
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
      "<RadialProfile type='SplineLog' intent='zero3'>"
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
      "<variable id='norm1'> <value>2.0</value> </variable>"
      "<variable id='norm2'> <value>2.0</value> </variable>"
      "<variable id='norm3'> <value>2.0</value> </variable>"
      "<variable id='zero2'> <value>2.0</value> </variable>"
      "<variable id='zero3'> <value>2.0</value> </variable>"
      "</PlaneProfile>"

      "<PlaneProfile type='Arm'>"
      "<name>ArmProfile</name>"
      "<armfunction>Exponential</armfunction>"
      "<variable id='norm'> <value>2.0</value> </variable>"
      "<variable id='a'> <value>4.0</value> </variable>"
      "<variable id='phiMin'> <value>1.0</value> </variable>"
      "<variable id='rMin'> <value>2.0</value> </variable>"
      "<variable id='rMax'> <value>20.0</value> </variable>"
      "<variable id='width'> <value>1.0</value> </variable>"
      "<RadialProfile type='Gaussian'>"
      "<name>GaussianProfile</name>"
      "<variable id='R0'> <value>2.5</value> </variable>"
      "<variable id='Roff'> <value>3.0</value> </variable>"
      "</RadialProfile>"
      "</PlaneProfile>"

      "<PlaneProfile type='SymmetricArms'>"
      "<name>SymmetricArmProfile</name>"
      "<numberOfArms>4</numberOfArms>"
      "<armfunction>Sech2</armfunction>"
      "<variable id='norm'> <value>2.0</value> </variable>"
      "<variable id='a'> <value>4.0</value> </variable>"
      "<variable id='phiMin'> <value>1.0</value> </variable>"
      "<variable id='rMin'> <value>2.0</value> </variable>"
      "<variable id='rMax'> <value>20.0</value> </variable>"
      "<variable id='width'> <value>1.0</value> </variable>"
      "<RadialProfile type='Gaussian'>"
      "<name>GaussianProfile</name>"
      "<variable id='R0'> <value>2.5</value> </variable>"
      "<variable id='Roff'> <value>3.0</value> </variable>"
      "</RadialProfile>"
      "</PlaneProfile>"
      "</someBranch>";

   //And the branch
   utl::ReaderStringInput xmlInput(xmlString);
   utl::Reader xmlReader(xmlInput);
   utl::Branch b = xmlReader.GetTopBranch().GetFirstChild();

   //Set up the parameter instance to create the profiles for comparison
   utl::Parameters pars;
   pars.setParameter("norm", "2.0");
   pars.setParameter("norm1", "2.0");
   pars.setParameter("norm2", "2.0");
   pars.setParameter("norm3", "2.0");
   pars.setParameter("zero2", "2.0");
   pars.setParameter("zero3", "2.0");
   pars.setParameter("a", "4.0");
   pars.setParameter("phiMin", "1.0");
   pars.setParameter("rMin", "2.0");
   pars.setParameter("rMax", "20.0");
   pars.setParameter("width", "1.0");

   //Loop through the profiles and compare them to manually set up profiles
   std::unique_ptr<PlaneProfile> xmlpr = PlaneProfile::createProfile(b);
   std::unique_ptr<RadialProfile> crp (new ConstantRadialProfile() );
   std::unique_ptr<PlaneProfile> cmppr (new PlaneProfile1Mode(pars, "norm", std::move(crp)));

   compareProfiles(cmppr, xmlpr, "1Mode profile");
   
   b = b.GetNextSibling();
   xmlpr.release();
   cmppr.release();

   xmlpr = PlaneProfile::createProfile(b);
   std::unique_ptr<RadialProfile> prp( new PulsarRadialProfile(2.0, 5.0, 1.0, 8.5, "Pfx"));
   std::unique_ptr<RadialProfile> erp( new ExponentialRadialProfile(2.5, 8.5, "Pfx_R0", "Pfx_RS"));
   std::unique_ptr<RadialProfile> spr( new Sech2RadialProfile(2.5, "Pfx_R0"));
   cmppr = std::unique_ptr<PlaneProfile>(new PlaneProfile2Mode(pars, "norm1", "norm2", std::move(prp), std::move(erp), std::move(spr)));

   compareProfiles(cmppr, xmlpr, "2Mode profile");
   
   b = b.GetNextSibling();
   xmlpr.release();
   cmppr.release();

   xmlpr = PlaneProfile::createProfile(b);
   std::unique_ptr<RadialProfile> grp( new GaussianRadialProfile(2.5, 3.0, "Pfx_R0", "Pfx_Roff"));
   std::unique_ptr<RadialProfile> grp2( new GaussianRadialProfile(2.5, 3.0, "Pfx_R0", "Pfx_Roff"));
   std::unique_ptr<RadialProfile> grp3( new GaussianRadialProfile(2.5, 3.0, "Pfx_R0", "Pfx_Roff"));
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
   std::unique_ptr<RadialProfile> splpr (new SplineRadialProfile(values, names, gsl_interp_cspline));
   //Change values for log spline
   values[3].second = 0.1;
   values[4].second = 1e-4;
   std::unique_ptr<RadialProfile> lsplpr ( new SplineLogRadialProfile (values, names, gsl_interp_cspline));
   cmppr = std::unique_ptr<PlaneProfile>(new PlaneProfile3Mode(pars, "norm1", "norm2", "norm3", std::move(grp), std::move(grp2), std::move(splpr), std::move(grp3), std::move(lsplpr)));

   compareProfiles(cmppr, xmlpr, "3Mode profile");
   
   b = b.GetNextSibling();
   xmlpr.release();
   cmppr.release();

   xmlpr = PlaneProfile::createProfile(b);
   std::vector<std::string> nnames(3),znames(2);
   nnames[0] = "norm1";
   nnames[1] = "norm2";
   nnames[2] = "norm3";
   znames[0] = "zero2";
   znames[1] = "zero3";
   std::vector<std::unique_ptr<RadialProfile> > nprof;
   std::vector<std::unique_ptr<RadialProfile> > zprof;
   nprof.push_back(std::unique_ptr<RadialProfile>( new GaussianRadialProfile(2.5, 3.0, "Pfx_R0", "Pfx_Roff")));
   nprof.push_back(std::unique_ptr<RadialProfile>( new GaussianRadialProfile(2.5, 3.0, "Pfx_R0", "Pfx_Roff")));
   nprof.push_back(std::unique_ptr<RadialProfile>( new GaussianRadialProfile(2.5, 3.0, "Pfx_R0", "Pfx_Roff")));
   //Values for regular spline
   values[3].second = 1.0;
   values[4].second = 0.1;
   zprof.push_back( std::unique_ptr<RadialProfile>( new SplineRadialProfile(values, names, gsl_interp_cspline)));
   //Change values for log spline
   values[3].second = 0.1;
   values[4].second = 1e-4;
   zprof.push_back(std::unique_ptr<RadialProfile>( new SplineLogRadialProfile (values, names, gsl_interp_cspline)));
   cmppr = std::unique_ptr<PlaneProfile>(new PlaneProfileNMode(pars, nnames, znames, std::move(nprof), std::move(zprof)));

   compareProfiles(cmppr, xmlpr, "NMode profile");
   
   b = b.GetNextSibling();
   xmlpr.release();
   cmppr.release();

   xmlpr = PlaneProfile::createProfile(b);
   nnames.resize(6);
   nnames[0]="norm";
   nnames[1]="a";
   nnames[2]="rMin";
   nnames[3]="phiMin";
   nnames[4]="width";
   nnames[5]="rMax";
   cmppr = std::unique_ptr<PlaneProfile>(new PlaneProfileArm(pars,
            std::unique_ptr<ArmFunction>(new ArmFunction(4.0, 2.0, 1.0, 20.0, 1.0, std::unique_ptr<ScaleFunction>(new ExpFunction))),
            std::unique_ptr<RadialProfile>( new GaussianRadialProfile(2.5, 3.0, "Pfx_R0", "Pfx_Roff")), 
            nnames, true));
   
   compareProfiles(cmppr, xmlpr, "Arm profile");
   
   b = b.GetNextSibling();
   xmlpr.release();
   cmppr.release();

   xmlpr = PlaneProfile::createProfile(b);
   std::vector<std::unique_ptr<ArmFunction> > arms; 
   for (size_t i(0); i < 4; ++i)
      arms.push_back(std::unique_ptr<ArmFunction>(new ArmFunction(4.0, 1.0, 2.0, 20.0, 1.0, std::unique_ptr<ScaleFunction>(new Sech2Function))));
   cmppr = std::unique_ptr<PlaneProfile>(new PlaneProfileSymmetricArms(pars, std::move(arms), 
            std::unique_ptr<RadialProfile>( new GaussianRadialProfile(2.5, 3.0, "Pfx_R0", "Pfx_Roff")), 
            nnames,true));

   compareProfiles(cmppr, xmlpr, "SymmetricArms profile");
   
}

void testPlaneProfile::compareProfiles(const std::unique_ptr<PlaneProfile> &expected, const std::unique_ptr<PlaneProfile> &test, const std::string &message)
{
   for (size_t i(0); i < 10; ++i) {
      const double R = (i+1)*1.5;
      for (size_t j(0); j < 10; ++j) {
         const double th = j*2*M_PI/10.;
         const double exp = (*expected)(R,th);
         const double t = (*test)(R,th);
         CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE( message, exp, t, 1e-8);
      }
   }
}
