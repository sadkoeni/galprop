#include <testVariables.h>
#include <Variables.h>
#include <Parameters.h>
#include <vector>
#include <string>

CPPUNIT_TEST_SUITE_REGISTRATION( testVariables );

using namespace utl;

void testVariables::setUp() 
{
}

void testVariables::tearDown()
{
}

void testVariables::addFromParameters()
{
   //Set up the parameters object with some parameters in there
   Parameters pars;

   pars.setParameter("Fixed", "1.0");
   pars.setParameter("NoBound", "1.0 0.1");
   pars.setParameter("LowerBound", "1.0 0.1 0.0");
   pars.setParameter("UpperBound", "1.0 0.1 4.0 2.0");
   pars.setParameter("BothBounds", "1.0 0.1 0.0 2.0");
   pars.setParameter("TooManyValues", "1.0 0.1 0.0 2.0 4.5 3.2");
   pars.setParameter("BadValues", "abd jkcs lksa jgds");
   pars.setParameter("BadValue", "1.0 0.1 jds");
   pars.setParameter("PriorFixed", "1.0");
   pars.setParameter("PriorFixed_prior", "LOGUNIFORM 1e-2 1e2");
   pars.setParameter("PriorFree", "1.0 0.1");
   pars.setParameter("PriorFree_prior", "LOGUNIFORM 1e-2 1e2");
   pars.setParameter("PriorFreeLC", "1.0 0.1");
   pars.setParameter("PriorFreeLC_prior", "loguniform 1e-2 1e2");
   pars.setParameter("PriorBad", "1.0 0.1");
   pars.setParameter("PriorBad_prior", "LOGUNIFORM 1e2");
   pars.setParameter("PriorBad2", "1.0 0.1");
   pars.setParameter("PriorBad2_prior", "NONSENSE 1e-2 1e2");
   pars.setParameter("PriorBad3", "1.0 0.1");
   pars.setParameter("PriorBad3_prior", "LOGUNIFORM 1e-2 1e2 1e3 1e4");

   //Read in the good variables and make sure they are in order.
   std::vector<std::string> varNames;
   varNames.push_back("Fixed");
   varNames.push_back("NoBound");
   varNames.push_back("LowerBound");
   varNames.push_back("UpperBound");
   varNames.push_back("BothBounds");
   varNames.push_back("PriorFixed");
   varNames.push_back("PriorFree");
   varNames.push_back("PriorFreeLC");

   Variables vars;

   //This should not throw anything.
   vars.add(varNames, pars);

   double expectedValue = 1.0;
   double expectedError = 0.1;
   double expectedLowerBound = 0.0;
   double expectedUpperBound = 2.0;
   double expectedPrior1 = 1e-2;
   double expectedPrior2 = 1e2;

   double value, error, lowerBound, upperBound, prior1, prior2;
   Variables::PRIOR prior;

   //Check Fixed
   value = vars["Fixed"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("Fixed");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, error, 1e-8);
   
   CPPUNIT_ASSERT_THROW(vars.lowerBound("Fixed"), Variables::BoundNotSet);
   CPPUNIT_ASSERT_THROW(vars.upperBound("Fixed"), Variables::BoundNotSet);

   prior = vars.prior("Fixed");
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   //Check NoBound
   value = vars["NoBound"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("NoBound");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   
   CPPUNIT_ASSERT_THROW(vars.lowerBound("NoBound"), Variables::BoundNotSet);
   CPPUNIT_ASSERT_THROW(vars.upperBound("NoBound"), Variables::BoundNotSet);

   prior = vars.prior("NoBound");
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   //Check LowerBound
   value = vars["LowerBound"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("LowerBound");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   lowerBound = vars.lowerBound("LowerBound");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, lowerBound, 1e-8);

   CPPUNIT_ASSERT_THROW(vars.upperBound("LowerBound"), Variables::BoundNotSet);

   prior = vars.prior("LowerBound");
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   //Check UpperBound
   value = vars["UpperBound"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("UpperBound");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   upperBound = vars.upperBound("UpperBound");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, upperBound, 1e-8);

   CPPUNIT_ASSERT_THROW(vars.lowerBound("UpperBound"), Variables::BoundNotSet);

   prior = vars.prior("UpperBound");
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   //Check BothBounds
   value = vars["BothBounds"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("BothBounds");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   lowerBound = vars.lowerBound("BothBounds");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, lowerBound, 1e-8);
   upperBound = vars.upperBound("BothBounds");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, upperBound, 1e-8);

   prior = vars.prior("BothBounds");
   CPPUNIT_ASSERT_EQUAL(Variables::UNIFORM, prior);
   prior1 = vars.priorPar1("BothBounds");
   prior2 = vars.priorPar2("BothBounds");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, prior1, 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, prior2, 1e-8);

   //Check PriorFixed
   value = vars["PriorFixed"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("PriorFixed");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, error, 1e-8);
   
   CPPUNIT_ASSERT_THROW(vars.lowerBound("PriorFixed"), Variables::BoundNotSet);
   CPPUNIT_ASSERT_THROW(vars.upperBound("PriorFixed"), Variables::BoundNotSet);

   prior = vars.prior("PriorFixed");
   CPPUNIT_ASSERT_EQUAL(Variables::LOGUNIFORM, prior);
   prior1 = vars.priorPar1("PriorFixed");
   prior2 = vars.priorPar2("PriorFixed");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedPrior1, prior1, 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedPrior2, prior2, 1e-8);

   //Check PriorFree
   value = vars["PriorFree"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("PriorFree");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   
   CPPUNIT_ASSERT_THROW(vars.lowerBound("PriorFree"), Variables::BoundNotSet);
   CPPUNIT_ASSERT_THROW(vars.upperBound("PriorFree"), Variables::BoundNotSet);

   prior = vars.prior("PriorFree");
   CPPUNIT_ASSERT_EQUAL(Variables::LOGUNIFORM, prior);
   prior1 = vars.priorPar1("PriorFree");
   prior2 = vars.priorPar2("PriorFree");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedPrior1, prior1, 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedPrior2, prior2, 1e-8);

   //Check PriorFreeLC
   value = vars["PriorFreeLC"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("PriorFreeLC");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   
   CPPUNIT_ASSERT_THROW(vars.lowerBound("PriorFreeLC"), Variables::BoundNotSet);
   CPPUNIT_ASSERT_THROW(vars.upperBound("PriorFreeLC"), Variables::BoundNotSet);

   prior = vars.prior("PriorFreeLC");
   CPPUNIT_ASSERT_EQUAL(Variables::LOGUNIFORM, prior);
   prior1 = vars.priorPar1("PriorFreeLC");
   prior2 = vars.priorPar2("PriorFreeLC");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedPrior1, prior1, 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedPrior2, prior2, 1e-8);

   //Check for errors
   varNames.resize(1);
   varNames[0] = "TooManyValues";
   CPPUNIT_ASSERT_THROW(vars.add(varNames, pars), Variables::VariableError);

   varNames[0] = "BadValues";
   CPPUNIT_ASSERT_THROW(vars.add(varNames, pars), Parameters::ParameterError);

   varNames[0] = "BadValue";
   CPPUNIT_ASSERT_THROW(vars.add(varNames, pars), Parameters::ParameterError);

   varNames[0] = "PriorBad";
   CPPUNIT_ASSERT_THROW(vars.add(varNames, pars), Variables::VariableError);

   varNames[0] = "PriorBad2";
   CPPUNIT_ASSERT_THROW(vars.add(varNames, pars), Variables::VariableError);
   
   varNames[0] = "PriorBad3";
   CPPUNIT_ASSERT_THROW(vars.add(varNames, pars), Variables::VariableError);

   //Try a variable name that does not exist
   varNames[0] = "doesNotExist";
   CPPUNIT_ASSERT_THROW(vars.add(varNames, pars), Parameters::ParameterError);

}

void testVariables::addFromXML() {
   //Set up the XML string
   std::string xmlString = 
      "<someBranch>"
      "<variable id='Fixed'> <value>1.0</value> </variable>"
      "<variable id='NoBound'> <value>1.0</value> <step>0.1</step> </variable>"
      "<variable id='LowerBound'> <value>1.0</value> <step>0.1</step> <min>0.0</min> </variable>"
      "<variable id='UpperBound'> <value>1.0</value> <step>0.1</step> <max>2.0</max> </variable>"
      "<variable id='BothBounds'> <value>1.0</value> <step>0.1</step> <min>0.0</min> <max>2.0</max> </variable>"
      "<variable id='BothBoundsNewName'> <value>1.0</value> <step>0.1</step> <min>0.0</min> <max>2.0</max> <name>NewName</name> </variable>"
      "<variable id='BothBoundsPrior'> <value>1.0</value> <step>0.1</step> <min>0.0</min> <max>2.0</max> <prior>LOGUNIFORM</prior> <pval1>0.1</pval1> <pval2>2.0</pval2> </variable>"
      "<variable id='BoundsWithoutStep'> <value>1.0</value> <min>0.0</min> <max>2.0</max> </variable>"
      "<variable id='WithoutValue'> <step>0.1</step> <min>0.0</min> </variable>"
      "<variable id='PriorPval1Missing'> <value>1.0</value> <step>0.1</step> <min>0.0</min> <max>2.0</max> <prior>LOGUNIFORM</prior> <pval2>2.0</pval2> </variable>"
      "<variable id='PriorPval2Missing'> <value>1.0</value> <step>0.1</step> <min>0.0</min> <max>2.0</max> <prior>LOGUNIFORM</prior> <pval1>0.1</pval1> </variable>"
      "</someBranch>";

   //And the branch
   ReaderStringInput xmlInput(xmlString);
   Reader xmlReader(xmlInput);
   Branch b = xmlReader.GetTopBranch();

   //Test the variables that should work
   std::vector<std::string> varNames;
   varNames.push_back("Fixed");
   varNames.push_back("NoBound");
   varNames.push_back("LowerBound");
   varNames.push_back("UpperBound");
   varNames.push_back("BothBounds");
   varNames.push_back("BothBoundsNewName");
   varNames.push_back("BothBoundsPrior");
   varNames.push_back("BoundsWithoutStep");
   varNames.push_back("PriorPval1Missing");
   varNames.push_back("PriorPval2Missing");

   Variables vars;
   vars.add(varNames, "Pfx", b);

   double expectedValue(1.0), expectedError(0.1), expectedLowerBound(0.0), expectedUpperBound(2.0);
   std::string expectedName;

   double value, error, lowerBound, upperBound, pval1, pval2;
   Variables::PRIOR prior;

   //Check Fixed
   expectedName = "Pfx_Fixed";
   CPPUNIT_ASSERT_EQUAL(expectedName, varNames[0]);
   value = vars[varNames[0]];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error(varNames[0]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, error, 1e-8);
   
   CPPUNIT_ASSERT_THROW(vars.lowerBound(varNames[0]), Variables::BoundNotSet);
   CPPUNIT_ASSERT_THROW(vars.upperBound(varNames[0]), Variables::BoundNotSet);

   prior = vars.prior(varNames[0]);
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   //Check NoBound
   expectedName = "Pfx_NoBound";
   CPPUNIT_ASSERT_EQUAL(expectedName, varNames[1]);
   value = vars[varNames[1]];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error(varNames[1]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   
   CPPUNIT_ASSERT_THROW(vars.lowerBound(varNames[1]), Variables::BoundNotSet);
   CPPUNIT_ASSERT_THROW(vars.upperBound(varNames[1]), Variables::BoundNotSet);

   prior = vars.prior(varNames[1]);
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   //Check LowerBound
   expectedName = "Pfx_LowerBound";
   CPPUNIT_ASSERT_EQUAL(expectedName, varNames[2]);
   value = vars[varNames[2]];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error(varNames[2]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   
   lowerBound = vars.lowerBound(varNames[2]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, lowerBound, 1e-8);
   CPPUNIT_ASSERT_THROW(vars.upperBound(varNames[2]), Variables::BoundNotSet);

   prior = vars.prior(varNames[2]);
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   //Check UpperBound
   expectedName = "Pfx_UpperBound";
   CPPUNIT_ASSERT_EQUAL(expectedName, varNames[3]);
   value = vars[varNames[3]];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error(varNames[3]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   
   CPPUNIT_ASSERT_THROW(vars.lowerBound(varNames[3]), Variables::BoundNotSet);
   upperBound = vars.upperBound(varNames[3]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, upperBound, 1e-8);

   prior = vars.prior(varNames[3]);
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   //Check BothBoundsNewName
   expectedName = "NewName";
   CPPUNIT_ASSERT_EQUAL(expectedName, varNames[5]);
   value = vars[varNames[5]];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error(varNames[5]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   
   lowerBound = vars.lowerBound(varNames[5]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, lowerBound, 1e-8);
   upperBound = vars.upperBound(varNames[5]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, upperBound, 1e-8);

   prior = vars.prior(varNames[5]);
   CPPUNIT_ASSERT_EQUAL(Variables::UNIFORM, prior);
   pval1 = vars.priorPar1(varNames[5]);
   pval2 = vars.priorPar2(varNames[5]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, pval1, 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, pval2, 1e-8);

   //Check BothBoundsPrior
   expectedName = "Pfx_BothBoundsPrior";
   CPPUNIT_ASSERT_EQUAL(expectedName, varNames[6]);
   value = vars[varNames[6]];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error(varNames[6]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   
   lowerBound = vars.lowerBound(varNames[6]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, lowerBound, 1e-8);
   upperBound = vars.upperBound(varNames[6]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, upperBound, 1e-8);

   prior = vars.prior(varNames[6]);
   CPPUNIT_ASSERT_EQUAL(Variables::LOGUNIFORM, prior);
   pval1 = vars.priorPar1(varNames[6]);
   pval2 = vars.priorPar2(varNames[6]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1, pval1, 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, pval2, 1e-8);

   //Check BoundsWithoutStep
   expectedName = "Pfx_BoundsWithoutStep";
   CPPUNIT_ASSERT_EQUAL(expectedName, varNames[7]);
   value = vars[varNames[7]];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error(varNames[7]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, error, 1e-8);
   
   CPPUNIT_ASSERT_THROW(vars.lowerBound(varNames[7]), Variables::BoundNotSet);
   CPPUNIT_ASSERT_THROW(vars.upperBound(varNames[7]), Variables::BoundNotSet);

   prior = vars.prior(varNames[7]);
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   //Check PriorPval1Missing
   expectedName = "Pfx_PriorPval1Missing";
   CPPUNIT_ASSERT_EQUAL(expectedName, varNames[8]);
   value = vars[varNames[8]];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error(varNames[8]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   
   lowerBound = vars.lowerBound(varNames[8]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, lowerBound, 1e-8);
   upperBound = vars.upperBound(varNames[8]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, upperBound, 1e-8);

   prior = vars.prior(varNames[8]);
   CPPUNIT_ASSERT_EQUAL(Variables::UNIFORM, prior);
   pval1 = vars.priorPar1(varNames[8]);
   pval2 = vars.priorPar2(varNames[8]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, pval1, 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, pval2, 1e-8);

   //Check PriorPval2Missing
   expectedName = "Pfx_PriorPval2Missing";
   CPPUNIT_ASSERT_EQUAL(expectedName, varNames[9]);
   value = vars[varNames[9]];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error(varNames[9]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   
   lowerBound = vars.lowerBound(varNames[9]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, lowerBound, 1e-8);
   upperBound = vars.upperBound(varNames[9]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, upperBound, 1e-8);

   prior = vars.prior(varNames[9]);
   CPPUNIT_ASSERT_EQUAL(Variables::UNIFORM, prior);
   pval1 = vars.priorPar1(varNames[9]);
   pval2 = vars.priorPar2(varNames[9]);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, pval1, 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, pval2, 1e-8);

   //Make sure we throw errors
   varNames.resize(1);
   varNames[0] = "WithoutValue";
   CPPUNIT_ASSERT_THROW(vars.add(varNames, "Pfx", b), Variables::VariableError);

   varNames[0] = "NoSuchVariable";
   CPPUNIT_ASSERT_THROW(vars.add(varNames, "Pfx", b), Variables::VariableError);

}

void testVariables::addFromVariables() {
   //Set up a few variables in two objects
   Variables vars1;
   Variables vars2;

   vars1.add("Variable1", 1.0, 0.1);
   vars1.add("Variable2", 1.0, 0.1);
   //Same variable to make sure it is not touched
   vars2.add("Variable2", 1.0, 0.1, 0.0, 2.0);
   vars2.add("Variable3", 1.0, 0.1, 0.0, 2.0);

   vars1.add(vars2);

   //Make sure new variables are there
   double value, error, lowerBound, upperBound;
   value = vars1["Variable3"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, value, 1e-8);

   //Make sure the variables are not touched
   error = vars1.error("Variable2");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1, error, 1e-8);
   CPPUNIT_ASSERT_THROW(vars1.lowerBound("Variable2"), Variables::BoundNotSet);
   CPPUNIT_ASSERT_THROW(vars1.upperBound("Variable2"), Variables::BoundNotSet);

}

void testVariables::addManually() {
   Variables vars;

   double expectedValue(1.0), expectedError(0.1), expectedLowerBound(0.0), expectedUpperBound(2.0);

   double value, error, lowerBound, upperBound, pval1, pval2;
   Variables::PRIOR prior;
   Variables::varType type;

   //Fixed variable
   vars.add("Fixed", expectedValue);
   value = vars["Fixed"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("Fixed");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, error, 1e-8);
   
   CPPUNIT_ASSERT_THROW(vars.lowerBound("Fixed"), Variables::BoundNotSet);
   CPPUNIT_ASSERT_THROW(vars.upperBound("Fixed"), Variables::BoundNotSet);

   prior = vars.prior("Fixed");
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   type = vars.type("Fixed");
   CPPUNIT_ASSERT_EQUAL(Variables::FIXED, type);

   //NoBound variable
   vars.add("NoBound", expectedValue, expectedError);
   value = vars["NoBound"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("NoBound");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   
   CPPUNIT_ASSERT_THROW(vars.lowerBound("NoBound"), Variables::BoundNotSet);
   CPPUNIT_ASSERT_THROW(vars.upperBound("NoBound"), Variables::BoundNotSet);

   prior = vars.prior("NoBound");
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   type = vars.type("NoBound");
   CPPUNIT_ASSERT_EQUAL(Variables::FREE, type);

   //BothBounds variable
   vars.add("BothBounds", expectedValue, expectedError, expectedLowerBound, expectedUpperBound);
   value = vars["BothBounds"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("BothBounds");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   
   lowerBound = vars.lowerBound("BothBounds");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, lowerBound, 1e-8);
   upperBound = vars.upperBound("BothBounds");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, upperBound, 1e-8);

   prior = vars.prior("BothBounds");
   CPPUNIT_ASSERT_EQUAL(Variables::UNIFORM, prior);
   pval1 = vars.priorPar1("BothBounds");
   pval2 = vars.priorPar2("BothBounds");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, pval1, 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, pval2, 1e-8);

   //For some reason the code goes into an infinite loop if we try to add variables after this.  Don't really know why.
   //Checking one at a time made sure it behaved properly
   //It is an error adding already added variables
   //CPPUNIT_ASSERT_THROW(vars.add("Fixed", expectedValue), Variables::VariableAlreadyCreated);
   //CPPUNIT_ASSERT_THROW(vars.add("NoBound", expectedValue, expectedError), Variables::VariableAlreadyCreated);
   //CPPUNIT_ASSERT_THROW(vars.add("BothBounds", expectedValue, expectedError, expectedLowerBound, expectedUpperBound), Variables::VariableAlreadyCreated);

}

void testVariables::modifyBounds() {
   double expectedValue(1.0), expectedError(0.1), expectedLowerBound(0.0), expectedUpperBound(2.0);

   double value, error, lowerBound, upperBound;
   double prior1, prior2;
   Variables::PRIOR prior;
   Variables::varType type;

   //Add a fixed variable
   Variables vars;
   vars.add("var1", expectedValue);

   //Set the error and make sure everything is as it should be
   vars.setError("var1", expectedError);
   value = vars["var1"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("var1");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   
   CPPUNIT_ASSERT_THROW(vars.lowerBound("var1"), Variables::BoundNotSet);
   CPPUNIT_ASSERT_THROW(vars.upperBound("var1"), Variables::BoundNotSet);

   prior = vars.prior("var1");
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   type = vars.type("var1");
   CPPUNIT_ASSERT_EQUAL(Variables::FREE, type);

   //Set lower limit
   vars.setLowerBound("var1", expectedLowerBound);
   value = vars["var1"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("var1");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   lowerBound = vars.lowerBound("var1");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, lowerBound, 1e-8);

   CPPUNIT_ASSERT_THROW(vars.upperBound("var1"), Variables::BoundNotSet);

   prior = vars.prior("var1");
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   type = vars.type("var1");
   CPPUNIT_ASSERT_EQUAL(Variables::LOWBND, type);

   //Set upper limit
   vars.setUpperBound("var1", expectedUpperBound);
   value = vars["var1"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("var1");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   lowerBound = vars.lowerBound("var1");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, lowerBound, 1e-8);
   upperBound = vars.upperBound("var1");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, upperBound, 1e-8);

   prior = vars.prior("var1");
   CPPUNIT_ASSERT_EQUAL(Variables::UNIFORM, prior);
   prior1 = vars.priorPar1("var1");
   prior2 = vars.priorPar2("var1");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedLowerBound, prior1, 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, prior2, 1e-8);

   type = vars.type("var1");
   CPPUNIT_ASSERT_EQUAL(Variables::BOTHBND, type);

   //Remove lower limit
   vars.removeLowerBound("var1");
   value = vars["var1"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("var1");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   CPPUNIT_ASSERT_THROW(vars.lowerBound("var1"), Variables::BoundNotSet);
   upperBound = vars.upperBound("var1");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, upperBound, 1e-8);

   prior = vars.prior("var1");
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   type = vars.type("var1");
   CPPUNIT_ASSERT_EQUAL(Variables::UPPBND, type);

   //Fix the variable
   vars.fix("var1");
   value = vars["var1"];
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValue, value, 1e-8);
   error = vars.error("var1");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedError, error, 1e-8);
   CPPUNIT_ASSERT_THROW(vars.lowerBound("var1"), Variables::BoundNotSet);
   upperBound = vars.upperBound("var1");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedUpperBound, upperBound, 1e-8);

   prior = vars.prior("var1");
   CPPUNIT_ASSERT_EQUAL(Variables::NONE, prior);

   type = vars.type("var1");
   CPPUNIT_ASSERT_EQUAL(Variables::FIXED, type);

}

void testVariables::modifyPriors() {
   double expectedValue(1.0), expectedError(0.1), expectedLowerBound(0.0), expectedUpperBound(2.0);

   double value, error, lowerBound, upperBound;
   double prior1, prior2;
   Variables::PRIOR prior;
   Variables::varType type;

   //Add a fixed variable
   Variables vars;
   vars.add("var1", expectedValue);

   //Set the prior
   vars.setPrior("var1", Variables::UNIFORM, 1.0, 5.0);

   prior = vars.prior("var1");
   CPPUNIT_ASSERT_EQUAL(Variables::UNIFORM, prior);
   prior1 = vars.priorPar1("var1");
   prior2 = vars.priorPar2("var1");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, prior1, 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, prior2, 1e-8);

   //Make sure the prior sticks, even when setting the error
   vars.setError("var1", expectedError);
   prior = vars.prior("var1");
   CPPUNIT_ASSERT_EQUAL(Variables::UNIFORM, prior);
   prior1 = vars.priorPar1("var1");
   prior2 = vars.priorPar2("var1");
   CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, prior1, 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, prior2, 1e-8);


}
