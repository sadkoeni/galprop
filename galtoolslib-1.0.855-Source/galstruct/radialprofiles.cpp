#include "radialprofiles.h"
#include <ErrorLogger.h>
#include <ReaderErrorReporter.h>
#include <Registry.h>

#include <cmath>
#include <cassert>
#include <stdexcept>
#include <iomanip>
#include <sstream>
#include <gsl/gsl_interp.h>


namespace GalacticStructure {

double sech ( double x )
{
   if (fabs(x) > 300.)
      return 0;

   const double expx = exp(x);
   return 2.0/ (expx + 1./expx);
}

double sech2 ( double x )
{
   const double s = sech(x);
   return s*s;
}

RadialProfile::RadialProfile( const utl::Branch &b ) {

   assert(b.GetBranchNameString() == "RadialProfile");

   utl::Branch nb = b.GetChild("name");
   if ( ! nb) {
      std::ostringstream os;
      os<<"A name is required for radial profiles. Please specify the name element.";
      FATAL(os.str());
      throw(std::runtime_error("No name element found"));
   }
   nb.GetData(profileName);

}

std::unique_ptr<RadialProfile> RadialProfile::createProfile(const utl::Branch &b ) {
   assert(b.GetBranchNameString() == "RadialProfile");

   const std::string name = b.GetAttributes()["type"];

   return utl::Registry1<RadialProfile, const utl::Branch&>::create(name,b);
}

#ifdef HAVE_OPENCL
   std::string RadialProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
   {
      ERROR("OpenCL implementation not ready for this class");
      throw(std::runtime_error("Not Implemented"));
   }

   size_t RadialProfile::getOpenCLNPars() const
   {
      ERROR("OpenCL implementation not ready for this class");
      throw(std::runtime_error("Not Implemented"));
   }

   std::vector<cl_float> RadialProfile::getOpenCLPars() const
   {
      ERROR("OpenCL implementation not ready for this class");
      throw(std::runtime_error("Not Implemented"));
   }
#endif

void RadialProfile::addToDOMhelper( xercesc::DOMNode *node, 
      const std::string &type,
      const std::map<std::string,std::string> &attributes,
      const std::map<std::string,std::string> &elements, 
      const std::vector<std::string> &varNames,
      const std::vector<std::string> &varIds
      ) const
{

   xercesc::DOMDocument *doc = node->getOwnerDocument();

   xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("RadialProfile").unicodeForm());
   node->appendChild(profileEl);
   profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr(type).unicodeForm());

   for ( auto it = attributes.begin(); it != attributes.end(); ++it )
      profileEl->setAttribute(utl::XStr(it->first).unicodeForm(), utl::XStr(it->second).unicodeForm());

   //Add the name
   {
      xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
      profileEl->appendChild(element);

      xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
      element->appendChild(text);
   }

   for ( auto it = elements.begin(); it != elements.end(); ++it ) {
      xercesc::DOMElement* element = doc->createElement(utl::XStr(it->first).unicodeForm());
      profileEl->appendChild(element);

      xercesc::DOMText* text = doc->createTextNode(utl::XStr(it->second).unicodeForm());
      element->appendChild(text);
   }

   fvars.addToDOM( profileEl, varNames, varIds, profileName );

}




static utl::Registry1<RadialProfile, const utl::Branch&>::Registrar<ConstantRadialProfile> registrarConstant("Constant");

ConstantRadialProfile::ConstantRadialProfile( const utl::Branch &b ) :
   RadialProfile(b)
{

   assert(b.GetAttributes()["type"] == "Constant");

}

void ConstantRadialProfile::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   std::map<std::string, std::string> emptyMap;
   std::vector<std::string> emptyVector;
   addToDOMhelper( node, "Constant", attributes, emptyMap, emptyVector, emptyVector );
}

#ifdef HAVE_OPENCL
std::string ConstantRadialProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   name = "ConstantRadialProfile_"+name;
   std::ostringstream codeStr;
   codeStr<<"inline float "<<name<<"(float R, __constant float *pars) { return 1; }\n";
   return codeStr.str();
}

size_t ConstantRadialProfile::getOpenCLNPars() const { return 0; }

std::vector<cl_float> ConstantRadialProfile::getOpenCLPars() const { return std::vector<cl_float>(); }

#endif



static utl::Registry1<RadialProfile, const utl::Branch&>::Registrar<PulsarRadialProfile> registrarPulsar("Pulsar");

PulsarRadialProfile::PulsarRadialProfile( const utl::Branch &b ) :
   RadialProfile(b)
{

   assert(b.GetAttributes()["type"] == "Pulsar");

   fvarNames.resize(4);
   fvarNames[0] = "alpha";
   fvarNames[1] = "beta";
   fvarNames[2] = "Roff";
   fvarNames[3] = "R0";

   fvars.add(fvarNames, profileName, b);

   updateVariableValues( fvars );
}

PulsarRadialProfile::PulsarRadialProfile( const utl::Parameters &pars, const std::string &prefix ) 
{
   profileName = prefix;

   //Create the parameter names
   fvarNames.resize(4);
   fvarNames[0] = prefix + "_alpha";
   fvarNames[1] = prefix + "_beta";
   fvarNames[2] = prefix + "_Roff";
   fvarNames[3] = prefix + "_R0";

   //Add to the variables object
   fvars.add(fvarNames, pars);

   //Update cached values
   updateVariableValues( fvars );
}

PulsarRadialProfile::PulsarRadialProfile ( double alpha, double beta, double Roff, double R0, const std::string &prefix )
{
   profileName = prefix;

   //Create the parameter names
   fvarNames.resize(4);
   fvarNames[0] = prefix + "_alpha";
   fvarNames[1] = prefix + "_beta";
   fvarNames[2] = prefix + "_Roff";
   fvarNames[3] = prefix + "_R0";

   //Add to the variables object
   fvars.add(fvarNames[0], alpha);
   fvars.add(fvarNames[1], beta);
   fvars.add(fvarNames[2], Roff);
   fvars.add(fvarNames[3], R0);

   //Update cached values
   updateVariableValues( fvars );
}

PulsarRadialProfile::PulsarRadialProfile( const utl::Parameters &pars, const std::vector<std::string> &varNames ) :
   fvarNames(varNames)
{

   //Add to the variables object
   fvars.add(fvarNames, pars);

   //Update cached values
   updateVariableValues( fvars );
}

double PulsarRadialProfile::operator() (double R) const 
{
   //Apply the offset
   const double x = R + Roff;

   //The small number is there to allow negative Roff

   return norm * pow(x*oneOver_x0, alpha) * exp( -beta * (x-x0)*oneOver_x0 );
}

void PulsarRadialProfile::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   std::vector<std::string> ids(4);
   ids[0] = "alpha";
   ids[1] = "beta";
   ids[2] = "Roff";
   ids[3] = "R0";
   std::map<std::string, std::string> emptyMap;
   addToDOMhelper( node, "Pulsar", attributes, emptyMap, fvarNames, ids );
}

#ifdef HAVE_OPENCL
std::string PulsarRadialProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   name = "PulsarRadialProfile_"+name;
   std::ostringstream codeStr;
   codeStr<<"inline float "<<name<<"(float R, __constant float *pars)\n"
      <<"{\n"
      <<"return pars["<<parIndex+5<<"]*pow((R + pars["<<parIndex+2<<"])*pars["<<parIndex+4<<"], pars["<<parIndex<<"])\n"
      <<" * exp( -pars["<<parIndex+1<<"] * (R + pars["<<parIndex+2<<"]-pars["<<parIndex+3<<"])*pars["<<parIndex+4<<"]);\n" 
      <<"}\n";
   return codeStr.str();
}

size_t PulsarRadialProfile::getOpenCLNPars() const { return 6; }

std::vector<cl_float> PulsarRadialProfile::getOpenCLPars() const { 
   std::vector<cl_float> pars(6); 

   pars[0] = alpha;
   pars[1] = beta;
   pars[2] = Roff;
   pars[3] = x0;
   pars[4] = oneOver_x0;
   pars[5] = norm;

   return pars;
}
#endif

void PulsarRadialProfile::updateVariableValues( const utl::Variables &vars ) 
{

   //Update the values of the variables
   fvars.setFromVariables( vars );

   //Set the cache variables
   alpha = vars[fvarNames[0]];
   beta = vars[fvarNames[1]];
   Roff = vars[fvarNames[2]];
   x0 = vars[fvarNames[3]] + Roff;
   oneOver_x0 = 1./x0;
   norm = pow(beta/alpha, alpha) * exp(alpha-beta);

}

void PulsarRadialProfile::setAlpha(double a) {
   alpha = a;
   norm = pow(beta/alpha, alpha) * exp(alpha-beta);
}

void PulsarRadialProfile::setBeta(double b) {
   beta = b;
   norm = pow(beta/alpha, alpha) * exp(alpha-beta);
}

void PulsarRadialProfile::setRoff(double r) {
   Roff = r;
   x0 = fvars[fvarNames[3]] + Roff;
   oneOver_x0 = 1./x0;
}

void PulsarRadialProfile::setR0(double r) {
   fvars[fvarNames[3]] = r;
   x0 = r + Roff;
   oneOver_x0 = 1./x0;
}



static utl::Registry1<RadialProfile, const utl::Branch&>::Registrar<ExponentialRadialProfile> registrarExponential("Exponential");

ExponentialRadialProfile::ExponentialRadialProfile( const utl::Branch &b ) :
   RadialProfile(b)
{

   assert(b.GetAttributes()["type"] == "Exponential");

   std::vector<std::string> varNames(2);
   varNames[0] = "R0";
   varNames[1] = "RS";

   fvars.add(varNames, profileName, b);

   fR0name = varNames[0];
   fRSname = varNames[1];
   oneOver_R0 = 1./fvars[fR0name];
   rs = fvars[fRSname];
}

ExponentialRadialProfile::ExponentialRadialProfile ( const utl::Parameters &pars, const std::string &R0name, const std::string& RSname ) :
  fR0name(R0name), fRSname(RSname) {

   //Set the variables and the cache
   std::vector<std::string> varNames(2);
   varNames[0] = R0name;
   varNames[1] = RSname;
   fvars.add(varNames, pars);

   oneOver_R0 = 1./fvars[fR0name];
   rs = fvars[fRSname];

}

ExponentialRadialProfile::ExponentialRadialProfile ( double R0, double RS, const std::string &R0name, const std::string& RSname ) :
  fR0name(R0name), fRSname(RSname) {
  fvars.add(R0name, R0);
  fvars.add(RSname, RS);
  oneOver_R0 = 1./R0;
  rs = RS;
}

double ExponentialRadialProfile::operator () ( double R ) const {

  return exp(-(fabs(R) - rs)*oneOver_R0);

}
  
void ExponentialRadialProfile::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   std::vector<std::string> ids(2), names(2);
   ids[0] = "R0";
   ids[1] = "RS";
   names[0] = fR0name;
   names[1] = fRSname;
   std::map<std::string, std::string> emptyMap;
   addToDOMhelper( node, "Exponential", attributes, emptyMap, names, ids );
}

#ifdef HAVE_OPENCL
std::string ExponentialRadialProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   name = "ExponentialRadialProfile_"+name;
   std::ostringstream codeStr;
   codeStr<<"inline float "<<name<<"(float R, __constant float *pars)\n"
      <<"{\n"
      <<"return exp(-(fabs(R)-pars["<<parIndex<<"]) *pars["<<parIndex+1<<"]);\n"
      <<"}\n";
   return codeStr.str();
}

size_t ExponentialRadialProfile::getOpenCLNPars() const { return 2; }

std::vector<cl_float> ExponentialRadialProfile::getOpenCLPars() const { 
   std::vector<cl_float> pars(2); 

   pars[0] = rs;
   pars[1] = oneOver_R0;

   return pars;
}
#endif

void ExponentialRadialProfile::updateVariableValues ( const utl::Variables &vars ) {

  fvars.setFromVariables(vars);
  oneOver_R0 = 1./fvars[fR0name];
  rs = fvars[fRSname];

}

void ExponentialRadialProfile::setR0(double r) {
   oneOver_R0 = 1./r;
}

void ExponentialRadialProfile::setRS(double r) {
  
  rs = r;

}  




static utl::Registry1<RadialProfile, const utl::Branch&>::Registrar<FreudenreichWarpRadialProfile> registrarFreudenreichWarp("FreudenreichWarp");

FreudenreichWarpRadialProfile::FreudenreichWarpRadialProfile( const utl::Branch &b ) :
   RadialProfile(b)
{

  assert(b.GetAttributes()["type"] == "FreudenreichWarp");

  std::vector<std::string> varNames(4);
  varNames[0] = "Rw";
  varNames[1] = "c1";
  varNames[2] = "c2";
  varNames[3] = "c3";

  fvars.add(varNames, profileName, b);

  fRwname = varNames[0];
  fc1name = varNames[1];
  fc2name = varNames[2];
  fc3name = varNames[3];
  
  Rw = fvars[fRwname];
  c1 = fvars[fc1name];
  c2 = fvars[fc2name];
  c3 = fvars[fc3name];

}

  FreudenreichWarpRadialProfile::FreudenreichWarpRadialProfile ( const utl::Parameters &pars, const std::string &Rwname, const std::string& c1name, const std::string& c2name, const std::string& c3name ) :
    fRwname(Rwname), fc1name(c1name), fc2name(c2name), fc3name(c3name) {

    //Set the variables and the cache
    std::vector<std::string> varNames(4);
    varNames[0] = Rwname;
    varNames[1] = c1name;
    varNames[2] = c2name;
    varNames[3] = c3name;
    fvars.add(varNames, pars);

    Rw = fvars[fRwname];
    c1 = fvars[fc1name];
    c2 = fvars[fc2name];
    c3 = fvars[fc3name];

}

  FreudenreichWarpRadialProfile::FreudenreichWarpRadialProfile ( double Rw_, double c1_, double c2_, double c3_, const std::string &Rwname, const std::string& c1name, const std::string& c2name, const std::string& c3name ) :
    fRwname(Rwname), fc1name(c1name), fc2name(c2name), fc3name(c3name) {
  fvars.add(Rwname, Rw_);
  fvars.add(c1name, c1_);
  fvars.add(c2name, c2_);
  fvars.add(c3name, c3_);
  Rw = Rw_;
  c1 = c1_;
  c2 = c2_;
  c3 = c3_;
}

double FreudenreichWarpRadialProfile::operator () ( double R ) const {

  const auto u = R - Rw;

  return (u > 0. ? c1*u + c2*u*u + c3*u*u*u : 0.);

}
  
void FreudenreichWarpRadialProfile::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   std::vector<std::string> ids(4), names(4);
   ids[0] = "Rw";
   ids[1] = "c1";
   ids[2] = "c2";
   ids[3] = "c3";
   names[0] = fRwname;
   names[1] = fc1name;
   names[2] = fc2name;
   names[3] = fc3name;
   std::map<std::string, std::string> emptyMap;
   addToDOMhelper( node, "FreudenreichWarp", attributes, emptyMap, names, ids );
}

#ifdef HAVE_OPENCL
std::string FreudenreichWarpRadialProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   name = "FreudenreichWarpRadialProfile_"+name;
   std::ostringstream codeStr;
   codeStr<<"inline float "<<name<<"(float R, __constant float *pars)\n"
	  <<"{\n"
          <<"if ((R - pars["<<parIndex<<"]) < 0) return 0;\n"
	  <<"return (R - pars["<<parIndex<<"])*pars["<<parIndex+1<<"]+(R - pars["<<parIndex<<"])*(R - pars["<<parIndex<<"])*pars["<<parIndex+2<<"]+(R - pars["<<parIndex<<"])*(R - pars["<<parIndex<<"])*(R - pars["<<parIndex<<"])*pars["<<parIndex+3<<"];\n"
	  <<"}\n";
   return codeStr.str();
}

size_t FreudenreichWarpRadialProfile::getOpenCLNPars() const { return 4; }

std::vector<cl_float> FreudenreichWarpRadialProfile::getOpenCLPars() const { 
   std::vector<cl_float> pars(4); 

   pars[0] = Rw;
   pars[1] = c1;
   pars[2] = c2;
   pars[3] = c3;

   return pars;
}
#endif

void FreudenreichWarpRadialProfile::updateVariableValues ( const utl::Variables &vars ) {

  fvars.setFromVariables(vars);
  Rw = fvars[fRwname];
  c1 = fvars[fc1name];
  c2 = fvars[fc2name];
  c3 = fvars[fc3name];

}

void FreudenreichWarpRadialProfile::setRw(double r) {
   Rw = r;
}

void FreudenreichWarpRadialProfile::setc1(double c) {
  
  c1 = c;

}  

void FreudenreichWarpRadialProfile::setc2(double c) {
  
  c2 = c;

}

void FreudenreichWarpRadialProfile::setc3(double c) {
  
  c3 = c;

}




static utl::Registry1<RadialProfile, const utl::Branch&>::Registrar<CutOffPolynomialRadialProfile> registrarCutOffPolynomial("CutOffPolynomial");

CutOffPolynomialRadialProfile::CutOffPolynomialRadialProfile( const utl::Branch &b ) :
   RadialProfile(b)
{

   assert(b.GetAttributes()["type"] == "CutOffPolynomial");

   utl::Branch tb = b.GetChild("polynomialDegree");

   if ( !tb ) {
      std::ostringstream os;
      os<<"polynomialDegree element of radial profile \""<<profileName<<"\" not found.";
      FATAL(os.str());
      throw(std::runtime_error("polynomialDegree element required for cutoffpolynomial profiles"));
   }
      
   size_t nPoints(0);
   tb.GetData(nPoints);

   ci.resize(nPoints+1);
   ciname.resize(nPoints+1);

   //Construct the variable names and add them
   std::vector<std::string> varNames(nPoints+2);
   for (size_t i(0); i < nPoints+1; ++i) {
      std::ostringstream os;
      os << "c" << std::setfill('0') << std::setw(2) << i;
      varNames[i] = os.str();
   }
   varNames.back() = "Rc";

   fvars.add(varNames, profileName, b);

   //Store the actual values of the variable names
   for (size_t i(0); i < ciname.size(); ++i)
      ciname[i] = varNames[i];

   fRcname = varNames.back();

   updateVariableValues(fvars);
}

CutOffPolynomialRadialProfile::CutOffPolynomialRadialProfile ( const utl::Parameters &pars, const std::string& Rcname, const std::vector<std::string> &cinames ) :
    fRcname(Rcname), 
    ciname(cinames) 
{

   if (ciname.size() == 0) {
      FATAL("Need at least a polynomial of degree 0");
      throw(std::runtime_error("Fix your code"));
   }

   std::vector<std::string> varNames(ciname.size()+1);
   for (size_t i(0); i < ciname.size(); ++i)
      varNames[i] = ciname[i];
   varNames.back() = Rcname;

   fvars.add(varNames, pars);

   ci.resize(ciname.size());

   updateVariableValues(fvars);
}

CutOffPolynomialRadialProfile::CutOffPolynomialRadialProfile ( double Rc, std::vector<double> cis, const std::string& Rcname, const std::vector<std::string>& cinames ) :
   Rcut(Rc),
   ci(cis),
   fRcname(Rcname),
   ciname(cinames)
{

   if (ci.size() != ciname.size()) {
      FATAL("Values and names should have the same size");
      throw(std::runtime_error("Fix your code"));
   }

   if (ciname.size() == 0) {
      FATAL("Need at least polynomial of degree 0");
      throw(std::runtime_error("Fix your code"));
   }

   fvars.add(fRcname, Rc);
   
   for (size_t i(0); i < ciname.size(); ++i)
      fvars.add(ciname[i], ci[i]);

   updateVariableValues(fvars);
}

double CutOffPolynomialRadialProfile::operator () ( double R ) const {

  const auto u = R - Rcut;

  if ( u > 0 ) {
     double sum(ci[0]);
     for (size_t i(1); i < ci.size(); ++i)
        sum += ci[i]*pow(u,i);
     return sum;
  } else {
     return ci[0];
  }

}
  
void CutOffPolynomialRadialProfile::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   std::vector<std::string> ids(ciname.size()+1), names(ciname.size()+1);

   for (size_t i(0); i < ciname.size(); ++i) {
      std::ostringstream os;
      os << "c" << std::setfill('0') << std::setw(2) << i;
      ids[i] = os.str();
   }
   ids.back() = "Rc";

   for (size_t i(0); i < ciname.size(); ++i) 
      names[i] = ciname[i];
   names.back() = fRcname;

   std::map<std::string, std::string> nodes;
   std::ostringstream os;
   os << ciname.size() - 1;
   nodes["polynomialDegree"] = os.str();

   addToDOMhelper( node, "CutOffPolynomial", attributes, nodes, names, ids );
}

#ifdef HAVE_OPENCL
std::string CutOffPolynomialRadialProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   name = "CutOffPolynomialRadialProfile_"+name;
   std::ostringstream codeStr;
   codeStr<<"inline float "<<name<<"(float R, __constant float *pars)\n"
	  <<"{\n"
          <<"if ((R - pars["<<parIndex<<"]) < 0) return pars["<<parIndex+1<<"];\n"
	  <<"return pars["<<parIndex+1<<"]";
   for (size_t i = 1; i < ci.size(); ++i)
      codeStr<<" + pars["<<parIndex+1+i<<"]*pown((R - pars["<<parIndex<<"]),"<<i<<")";
   codeStr<<";\n"
	  <<"}\n";
   return codeStr.str();
}

size_t CutOffPolynomialRadialProfile::getOpenCLNPars() const { return ciname.size() + 1; }

std::vector<cl_float> CutOffPolynomialRadialProfile::getOpenCLPars() const { 
   std::vector<cl_float> pars(ciname.size()+1); 

   pars[0] = Rcut;
   for (size_t i = 1; i < pars.size(); ++i)
      pars[i] = ci[i-1];

   return pars;
}
#endif

void CutOffPolynomialRadialProfile::updateVariableValues ( const utl::Variables &vars ) {

  fvars.setFromVariables(vars);
  Rcut = fvars[fRcname];
  for (size_t i = 0; i < ci.size(); ++i)
     ci[i] = fvars[ciname[i]];

}

void CutOffPolynomialRadialProfile::setRc(double r) {
   Rcut = r;

   fvars[fRcname] = Rcut;
}

void CutOffPolynomialRadialProfile::setci(std::vector<double> c) {
  
   if (c.size() != ci.size()) {
      FATAL("Size of vector should match polynomial degree");
      throw(std::runtime_error("Fix your code"));
   }

   for (size_t i(0); i < ci.size(); ++i) {
      ci[i] = c[i];
      fvars[ciname[i]] = ci[i];
   }
}  



static utl::Registry1<RadialProfile, const utl::Branch&>::Registrar<ExpHoleRadialProfile> registrarExpHole("ExpHole");

ExpHoleRadialProfile::ExpHoleRadialProfile( const utl::Branch &b ) :
   RadialProfile(b)
{

   assert(b.GetAttributes()["type"] == "ExpHole");

   std::vector<std::string> varNames(4);
   varNames[0] = "R0";
   varNames[1] = "Rs";
   varNames[2] = "Rh";
   varNames[3] = "hi";

   fvars.add(varNames, profileName, b);

   fR0name = varNames[0];
   fRsname = varNames[1];
   fRhname = varNames[2];
   fhiname = varNames[3];
   oneOver_R0 = 1./fvars[fR0name];
   Rs = fvars[fRsname];
   oneOver_Rh = 1./fvars[fRhname];
   hi = fvars[fhiname];

}

  ExpHoleRadialProfile::ExpHoleRadialProfile ( const utl::Parameters &pars, const std::string &R0name, const std::string &Rsname, const std::string &Rhname, const std::string &hiname) :
   fR0name(R0name),
   fRsname(Rsname),
   fRhname(Rhname),
   fhiname(hiname)
{

   //Set the variables and the cache
   std::vector<std::string> varNames(4);
   varNames[0] = R0name;
   varNames[1] = Rsname;
   varNames[2] = Rhname;
   varNames[3] = hiname;

   fvars.add(varNames, pars);

   oneOver_R0 = 1./fvars[fR0name];
   Rs = fvars[fRsname];
   oneOver_Rh = 1./fvars[fRhname];
   hi = fvars[fhiname];
}

  ExpHoleRadialProfile::ExpHoleRadialProfile ( double R0, double Rs_, double Rh, double hi_, const std::string &R0name, const std::string &Rsname, const std::string &Rhname, const std::string &hiname ) :
   fR0name(R0name),
   fRsname(Rsname),
   fRhname(Rhname),
   fhiname(hiname)
{
   fvars.add(R0name, R0);
   fvars.add(Rsname, Rs_);
   fvars.add(Rhname, Rh);
   fvars.add(hiname, hi_);
   oneOver_R0 = 1./R0;
   Rs = Rs_;
   oneOver_Rh = 1./Rh;
   hi = hi_;
}

double ExpHoleRadialProfile::operator () ( double R ) const 
{
  
   return exp(-(R-Rs)*oneOver_R0) * (1 - exp(-pow(R*oneOver_Rh, hi)));
}

void ExpHoleRadialProfile::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   std::vector<std::string> ids(4), names(4);
   ids[0] = "R0";
   ids[1] = "Rs";
   ids[2] = "Rh";
   ids[3] = "hi";
   names[0] = fR0name;
   names[1] = fRsname;
   names[2] = fRhname;
   names[3] = fhiname;
   std::map<std::string, std::string> emptyMap;
   addToDOMhelper( node, "ExpHole", attributes, emptyMap, names, ids );
}

#ifdef HAVE_OPENCL
std::string ExpHoleRadialProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   name = "ExpHoleRadialProfile_"+name;
   std::ostringstream codeStr;
   codeStr<<"inline float "<<name<<"(float R, __constant float *pars)\n"
      <<"{\n"
      <<"return exp(-(R-pars["<<parIndex<<"]) * pars["<<parIndex+1<<"]) * (1 - exp(-pow(R*pars["<<parIndex+2<<"], pars["<<parIndex+3<<"])));\n"
      <<"}\n";
   return codeStr.str();
}

size_t ExpHoleRadialProfile::getOpenCLNPars() const { return 4; }

std::vector<cl_float> ExpHoleRadialProfile::getOpenCLPars() const { 
   std::vector<cl_float> pars(4); 

   pars[0] = Rs;
   pars[1] = oneOver_R0;
   pars[2] = oneOver_Rh;
   pars[3] = hi;

   return pars;
}
#endif 

void ExpHoleRadialProfile::updateVariableValues ( const utl::Variables &vars ) 
{
   fvars.setFromVariables(vars);
   oneOver_R0 = 1./fvars[fR0name];
   Rs = fvars[fRsname];
   oneOver_Rh = 1./fvars[fRhname];
   hi = fvars[fhiname];
}

void ExpHoleRadialProfile::setR0(double r) {
   oneOver_R0 = 1./r;
}

void ExpHoleRadialProfile::setRs(double r) {
   Rs = r;
}

void ExpHoleRadialProfile::setRh(double r) {
   oneOver_Rh = 1./r;
}

void ExpHoleRadialProfile::sethi(double r) {
   hi = r;
}




static utl::Registry1<RadialProfile, const utl::Branch&>::Registrar<ExpGaussianHoleRadialProfile> registrarExpGaussianHole("ExpGaussianHole");

ExpGaussianHoleRadialProfile::ExpGaussianHoleRadialProfile( const utl::Branch &b ) :
   RadialProfile(b)
{

   assert(b.GetAttributes()["type"] == "ExpGaussianHole");

   std::vector<std::string> varNames(3);
   varNames[0] = "Sigma0";
   varNames[1] = "Mu0";
   varNames[2] = "Rh";

   fvars.add(varNames, profileName, b);

   fSigma0name = varNames[0];
   fMu0name = varNames[1];
   fRhname = varNames[2];

   updateVariableValues(fvars);

   /*oneOver_Rh = 1./fvars[fRhname];
   oneOver_sigma0 = 1./fvars[fSigma0name];

   rSmooth = fvars[fSigma0name]*fvars[fSigma0name]*oneOver_Rh + fvars[fMu0name];
   f0 = std::exp(-rSmooth*oneOver_Rh)*std::pow(std::exp(-0.5*std::pow((rSmooth - fvars[fMu0name])*oneOver_sigma0, 2.)), -1.);
   */
}

  ExpGaussianHoleRadialProfile::ExpGaussianHoleRadialProfile ( const utl::Parameters &pars, const std::string &Sigma0name, const std::string &Mu0name, const std::string &Rhname) :
   fMu0name(Mu0name),
   fSigma0name(Sigma0name),
   fRhname(Rhname) {

   //Set the variables and the cache
   std::vector<std::string> varNames(3);
   varNames[0] = Sigma0name;
   varNames[1] = Mu0name;
   varNames[2] = Rhname;

   fvars.add(varNames, pars);

   updateVariableValues(fvars);

   /*oneOver_Rh = 1./fvars[fRhname];
   oneOver_sigma0 = 1./fvars[fSigma0name];

   rSmooth = fvars[fSigma0name]*fvars[fSigma0name]*oneOver_Rh + fvars[fMu0name];
   f0 = std::exp(-rSmooth*oneOver_Rh)*std::pow(std::exp(-0.5*std::pow((rSmooth - fvars[fMu0name])*oneOver_sigma0, 2.)), -1.);
   */
}

  ExpGaussianHoleRadialProfile::ExpGaussianHoleRadialProfile ( double Sigma0, double Mu0, double Rh, const std::string &Sigma0name, const std::string &Mu0name, const std::string &Rhname) :
   fMu0name(Mu0name),
   fSigma0name(Sigma0name),
   fRhname(Rhname) {

   fvars.add(Sigma0name, Sigma0);
   fvars.add(Mu0name, Mu0);
   fvars.add(Rhname, Rh);

   updateVariableValues(fvars);

   /*oneOver_Rh = 1./Rh;
   oneOver_sigma0 = 1./Sigma0;

   rSmooth = fvars[fSigma0name]*fvars[fSigma0name]*oneOver_Rh + fvars[fMu0name];
   f0 = std::exp(-rSmooth*oneOver_Rh)*std::pow(std::exp(-0.5*std::pow((rSmooth - fvars[fMu0name])*oneOver_sigma0, 2.)), -1.);
   */
}

double ExpGaussianHoleRadialProfile::operator () ( double R ) const 
{

  //std::cout << R << " " << 1./oneOver_Rh << " " << sigma0 << " " << mu0 << " " << f0 << " " << std::exp(-0.5*std::pow((R - mu0)*oneOver_sigma0, 2.)) << " " ;// << std::endl;
  
  return (R < rSmooth ? f0*std::exp(-0.5*std::pow((R - mu0)*oneOver_sigma0, 2.)) : std::exp(-R*oneOver_Rh));

}

void ExpGaussianHoleRadialProfile::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   std::vector<std::string> ids(3), names(3);
   ids[0] = "Sigma0";
   ids[1] = "Mu0";
   ids[2] = "Rh";
   names[0] = fSigma0name;
   names[1] = fMu0name;
   names[2] = fRhname;
   std::map<std::string, std::string> emptyMap;
   addToDOMhelper( node, "ExpGaussianHole", attributes, emptyMap, names, ids );
}

#ifdef HAVE_OPENCL
std::string ExpGaussianHoleRadialProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   name = "ExpGaussianHoleRadialProfile_"+name;
   std::ostringstream codeStr;
   codeStr<<"inline float "<<name<<"(float R, __constant float *pars)\n"
      <<"{\n"
      <<"if ( R < pown(pars["<<parIndex+1<<"],2)*pars["<<parIndex+2<<"] + pars["<<parIndex<<"] ) \n"
      <<"return exp( float(0.5)*( pown( pars["<<parIndex+1<<"]*pars["<<parIndex+2<<"], 2) - \n"
      <<"                  pown( ( R - pars["<<parIndex<<"] ) / pars["<<parIndex+1<<"], 2) ) - \n"
      <<"            (pown(pars["<<parIndex+1<<"],2)*pars["<<parIndex+2<<"] + pars["<<parIndex<<"])*pars["<<parIndex+2<<"]);\n"
      <<"else\n;"
      <<"return exp(-R*pars["<<parIndex+2<<"]);\n"
      <<"}\n";
   return codeStr.str();
}

size_t ExpGaussianHoleRadialProfile::getOpenCLNPars() const { return 3; }

std::vector<cl_float> ExpGaussianHoleRadialProfile::getOpenCLPars() const { 
   std::vector<cl_float> pars(3); 

   pars[0] = mu0;
   pars[1] = sigma0;
   pars[2] = oneOver_Rh;

   return pars;
}
#endif 


void ExpGaussianHoleRadialProfile::updateVariableValues ( const utl::Variables &vars ) 
{
   fvars.setFromVariables(vars);
   mu0 = fvars[fMu0name];
   sigma0 = fvars[fSigma0name];
   oneOver_sigma0 = 1./sigma0;
   oneOver_Rh = 1./fvars[fRhname];

   rSmooth = sigma0*sigma0*oneOver_Rh + mu0;
   f0 = std::exp(-rSmooth*oneOver_Rh)*std::pow(std::exp(-0.5*std::pow((rSmooth - mu0)*oneOver_sigma0, 2.)), -1.);

}

void ExpGaussianHoleRadialProfile::setMu0(double mu0_) {
   mu0 = mu0_;
   rSmooth = sigma0*sigma0*oneOver_Rh + mu0;
   f0 = std::exp(-rSmooth*oneOver_Rh)*std::pow(std::exp(-0.5*std::pow((rSmooth - mu0)*oneOver_sigma0, 2.)), -1.);
}

void ExpGaussianHoleRadialProfile::setSigma0(double sigma0_) {
   sigma0 = sigma0_;
   oneOver_sigma0 = 1./sigma0;
   rSmooth = sigma0*sigma0*oneOver_Rh + mu0;
   f0 = std::exp(-rSmooth*oneOver_Rh)*std::pow(std::exp(-0.5*std::pow((rSmooth - mu0)*oneOver_sigma0, 2.)), -1.);
}

void ExpGaussianHoleRadialProfile::setRh(double r) {
   oneOver_Rh = 1./r;
   rSmooth = sigma0*sigma0*oneOver_Rh + mu0;
   f0 = std::exp(-rSmooth*oneOver_Rh)*std::pow(std::exp(-0.5*std::pow((rSmooth - mu0)*oneOver_sigma0, 2.)), -1.);
}




static utl::Registry1<RadialProfile, const utl::Branch&>::Registrar<Sech2RadialProfile> registrarSech2("Sech2");

Sech2RadialProfile::Sech2RadialProfile( const utl::Branch &b ) :
   RadialProfile(b)
{

   assert(b.GetAttributes()["type"] == "Sech2");

   std::vector<std::string> varNames(1);
   varNames[0] = "R0";

   fvars.add(varNames, profileName, b);

   fR0name = varNames[0];
   oneOver_R0 = 1./fvars[fR0name];
}

Sech2RadialProfile::Sech2RadialProfile ( const utl::Parameters &pars, const std::string &R0name ) :
   fR0name(R0name)
{

   //Set the variables and the cache
   std::vector<std::string> varNames(1);
   varNames[0] = R0name;
   fvars.add(varNames, pars);

   oneOver_R0 = 1./fvars[fR0name];
}

Sech2RadialProfile::Sech2RadialProfile ( double R0, const std::string &R0name ) :
   fR0name(R0name)
{
   fvars.add(R0name, R0);
   oneOver_R0 = 1./R0;
}

double Sech2RadialProfile::operator () ( double R ) const 
{
   return sech2(R*oneOver_R0);
}

void Sech2RadialProfile::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   std::vector<std::string> ids(1), names(1);
   ids[0] = "R0";
   names[0] = fR0name;
   std::map<std::string, std::string> emptyMap;
   addToDOMhelper( node, "Sech2", attributes, emptyMap, names, ids );
}

#ifdef HAVE_OPENCL
std::string Sech2RadialProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   name = "Sech2RadialProfile_"+name;
   std::ostringstream codeStr;
   codeStr<<"inline float "<<name<<"(float R, __constant float *pars)\n"
      <<"{\n"
      <<"return 4./(exp(2*R*pars["<<parIndex<<"]) + 2 + exp(-2*R*pars["<<parIndex<<"]));\n"
      <<"}\n";
   return codeStr.str();
}

size_t Sech2RadialProfile::getOpenCLNPars() const { return 1; }

std::vector<cl_float> Sech2RadialProfile::getOpenCLPars() const { 
   std::vector<cl_float> pars(1); 

   pars[0] = oneOver_R0;

   return pars;
}
#endif 

void Sech2RadialProfile::updateVariableValues ( const utl::Variables &vars ) 
{
   fvars.setFromVariables(vars);
   oneOver_R0 = 1./fvars[fR0name];
}

void Sech2RadialProfile::setR0(double r) {
   oneOver_R0 = 1./r;
}



static utl::Registry1<RadialProfile, const utl::Branch&>::Registrar<GaussianRadialProfile> registrarGaussian("Gaussian");

GaussianRadialProfile::GaussianRadialProfile( const utl::Branch &b ) :
   RadialProfile(b)
{

   assert(b.GetAttributes()["type"] == "Gaussian");

   std::vector<std::string> varNames(2);
   varNames[0] = "R0";
   varNames[1] = "Roff";

   fvars.add(varNames, profileName, b);

   fR0name = varNames[0];
   fRoffname = varNames[1];
   oneOver_R02 = 1./pow(fvars[fR0name],2);
   Roff = fvars[fRoffname];
}

GaussianRadialProfile::GaussianRadialProfile ( const utl::Parameters &pars, const std::string &R0name, const std::string &Roffname ) :
   fR0name(R0name),
   fRoffname(Roffname)
{

   //Set the variables and the cache
   std::vector<std::string> varNames(2);
   varNames[0] = R0name;
   varNames[1] = Roffname;
   fvars.add(varNames, pars);

   oneOver_R02 = 1./pow(fvars[fR0name],2);
   Roff = fvars[fRoffname];
}

GaussianRadialProfile::GaussianRadialProfile ( double R0, double Roff_, const std::string &R0name, const std::string &Roffname ) :
   fR0name(R0name),
   fRoffname(Roffname)
{

   //Set the variables and the cache
   fvars.add(R0name, R0);
   fvars.add(Roffname, Roff_);

   oneOver_R02 = 1./pow(R0,2);
   Roff = Roff_;
}

double GaussianRadialProfile::operator () ( double R ) const 
{
   const double x = R + Roff;
   return exp(-x*x*oneOver_R02);
}

void GaussianRadialProfile::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   std::vector<std::string> ids(2), names(2);
   ids[0] = "R0";
   ids[1] = "Roff";
   names[0] = fR0name;
   names[1] = fRoffname;
   std::map<std::string, std::string> emptyMap;
   addToDOMhelper( node, "Gaussian", attributes, emptyMap, names, ids );
}

#ifdef HAVE_OPENCL
std::string GaussianRadialProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   name = "GaussianRadialProfile_"+name;
   std::ostringstream codeStr;
   codeStr<<"inline float "<<name<<"(float R, __constant float *pars)\n"
      <<"{\n"
      <<"return exp(- (R+pars["<<parIndex<<"])*(R+pars["<<parIndex<<"])*pars["<<parIndex+1<<"]);\n"
      <<"}\n";
   return codeStr.str();
}

size_t GaussianRadialProfile::getOpenCLNPars() const { return 2; }

std::vector<cl_float> GaussianRadialProfile::getOpenCLPars() const { 
   std::vector<cl_float> pars(2); 

   pars[0] = Roff;
   pars[1] = oneOver_R02;

   return pars;
}
#endif 

void GaussianRadialProfile::updateVariableValues ( const utl::Variables &vars ) 
{
   fvars.setFromVariables(vars);
   oneOver_R02 = 1./pow(fvars[fR0name],2);
   Roff = fvars[fRoffname];
}

void GaussianRadialProfile::setR0(double r) {
   oneOver_R02 = 1./(r*r);
}

void GaussianRadialProfile::setRoff(double r) {
   Roff = r;
}



static utl::Registry1<RadialProfile, const utl::Branch&>::Registrar<MultipleGaussianRadialProfile> registrarMultipleGaussian("MultipleGaussian");

MultipleGaussianRadialProfile::MultipleGaussianRadialProfile( const utl::Branch &b ) :
   RadialProfile(b)
{

   assert(b.GetAttributes()["type"] == "MultipleGaussian");

   utl::Branch tb = b.GetChild("numberOfGaussians");

   if ( !tb ) {
      std::ostringstream os;
      os<<"numberOfGaussians element of radial profile \""<<profileName<<"\" not found.";
      FATAL(os.str());
      throw(std::runtime_error("numberOfGaussians element required for MultipleGaussian profiles"));
   }
      
   size_t nGaussians(0);
   tb.GetData(nGaussians);

   //The variable names
   fnnames.resize(nGaussians);
   fR0names.resize(nGaussians);
   fRoffnames.resize(nGaussians);
   for (size_t i(0); i < nGaussians; ++i) {
      std::ostringstream os;
      os << "n_"<< std::setfill('0') << std::setw(2) << i;
      fnnames[i] = os.str();
      os.str("");
      os << "R0_"<< std::setfill('0') << std::setw(2) << i;
      fR0names[i] = os.str();
      os.str("");
      os << "Roff_"<< std::setfill('0') << std::setw(2) << i;
      fRoffnames[i] = os.str();
   }

   //Resize the parameter storage
   ni.resize(nGaussians);
   oneOver_R02.resize(nGaussians);
   Roff.resize(nGaussians);

   fvars.add(fnnames, profileName, b);
   fvars.add(fR0names, profileName, b);
   fvars.add(fRoffnames, profileName, b);

   updateVariableValues(fvars);
}

MultipleGaussianRadialProfile::MultipleGaussianRadialProfile ( const utl::Parameters &pars,
            const std::vector<std::string> &nnames, 
            const std::vector<std::string> &R0names, 
            const std::vector<std::string> &Roffnames ):
   fnnames(nnames),
   fR0names(R0names),
   fRoffnames(Roffnames)
{
   //Make sure they are off the same size
   assert( nnames.size() == R0names.size() );
   assert( nnames.size() == Roffnames.size() );

   fvars.add(fnnames, pars);
   fvars.add(fR0names, pars);
   fvars.add(fRoffnames, pars);

   //Resize the parameter storage
   ni.resize(fnnames.size());
   oneOver_R02.resize(fnnames.size());
   Roff.resize(fnnames.size());

   updateVariableValues(fvars);
}

MultipleGaussianRadialProfile::MultipleGaussianRadialProfile (  const utl::Parameters &pars, const std::string &prefix, const size_t nGaussians ) 
{
   profileName = prefix;

   //The variable names
   fnnames.resize(nGaussians);
   fR0names.resize(nGaussians);
   fRoffnames.resize(nGaussians);
   for (size_t i(0); i < nGaussians; ++i) {
      std::ostringstream os;
      os << prefix << "_n_"<< std::setfill('0') << std::setw(2) << i;
      fnnames[i] = os.str();
      os.str("");
      os << prefix << "_R0_"<< std::setfill('0') << std::setw(2) << i;
      fR0names[i] = os.str();
      os.str("");
      os << prefix << "_Roff_"<< std::setfill('0') << std::setw(2) << i;
      fRoffnames[i] = os.str();
   }

   //Resize the parameter storage
   ni.resize(nGaussians);
   oneOver_R02.resize(nGaussians);
   Roff.resize(nGaussians);

   fvars.add(fnnames, pars);
   fvars.add(fR0names, pars);
   fvars.add(fRoffnames, pars);

   updateVariableValues(fvars);
}

double MultipleGaussianRadialProfile::operator () ( double R ) const 
{
   double out(0);
   for (size_t i(0); i < ni.size(); ++i) {
      const double x = R + Roff[i];
      out += ni[i] * exp(-x*x*oneOver_R02[i]);
   }
   return out;
}

void MultipleGaussianRadialProfile::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   std::map<std::string, std::string> elements;
   std::ostringstream oss;
   oss << fnnames.size();
   elements["numberOfGaussians"] = oss.str();
   oss.str("");

   std::vector<std::string> ids, names;
   ids.reserve(3*fnnames.size());
   names.reserve(3*fnnames.size());
   for (size_t i = 0; i < fnnames.size(); ++i) {
      oss << "n_"<< std::setfill('0') << std::setw(2) << i;
      ids.push_back(oss.str());
      oss.str("");

      oss << "R0_"<< std::setfill('0') << std::setw(2) << i;
      ids.push_back(oss.str());
      oss.str("");

      oss << "Roff_"<< std::setfill('0') << std::setw(2) << i;
      ids.push_back(oss.str());
      oss.str("");

      names.push_back(fnnames[i]);
      names.push_back(fR0names[i]);
      names.push_back(fR0names[i]);
   }

   addToDOMhelper( node, "MultipleGaussian", attributes, elements, names, ids );
}

#ifdef HAVE_OPENCL
std::string MultipleGaussianRadialProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   name = "MultipleGaussianRadialProfile_"+name;
   std::ostringstream codeStr;
   codeStr<<"inline float "<<name<<"(float R, __constant float *pars)\n"
      <<"{\n"
      <<"return ";
   for (size_t i(0); i < fnnames.size(); ++i)
      codeStr<<"pars["<<parIndex+3*i+0<<"]*exp(- (R+pars["<<parIndex+3*i+2<<"])*(R+pars["<<parIndex+3*i+2<<"])*pars["<<parIndex+3*i+1<<"]) + ";
   codeStr<<"0;\n"
      <<"}\n";
   return codeStr.str();
}

size_t MultipleGaussianRadialProfile::getOpenCLNPars() const { return 3*fnnames.size(); }

std::vector<cl_float> MultipleGaussianRadialProfile::getOpenCLPars() const { 
   std::vector<cl_float> pars(3*fnnames.size()); 

   for (size_t i(0); i < fnnames.size(); ++i) {
      pars[3*i+0] = ni[i];
      pars[3*i+1] = oneOver_R02[i];
      pars[3*i+2] = Roff[i];
   }

   return pars;
}
#endif 

void MultipleGaussianRadialProfile::updateVariableValues ( const utl::Variables &vars ) 
{
   fvars.setFromVariables(vars);
   for (size_t i(0); i < fnnames.size(); ++i) {
      ni[i] = fvars[fnnames[i]];
      oneOver_R02[i] = 1./pow(fvars[fR0names[i]],2);
      Roff[i] = fvars[fRoffnames[i]];
   }
}



static utl::Registry1<RadialProfile, const utl::Branch&>::Registrar<ConstantCoreRadialProfile> registrarConstantCore("ConstantCore");

ConstantCoreRadialProfile::ConstantCoreRadialProfile( const utl::Branch &b ) :
   RadialProfile(b)
{

   assert(b.GetAttributes()["type"] == "ConstantCore");

   std::vector<std::string> varNames(1);
   varNames[0] = "Rcore";

   fvars.add(varNames, profileName, b);

   fRcorename = varNames[0];
   Rcore = fvars[fRcorename];

   //Loop over the children to find the radial profiles
   for ( utl::Branch rp = b.GetFirstChild();  rp; rp = rp.GetNextSibling() ) {

      if ( rp.GetBranchNameString() == "RadialProfile" ) {

         fg = RadialProfile::createProfile( rp );
         break;
      }
   }

   //Make sure we found one
   if ( fg == 0 ) {
      FATAL("Missing radial profile for constant core profile");
      throw(std::runtime_error("Radial profile missing"));
   }

   fvars.add(fg->getVariables());
}

ConstantCoreRadialProfile::ConstantCoreRadialProfile (std::unique_ptr<RadialProfile> g, const utl::Parameters &pars, const std::string &Rcorename ) :
   fg(std::move(g)),
   fRcorename(Rcorename)
{

   //Set the variables and the cache
   std::vector<std::string> varNames(1);
   varNames[1] = Rcorename;
   fvars.add(varNames, pars);

   Rcore = fvars[fRcorename];

   fvars.add(fg->getVariables());
}

ConstantCoreRadialProfile::ConstantCoreRadialProfile ( std::unique_ptr<RadialProfile> g, double Rcore_, const std::string &Rcorename ) :
   fg(std::move(g)),
   fRcorename(Rcorename)
{

   //Set the variables and the cache
   fvars.add(Rcorename, Rcore_);

   Rcore = Rcore_;

   fvars.add(fg->getVariables());
}

double ConstantCoreRadialProfile::operator () ( double R ) const 
{
   return R < Rcore ? 1 : (*fg)(R-Rcore);
}

void ConstantCoreRadialProfile::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   std::vector<std::string> ids(1), names(1);
   ids[0] = "Rcore";
   names[0] = fRcorename;
   std::map<std::string, std::string> emptyMap;
   addToDOMhelper( node, "ConstantCore", attributes, emptyMap, names, ids );

   fg->addToDOM( node->getLastChild(), emptyMap );
}

#ifdef HAVE_OPENCL
std::string ConstantCoreRadialProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   name = "ConstantCoreRadialProfile_"+name;
   //Add the other function first
   std::string fname = name;
   std::ostringstream codeStr;
   codeStr<<fg->getOpenCLFunction(fname, parIndex);
   codeStr<<"inline float "<<name<<"(float R, __constant float *pars)\n"
      <<"{\n"
      <<"return R < pars["<<parIndex+fg->getOpenCLNPars()<<"] ? 1 : "<<fname<<"(R-pars["<<parIndex+fg->getOpenCLNPars()<<"],pars);\n"
      <<"}\n";
   return codeStr.str();
}

size_t ConstantCoreRadialProfile::getOpenCLNPars() const { return 1+fg->getOpenCLNPars(); }

std::vector<cl_float> ConstantCoreRadialProfile::getOpenCLPars() const { 
   std::vector<cl_float> pars = fg->getOpenCLPars(); 

   //Have our parameter last
   pars.push_back(Rcore);

   return pars;
}
#endif 

void ConstantCoreRadialProfile::updateVariableValues ( const utl::Variables &vars ) 
{
   fg->updateVariableValues(vars);
   fvars.setFromVariables(vars);
   Rcore = fvars[fRcorename];
}

void ConstantCoreRadialProfile::setRcore(double r) {
   Rcore = r;
}

void ConstantCoreRadialProfile::setg(std::unique_ptr<RadialProfile> g) {
   fg = std::move(g);
}



static utl::Registry1<RadialProfile, const utl::Branch&>::Registrar<SplineRadialProfile> registrarSpline("Spline");

SplineRadialProfile::SplineRadialProfile( const utl::Branch &b ) :
   RadialProfile(b),
   spline(0), 
   acc(0)
{

   //Need the substring because of SplineLog subclass
   assert(b.GetAttributes()["type"].substr(0,6) == "Spline");

   utl::Branch tb = b.GetChild("numberOfPoints");

   if ( !tb ) {
      std::ostringstream os;
      os<<"numberOfPoints element of radial profile \""<<profileName<<"\" not found.";
      FATAL(os.str());
      throw(std::runtime_error("numberOfPoints element required for Spline profiles"));
   }
      
   size_t nPoints(0);
   tb.GetData(nPoints);

   tb = b.GetChild("interpolationType");

   if ( !tb ) {
      std::ostringstream os;
      os<<"interpolationType element of radial profile \""<<profileName<<"\" not found.";
      FATAL(os.str());
      throw(std::runtime_error("interpolationType element required for Spline profiles"));
   }

   std::string upper;
   tb.GetData(upper);

   // explicit cast needed to resolve ambiguity
   std::transform(upper.begin(), upper.end(), upper.begin(), (int(*)(int)) std::toupper);

   if (upper == "LINEAR")
      type = gsl_interp_linear;
   else if (upper == "POLYNOMIAL")
      type = gsl_interp_polynomial;
   else if (upper == "CSPLINE")
      type = gsl_interp_cspline;
   else if (upper == "AKIMA")
      type = gsl_interp_akima;
   else {
      std::ostringstream os;
      os<<"interpolationType \""<<upper<<"\" of radial profile \""<<profileName<<"\" unknown";
      FATAL(os.str());
      throw(std::runtime_error("Unknown interpolation type in XML"));
   }

   const unsigned int minNumPoints = gsl_interp_type_min_size (type);

   if (nPoints < minNumPoints) {
      std::ostringstream os;
      os<<"Need at least "<<minNumPoints<<" points for interpolation type \""<<upper<<"\""<<std::endl;
      FATAL(os.str());
      throw(std::runtime_error("Not enough points for spline calculations"));
   }

   //The variable names
   std::vector<std::string> varNames(2*nPoints);
   for (size_t i(0); i < nPoints; ++i) {
      std::ostringstream os;
      os << "radius_"<< std::setfill('0') << std::setw(2) << i;
      varNames[2*i] = os.str();
      os.str("");
      os << "value_"<< std::setfill('0') << std::setw(2) << i;
      varNames[2*i+1] = os.str();
   }

   fvars.add(varNames, profileName, b);

   //Need to store the actual variable names
   fvarNames.resize(nPoints);
   for (size_t i(0); i < nPoints; ++i) {
      fvarNames[i].first = varNames[2*i];
      fvarNames[i].second = varNames[2*i+1];
   }

   updateVariableValues(fvars);

}

SplineRadialProfile::SplineRadialProfile( const utl::Parameters &pars, const std::string &prefix, size_t numPoints, const gsl_interp_type *interpType ) :
   spline(0), acc(0), type(interpType)
{
   profileName = prefix;

   const unsigned int minNumPoints = gsl_interp_type_min_size (type);

   if (numPoints < minNumPoints) {
      std::cerr<<"Need at least "<<minNumPoints<<" points for spline"<<std::endl;
      throw(std::runtime_error("Not enough points for spline calculations"));
   }

   //Create the parameter names
   fvarNames.resize(numPoints);

   std::vector<std::string> names(2*numPoints);

   for (size_t i(0); i < numPoints; ++i) {
      std::ostringstream ii;
      ii << std::setw(2) << std::setfill('0') << i;
      fvarNames[i] = std::pair<std::string,std::string>(prefix+"_radius_"+ii.str(), prefix+"_value_"+ii.str());
      names[2*i] = fvarNames[i].first;
      names[2*i+1] = fvarNames[i].second;
   }

   //Add to the variables object
   fvars.add(names, pars);

   //Update cached values
   updateVariableValues( fvars );
}

SplineRadialProfile::SplineRadialProfile( const utl::Parameters &pars, const std::vector<std::pair<std::string,std::string> > &varNames, const gsl_interp_type *interpType ) :
   spline(0),
   acc(0),
   type(interpType),
   fvarNames(varNames)
{

   const unsigned int minNumPoints = gsl_interp_type_min_size (type);

   if (varNames.size() < minNumPoints) {
      std::cerr<<"Need at least "<<minNumPoints<<" points for spline"<<std::endl;
      throw(std::runtime_error("Not enough points for spline calculations"));
   }

   std::vector<std::string> names(2*varNames.size());

   for (size_t i(0); i < varNames.size(); ++i) {
      names[2*i] = varNames[i].first;
      names[2*i+1] = varNames[i].second;
   }

   //Add to the variables object
   fvars.add(names, pars);

   //Update cached values
   updateVariableValues( fvars );
}

SplineRadialProfile::SplineRadialProfile( const std::vector< std::pair<double,double> > &varValues, const std::vector<std::pair<std::string,std::string> > &varNames, const gsl_interp_type *interpType ) :
   spline(0),
   acc(0),
   type(interpType),
   fvarNames(varNames)
{

   const unsigned int minNumPoints = gsl_interp_type_min_size (type);

   if (varNames.size() != varValues.size()){
      FATAL("Values and names should have same size");
      throw(std::runtime_error("Fix your code"));
   }

   if (varNames.size() < minNumPoints) {
      std::ostringstream os;
      os<<"Need at least "<<minNumPoints<<" points for spline";
      FATAL(os.str());
      throw(std::runtime_error("Not enough points for spline calculations"));
   }

   for (size_t i(0); i < varNames.size(); ++i) {
      fvars.add(varNames[i].first, varValues[i].first);
      fvars.add(varNames[i].second, varValues[i].second);
   }

   //Update cached values
   updateVariableValues( fvars );
}

SplineRadialProfile::~SplineRadialProfile() 
{
   if (spline != 0)
      gsl_spline_free(spline);
}

void SplineRadialProfile::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   addToDOMinternal( node, attributes, "Spline" );
}

void SplineRadialProfile::addToDOMinternal ( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes,
      const std::string &profileType ) const
{
   std::map<std::string, std::string> elements;
   std::ostringstream oss;
   oss << fvarNames.size();
   elements["numberOfPoints"] = oss.str();
   oss.str("");
   if (type == gsl_interp_linear)
      elements["interpolationType"] = "LINEAR";
   else if ( type == gsl_interp_polynomial )
      elements["interpolationType"] = "POLYNOMIAL";
   else if ( type == gsl_interp_cspline )
      elements["interpolationType"] = "CSPLINE";
   else if ( type == gsl_interp_akima )
      elements["interpolationType"] = "AKIMA";

   std::vector<std::string> ids(2*fvarNames.size()), names(2*fvarNames.size());
   for (size_t i(0); i < fvarNames.size(); ++i) {
      oss << "radius_"<< std::setfill('0') << std::setw(2) << i;
      ids[2*i] = oss.str();
      oss.str("");
      oss << "value_"<< std::setfill('0') << std::setw(2) << i;
      ids[2*i+1] = oss.str();
      oss.str("");

      names[2*i] = fvarNames[i].first;
      names[2*i+1] = fvarNames[i].second;
   }

   addToDOMhelper( node, profileType, attributes, elements, names, ids );
}

#ifdef HAVE_OPENCL
//Need to implement the eval function of GSL in OpenCL
//Pack the pars with x, y, and then state

//Start with the index search
std::string SplineRadialProfile::getOpenCLBSearchFun(std::string &name, size_t parIndex) const
{
   name += "bsearch";
   std::ostringstream codeStr;
   codeStr<<"inline uint "<<name<<"(float x, __constant float* pars)\n"
      <<"{\n"
      <<"uint ilo = 0;\n"
      <<"uint ihi = "<<spline->size-1<<";\n"
      <<"while(ihi > ilo + 1) {\n"
      <<"  uint i = (ihi + ilo)/2;\n"
      <<"  if(pars["<<parIndex+2<<"+i] > x)\n"
      <<"    ihi = i;\n"
      <<"  else\n"
      <<"    ilo = i;\n"
      <<"}\n"
      <<"return ilo;\n"
      <<"}\n";
   return codeStr.str();
}

//The linear function, no state in this case
std::string SplineRadialProfile::getOpenCLLinearFun(std::string &name, size_t parIndex) const
{
   std::ostringstream codeStr;
   
   //Start by adding the binary search
   std::string bname = name;
   codeStr<<getOpenCLBSearchFun(bname, parIndex);

   //Then the linear interpolation
   name += "SplineLin";
   codeStr<<"inline float "<<name<<"(float R, __constant float* pars)\n"
      <<"{\n"
      <<"  if ( R < pars["<<parIndex<<"] || R > pars["<<parIndex+1<<"] )\n"
      <<"    return 0;\n"
      <<"  const uint index = "<<bname<<"(R,pars);\n"
      <<"  return pars["<<parIndex+2+spline->size<<"+index] + \n"
      <<"    (R-pars["<<parIndex+2<<"+index])*\n"
      <<"    (pars["<<parIndex+2+spline->size<<"+index+1]-pars["<<parIndex+2+spline->size<<"+index])/\n"
      <<"    (pars["<<parIndex+2<<"+index+1]-pars["<<parIndex+2<<"+index]);\n"
      <<"}\n";
   return codeStr.str();
}

//The CSpline function, pack state->c after y.  It has same size.
std::string SplineRadialProfile::getOpenCLCSplineFun(std::string &name, size_t parIndex) const
{
   std::ostringstream codeStr;
   
   //Start by adding the binary search
   std::string bname = name;
   codeStr<<getOpenCLBSearchFun(bname, parIndex);

   //Then the linear interpolation
   name += "SplineCSP";
   codeStr<<"inline float "<<name<<"(float R, __constant float* pars)\n"
      <<"{\n"
      <<"  if ( R < pars["<<parIndex<<"] || R > pars["<<parIndex+1<<"] )\n"
      <<"    return 0;\n"
      <<"  const uint index = "<<bname<<"(R,pars);\n"
      <<"  return pars["<<parIndex+2+spline->size<<"+index] + \n"
      <<"    (R-pars["<<parIndex+2<<"+index]) * (((\n"
      <<"          (pars["<<parIndex+2+spline->size<<"+index+1]-pars["<<parIndex+2+spline->size<<"+index]) \n"
      <<"          / (pars["<<parIndex+2<<"+index+1]-pars["<<parIndex+2<<"+index]) ) - \n"
      <<"        (pars["<<parIndex+2<<"+index+1]-pars["<<parIndex+2<<"+index]) * \n"
      <<"        (pars["<<parIndex+2+2*spline->size<<"+index+1] + \n"
      <<"          2.0 * pars["<<parIndex+2+2*spline->size<<"+index])\n"
      <<"        / 3.0) + \n"
      <<"      (R-pars["<<parIndex+2<<"+index]) * (pars["<<parIndex+2+2*spline->size<<"+index] + \n"
      <<"        (R-pars["<<parIndex+2<<"+index]) * (\n"
      <<"          (pars["<<parIndex+2+2*spline->size<<"+index+1] - pars["<<parIndex+2+2*spline->size<<"+index])\n"
      <<"          / (3.0 * (pars["<<parIndex+2<<"+index+1]-pars["<<parIndex+2<<"+index])))));\n"
      <<"}\n";
   return codeStr.str();
}

// The akima spline function, pack state->b, c, d after y.
std::string SplineRadialProfile::getOpenCLAkimaFun(std::string &name, size_t parIndex) const
{
   std::ostringstream codeStr;
   
   //Start by adding the binary search
   std::string bname = name;
   codeStr<<getOpenCLBSearchFun(bname, parIndex);

   //Then the linear interpolation
   name += "SplineAkima";
   codeStr<<"inline float "<<name<<"(float R, __constant float* pars)\n"
      <<"{\n"
      <<"  if ( R < pars["<<parIndex<<"] || R > pars["<<parIndex+1<<"] )\n"
      <<"    return 0;\n"
      <<"  const uint index = "<<bname<<"(R,pars);\n"
      <<"  return pars["<<parIndex+2+spline->size<<"+index] + \n"
      <<"    (R - pars["<<parIndex+2<<"+index]) * (pars["<<parIndex+2+2*spline->size<<"+index] + \n"
      <<"      (R - pars["<<parIndex+2<<"+index]) * (pars["<<parIndex+2+3*spline->size<<"+index] + \n"
      <<"        (R - pars["<<parIndex+2<<"+index]) * pars["<<parIndex+2+4*spline->size<<"+index]));\n"
      <<"}\n";
   return codeStr.str();
}

std::string SplineRadialProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   //Select the correct function to return
   if (type == gsl_interp_linear)
      return getOpenCLLinearFun(name, parIndex);
   else if (type == gsl_interp_cspline)
      return getOpenCLCSplineFun(name, parIndex);
   else if (type == gsl_interp_akima)
      return getOpenCLAkimaFun(name, parIndex);
   else {
      ERROR("Implementation not available.");
      throw(std::runtime_error("OPENCL error"));
   }
}

size_t SplineRadialProfile::getOpenCLNPars() const
{
   //The size depends on the type, Rmin and Rmax added.
   if (type == gsl_interp_linear)
      return 2*spline->size+2;
   else if (type == gsl_interp_cspline)
      return 3*spline->size+2;
   else if (type == gsl_interp_akima)
      return 5*spline->size+2;
   else {
      ERROR("Implementation not available.");
      throw(std::runtime_error("OPENCL error"));
   }
}

//From gsl.  Hopefully this code does not change
typedef struct
{
   double * c;
   double * g;
   double * diag;
   double * offdiag;
} cspline_state_t;

typedef struct
{
   double * b;
   double * c;
   double * d;
   double * _m;
} akima_state_t;

std::vector<cl_float> SplineRadialProfile::getOpenCLPars() const
{
   std::vector<cl_float> pars;

   //Start with Rmin and Rmax
   pars.push_back(Rmin);
   pars.push_back(Rmax);

   //Always add the x and y
   for (size_t i = 0; i < spline->size; ++i)
      pars.push_back(spline->x[i]);
   for (size_t i = 0; i < spline->size; ++i)
      pars.push_back(spline->y[i]);

   //The rest depends on the type
   if (type == gsl_interp_linear)
      return pars;
   else if (type == gsl_interp_cspline) {
      const cspline_state_t *s = (const cspline_state_t *) spline->interp->state;
      for (size_t i = 0; i < spline->size; ++i)
         pars.push_back(s->c[i]);
      return pars;
   }
   else if (type == gsl_interp_akima) {
      const akima_state_t *s = (const akima_state_t *) spline->interp->state;
      for (size_t i = 0; i < spline->size; ++i)
         pars.push_back(s->b[i]);
      for (size_t i = 0; i < spline->size; ++i)
         pars.push_back(s->c[i]);
      for (size_t i = 0; i < spline->size; ++i)
         pars.push_back(s->d[i]);
      return pars;
   }
   else {
      ERROR("Implementation not available.");
      throw(std::runtime_error("OPENCL error"));
   }
}

#endif

double SplineRadialProfile::operator() (double R) const 
{
   if ( R < Rmin || R > Rmax ) 
      return 0;

   return gsl_spline_eval(spline, R, acc);
}

void SplineRadialProfile::updateVariableValues( const utl::Variables &vars ) 
{

   //Update the values of the variables
   fvars.setFromVariables( vars );

   //Free the spline if already set
   if (spline != 0)
      gsl_spline_free(spline);

   //Set up the spline
   spline = gsl_spline_alloc(type, fvarNames.size());

   //Get the values for the spline, assert that the order is preserved
   double x[fvarNames.size()], y[fvarNames.size()];

   for (size_t i(0); i < fvarNames.size(); ++i) {
      x[i] = vars[fvarNames[i].first];
      y[i] = vars[fvarNames[i].second];

      if (i > 0 && x[i] <= x[i-1]) {
         std::cerr<<"Radius not strictly ascending in spline profile"<<std::endl;
         throw(std::runtime_error("Could not calculate spline"));
      }
   }

   Rmin = x[0];
   Rmax = x[fvarNames.size()-1];

   gsl_spline_init( spline, x, y, fvarNames.size());

}

void SplineRadialProfile::setRadii(const std::vector<double> &r) {

   if (r.size() != fvarNames.size()) {
      std::cerr<<"Number of radii in spline incorrect"<<std::endl;
      throw(std::runtime_error("Could not update radii values"));
   }

   for (size_t i(0); i < fvarNames.size(); ++i) 
      fvars[fvarNames[i].first] = r[i];

   updateVariableValues(fvars);

}

void SplineRadialProfile::setValues(const std::vector<double> &v) {

   if (v.size() != fvarNames.size()) {
      std::cerr<<"Number of values in spline incorrect"<<std::endl;
      throw(std::runtime_error("Could not update values"));
   }

   for (size_t i(0); i < fvarNames.size(); ++i) 
      fvars[fvarNames[i].second] = v[i];

   updateVariableValues(fvars);

}



static utl::Registry1<RadialProfile, const utl::Branch&>::Registrar<SplineLogRadialProfile> registrarSplineLog("SplineLog");

SplineLogRadialProfile::SplineLogRadialProfile( const utl::Branch &b ) :
   SplineRadialProfile(b)
{

   assert(b.GetAttributes()["type"] == "SplineLog");

   //We need to explicitly update the variables
   updateVariableValues(fvars);

}

SplineLogRadialProfile::SplineLogRadialProfile( const utl::Parameters &pars, const std::string &prefix, size_t numPoints, const gsl_interp_type *interpType ) :
   SplineRadialProfile( pars, prefix, numPoints, interpType )
{

   //We need to explicitly update the variables
   updateVariableValues(fvars);

}


SplineLogRadialProfile::SplineLogRadialProfile( const utl::Parameters &pars, const std::vector<std::pair<std::string,std::string> > &varNames, const gsl_interp_type *interpType ) :
   SplineRadialProfile( pars, varNames, interpType)
{

   //We need to explicitly update the variables
   updateVariableValues(fvars);

}

SplineLogRadialProfile::SplineLogRadialProfile( const std::vector< std::pair<double,double> > &varValues, const std::vector<std::pair<std::string,std::string> > &varNames, const gsl_interp_type *interpType ) :
   SplineRadialProfile(varValues, varNames, interpType)
{

   //We need to explicitly update the variables because virtual functions do not work from constructors and destructors
   updateVariableValues(fvars);

}

void SplineLogRadialProfile::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   addToDOMinternal( node, attributes, "SplineLog" );
}


#ifdef HAVE_OPENCL

std::string SplineLogRadialProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   std::ostringstream codeStr;
   
   //Use the Spline function but wrap it
   std::string sname = name;
   codeStr<<SplineRadialProfile::getOpenCLFunction(sname, parIndex);

   name += "LogSpline";
   codeStr<<"inline float "<<name<<"(float R, __constant float* pars)\n"
      <<"{\n"
      <<"  if ( R < pars["<<parIndex<<"] || R > pars["<<parIndex+1<<"] )\n"
      <<"    return 0;\n"
      <<"  return exp("<<sname<<"(R, pars));\n"
      <<"}\n";
   return codeStr.str();
}

#endif

double SplineLogRadialProfile::operator() (double R) const 
{
   if ( R < Rmin || R > Rmax ) 
      return 0;

   return exp(gsl_spline_eval(spline, R, acc));
}

void SplineLogRadialProfile::updateVariableValues( const utl::Variables &vars ) 
{

   //Update the values of the variables
   fvars.setFromVariables( vars );

   //Free the spline if already set
   if (spline != 0)
      gsl_spline_free(spline);

   //Set up the spline
   spline = gsl_spline_alloc(type, fvarNames.size());

   //Get the values for the spline, assert that the order is preserved
   double x[fvarNames.size()], y[fvarNames.size()];

   for (size_t i(0); i < fvarNames.size(); ++i) {
      x[i] = vars[fvarNames[i].first];
      y[i] = vars[fvarNames[i].second];

      if (y[i] <= 0) {
         std::cerr<<"Need strictly positive values for logarithmic spline"<<std::endl;
         std::cerr<<fvarNames[i].second<<" is <=0!"<<std::endl;
         std::cerr<<vars<<std::endl;
         throw(std::runtime_error("Negative or zero values deteced"));
      }
      y[i] = log(y[i]);

      if (i > 0 && x[i] <= x[i-1]) {
         std::cerr<<"Radius not strictly ascending in spline profile"<<std::endl;
         throw(std::runtime_error("Could not calculate spline"));
      }
   }

   Rmin = x[0];
   Rmax = x[fvarNames.size()-1];

   gsl_spline_init( spline, x, y, fvarNames.size());

}

}
