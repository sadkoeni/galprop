#include "planeprofiles.h"

#include <ErrorLogger.h>
#include <Registry.h>
#include <ReaderErrorReporter.h>

#include <cmath>
#include <stdexcept>

namespace GalacticStructure {

PlaneProfile::PlaneProfile( const utl::Branch &b ) {

   assert(b.GetBranchNameString() == "PlaneProfile");

   utl::Branch nb = b.GetChild("name");
   if ( ! nb) {
      std::ostringstream os;
      os<<"A name is required for plane profiles. Please specify the name element.";
      FATAL(os.str());
      throw(std::runtime_error("No name element found"));
   }
   nb.GetData(profileName);

}

std::unique_ptr<PlaneProfile> PlaneProfile::createProfile(const utl::Branch &b ) {
   assert(b.GetBranchNameString() == "PlaneProfile");

   const std::string name = b.GetAttributes()["type"];

   return utl::Registry1<PlaneProfile, const utl::Branch&>::create(name,b);
}

#ifdef HAVE_OPENCL
   std::string PlaneProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
   {
      ERROR("OpenCL implementation not ready for this class");
      throw(std::runtime_error("Not Implemented"));
   }

   size_t PlaneProfile::getOpenCLNPars() const
   {
      ERROR("OpenCL implementation not ready for this class");
      throw(std::runtime_error("Not Implemented"));
   }

   std::vector<cl_float> PlaneProfile::getOpenCLPars() const
   {
      ERROR("OpenCL implementation not ready for this class");
      throw(std::runtime_error("Not Implemented"));
   }
#endif

static utl::Registry1<PlaneProfile, const utl::Branch&>::Registrar<ExpDiscWithEccentricHoleProfile> registrarExpDiscWithEccentricHole("ExpDiscWithEccentricHole");

  ExpDiscWithEccentricHoleProfile::ExpDiscWithEccentricHoleProfile ( const utl::Branch &b ) : PlaneProfile(b) {

    assert(b.GetAttributes()["type"] == "ExpDiscWithEccentricHole");

    std::vector<std::string> varNames(9);
    varNames[0] = "R0";
    varNames[1] = "Rs";
    varNames[2] = "Rh";
    varNames[3] = "hi";
    varNames[4] = "phiOffset";
    varNames[5] = "eccentricity";
    varNames[6] = "norm";
    varNames[7] = "pitchAngle";
    varNames[8] = "Rmax";

    fvars.add(varNames, profileName, b);

    fR0name = varNames[0];
    fRsname = varNames[1];
    fRhname = varNames[2];
    fhiname = varNames[3];
    fphiOffsetname = varNames[4];
    feccentricityname = varNames[5];
    fnormname = varNames[6];
    fpitchAnglename = varNames[7];
    fRmaxname = varNames[8];

    updateVariableValues( fvars );
   
  }

  ExpDiscWithEccentricHoleProfile::ExpDiscWithEccentricHoleProfile ( double R0, double Rs_, double Rh, double hi_, double phiOffset, double pitchAngle, double holeEccentricity, double norm, double Rmax, const std::string &R0name, const std::string &Rsname, const std::string &Rhname, const std::string &hiname, const std::string &phiOffsetname, const std::string& pitchAnglename, const std::string &holeEccentricityname, const std::string &normname, const std::string& Rmaxname ) :
    fR0name(R0name),
    fRsname(Rsname),
    fRhname(Rhname),
    fhiname(hiname),
    fphiOffsetname(phiOffsetname),
    fpitchAnglename(pitchAnglename),
    feccentricityname(holeEccentricityname),
    fnormname(normname),
    fRmaxname(Rmaxname)
  {
     fvars.add(R0name, R0);
     fvars.add(Rsname, Rs_);
     fvars.add(Rhname, Rh);
     fvars.add(hiname, hi_);
     oneOver_R0 = 1./R0;
     Rs = Rs_;
     oneOver_Rh = 1./Rh;
     hi = hi_;
     fvars.add(phiOffsetname, phiOffset);
     fvars.add(pitchAnglename, pitchAngle);
     fvars.add(holeEccentricityname, holeEccentricity);
     fvars.add(normname, norm);
     fvars.add(Rmaxname, Rmax);
     
     updateVariableValues( fvars );
  }

  ExpDiscWithEccentricHoleProfile::ExpDiscWithEccentricHoleProfile ( const utl::Parameters &pars, const std::string &R0name, const std::string &Rsname, const std::string &Rhname, const std::string &hiname, const std::string &phiOffsetname, const std::string& pitchAnglename, const std::string &holeEccentricityname, const std::string &normname, const std::string& Rmaxname) :
    fR0name(R0name),
    fRsname(Rsname),
    fRhname(Rhname),
    fhiname(hiname),
    fphiOffsetname(phiOffsetname),
    fpitchAnglename(pitchAnglename),
    feccentricityname(holeEccentricityname),
    fnormname(normname),
    fRmaxname(Rmaxname)
  {
    std::vector<std::string> varNames(9);
    varNames[0] = R0name;
    varNames[1] = Rsname;
    varNames[2] = Rhname;
    varNames[3] = hiname;
    varNames[4] = phiOffsetname;
    varNames[5] = holeEccentricityname;
    varNames[6] = normname;
    varNames[7] = pitchAnglename;
    varNames[8] = Rmaxname;
    
    fvars.add(varNames, pars);
    
    updateVariableValues( fvars );
  }
  


void ExpDiscWithEccentricHoleProfile::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   xercesc::DOMDocument *doc = node->getOwnerDocument();

   xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("PlaneProfile").unicodeForm());
   node->appendChild(profileEl);
   profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("ExpDiscWithEccentricHole").unicodeForm());

   for ( auto it = attributes.begin(); it != attributes.end(); ++it )
      profileEl->setAttribute(utl::XStr(it->first).unicodeForm(), utl::XStr(it->second).unicodeForm());

   //Add the name
   {
      xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
      profileEl->appendChild(element);

      xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
      element->appendChild(text);
   }

   std::vector<std::string> varNames(9), varIds(9);

   varIds[0] = "R0";
   varIds[1] = "Rs";
   varIds[2] = "Rh";
   varIds[3] = "hi";
   varIds[4] = "phiOffset";
   varIds[5] = "eccentricity";
   varIds[6] = "norm";
   varIds[7] = "pitchAngle";
   varIds[8] = "Rmax";

   varNames[0] = fR0name;
   varNames[1] = fRsname;
   varNames[2] = fRhname;
   varNames[3] = fhiname;
   varNames[4] = fphiOffsetname;
   varNames[5] = feccentricityname;
   varNames[6] = fnormname;
   varNames[7] = fpitchAnglename;
   varNames[8] = fRmaxname;

   fvars.addToDOM( profileEl, varNames, varIds, profileName );

}


  void ExpDiscWithEccentricHoleProfile::updateVariableValues( const utl::Variables &vars ) {
      
    fvars.setFromVariables(vars);
    
    setR0(fvars[fR0name]);
    setRs(fvars[fRsname]);
    setRh(fvars[fRhname]);
    sethi(fvars[fhiname]);
    setphiOffset(fvars[fphiOffsetname]);
    setpitchAngle(fvars[fpitchAnglename]);
    setholeEccentricity(fvars[feccentricityname]);
    setnorm(fvars[fnormname]);
    setRmax(fvars[fRmaxname]);

  }

  double ExpDiscWithEccentricHoleProfile::operator() (double r, double theta) const
  {
     const auto x = r*std::cos(theta), y = r*std::sin(theta);
     // Coordinates rotated into bar coordinate frame -- need rotation through
     // the phiOffset and bar pitch angle. Rotate first through y-axis, then 
     // about z-axis. The matrix elements for the z-rotation are changed for
     // the off-diagonal elements because of how phiOffset is defined.
     const double xp = cosPhiOffset*cosPitchAngle*x + sinPhiOffset*y;// - cosPhiOffset*sinPitchAngle*z;
     const double yp = -sinPhiOffset*cosPitchAngle*x + cosPhiOffset*y;// + sinPhiOffset*sinPitchAngle*z;
     //const auto xp = cosPhiOffset*x + sinPhiOffset*y;
     //const auto yp = -sinPhiOffset*x + cosPhiOffset*y;

     const auto rhole2 = xp*xp + eccentricity*eccentricity*yp*yp, rholeratio2 = rhole2*oneOver_Rh*oneOver_Rh;

     const auto holefn = 1. - std::exp(-std::pow(rholeratio2, 0.5*hi));

     const auto discfn = (r < Rmax ? std::exp(-(r-Rs)*oneOver_R0) : std::exp(-(r-Rs)*2.)); // Truncated disc following Freudenreich scheme with 0.5 kpc scale-length beyond Rmax

     return norm * discfn * holefn;// * (1 - std::exp(-std::pow(rhole*oneOver_Rh, hi)));

   }

  void ExpDiscWithEccentricHoleProfile::setR0(double r) {
    oneOver_R0 = 1./r;
  }
  
  void ExpDiscWithEccentricHoleProfile::setRs(double r) {
    Rs = r;
  }
  
  void ExpDiscWithEccentricHoleProfile::setRh(double r) {
    oneOver_Rh = 1./r;
  }
  
  void ExpDiscWithEccentricHoleProfile::sethi(double r) {
    hi = r;
  }
  
  void ExpDiscWithEccentricHoleProfile::setphiOffset(double phi) {
    cosPhiOffset = std::cos(phi);
    sinPhiOffset = std::sin(phi);
  }

  void ExpDiscWithEccentricHoleProfile::setpitchAngle(double pa) {

    cosPitchAngle = std::cos(pa);
    sinPitchAngle = std::sin(pa);

  }
  
  void ExpDiscWithEccentricHoleProfile::setholeEccentricity(double q) {
    eccentricity = q;
  }
      
  void ExpDiscWithEccentricHoleProfile::setnorm(double a) {
    norm = a;
  }
  
 void ExpDiscWithEccentricHoleProfile::setRmax(double r) {
    Rmax = r;
  }

static utl::Registry1<PlaneProfile, const utl::Branch&>::Registrar<PlaneProfile1Mode> registrar1Mode("1Mode");

PlaneProfile1Mode::PlaneProfile1Mode( const utl::Branch &b ) :
   PlaneProfile(b)
{

   assert(b.GetAttributes()["type"] == "1Mode");

   std::vector<std::string> varNames(1);
   varNames[0] = "norm";

   fvars.add(varNames, profileName, b);
   fnormName = varNames[0];

   //Loop over the children to find the radial profiles
   for ( utl::Branch rp = b.GetFirstChild();  rp; rp = rp.GetNextSibling() ) {
      if ( rp.GetBranchNameString() == "RadialProfile" ) {
         rProfile = RadialProfile::createProfile( rp );
         break;
      }
   }

   fvars.add(rProfile->getVariables());

   norm = fvars[fnormName];
}

PlaneProfile1Mode::PlaneProfile1Mode ( const utl::Parameters &pars, const std::string &normName, std::unique_ptr<RadialProfile> profile ):
   rProfile(std::move(profile)),
   fnormName(normName)
{

   //Add the variable
   std::vector<std::string> varNames(1);
   varNames[0] = normName;

   fvars.add(varNames, pars);
   fvars.add(rProfile->getVariables());

   norm = fvars[fnormName];

}


void PlaneProfile1Mode::updateVariableValues ( const utl::Variables &vars ) 
{

   rProfile->updateVariableValues(vars);

   fvars.setFromVariables(vars);

   norm = fvars[fnormName];

}

double PlaneProfile1Mode::operator () ( double R, double theta ) const {
   return norm * (*rProfile)(R);
}

void PlaneProfile1Mode::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   xercesc::DOMDocument *doc = node->getOwnerDocument();

   xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("PlaneProfile").unicodeForm());
   node->appendChild(profileEl);
   profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("1Mode").unicodeForm());

   for ( auto it = attributes.begin(); it != attributes.end(); ++it )
      profileEl->setAttribute(utl::XStr(it->first).unicodeForm(), utl::XStr(it->second).unicodeForm());

   //Add the name
   {
      xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
      profileEl->appendChild(element);

      xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
      element->appendChild(text);
   }

   std::vector<std::string> varNames(1), varIds(1);

   varIds[0] = "norm";

   varNames[0] = fnormName;

   fvars.addToDOM( profileEl, varNames, varIds, profileName );

   std::map<std::string,std::string> emptyMap;
   rProfile->addToDOM( profileEl, emptyMap );

}

#ifdef HAVE_OPENCL
std::string PlaneProfile1Mode::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   //Update the name
   name = "PlaneProfile1Mode_"+profileName+"_"+name;

   //The name for the radial profile
   std::string rProfileName = name;

   //Use an ostringstream object to store the function name, start with the external function
   //Also keep track of parIndex, store our parameters at the back
   std::ostringstream codeStr;

   codeStr<<rProfile->getOpenCLFunction(rProfileName, parIndex);
   parIndex += rProfile->getOpenCLNPars();

   codeStr<<"inline float "<<name<<"( float R, float theta, __constant float *pars )\n"
      <<"{\n"
      <<"return pars["<<parIndex<<"] * "<<rProfileName<<"(R,pars);\n"
      <<"}\n";

   return codeStr.str();
}

size_t PlaneProfile1Mode::getOpenCLNPars() const
{
   return rProfile->getOpenCLNPars() + 1;
}

std::vector<cl_float> PlaneProfile1Mode::getOpenCLPars() const
{
   std::vector<cl_float> pars = rProfile->getOpenCLPars();
   pars.push_back(norm);
   return pars;
}

#endif


static utl::Registry1<PlaneProfile, const utl::Branch&>::Registrar<PlaneProfile2Mode> registrar2Mode("2Mode");

PlaneProfile2Mode::PlaneProfile2Mode( const utl::Branch &b ) :
   PlaneProfile(b)
{

   assert(b.GetAttributes()["type"] == "2Mode");

   std::vector<std::string> varNames(2);
   varNames[0] = "norm1";
   varNames[1] = "norm2";

   fvars.add(varNames, profileName, b);
   fnorm1Name = varNames[0];
   fnorm2Name = varNames[1];

   //Loop over the children to find the radial profiles
   for ( utl::Branch rp = b.GetFirstChild();  rp; rp = rp.GetNextSibling() ) {

      if ( rp.GetBranchNameString() == "RadialProfile" ) {

         //Look for the intent attribute
         std::string upper = rp.GetAttributes()["intent"];
         // explicit cast needed to resolve ambiguity
         std::transform(upper.begin(), upper.end(), upper.begin(), (int(*)(int)) std::toupper);

         if ( upper == "NORM1" ) {
            fnorm1Profile = RadialProfile::createProfile( rp );
         } else if ( upper == "NORM2" ) {
            fnorm2Profile = RadialProfile::createProfile( rp );
         } else if ( upper == "ZERO2" ) {
            fzero2Profile = RadialProfile::createProfile( rp );
         }
      }
   }

   //Make sure they are all there
   if ( fnorm1Profile.get() == nullptr ) {
      FATAL("Missing normalization radial profile for mode 1 in mode2 profile");
      throw(std::runtime_error("Radial profile missing"));
   }
   if ( fnorm2Profile.get() == nullptr ) {
      FATAL("Missing normalization radial profile for mode 2 in mode2 profile");
      throw(std::runtime_error("Radial profile missing"));
   }
   if ( fzero2Profile.get() == nullptr ) {
      FATAL("Missing zero point radial profile for mode 2 in mode2 profile");
      throw(std::runtime_error("Radial profile missing"));
   }


   fvars.add(fnorm1Profile->getVariables());
   fvars.add(fnorm2Profile->getVariables());
   fvars.add(fzero2Profile->getVariables());

   norm1 = fvars[fnorm1Name];
   norm2 = fvars[fnorm2Name];
}

PlaneProfile2Mode::PlaneProfile2Mode ( const utl::Parameters &pars, 
      const std::string &norm1Name, 
      const std::string &norm2Name, 
      std::unique_ptr<RadialProfile> norm1profile,
      std::unique_ptr<RadialProfile> norm2profile, 
      std::unique_ptr<RadialProfile> zero2profile ):
   fnorm1Profile(std::move(norm1profile)),
   fnorm2Profile(std::move(norm2profile)),
   fzero2Profile(std::move(zero2profile)),
   fnorm1Name(norm1Name),
   fnorm2Name(norm2Name)
{

   //Add the variable
   std::vector<std::string> varNames(2);
   varNames[0] = norm1Name;
   varNames[1] = norm2Name;

   fvars.add(varNames, pars);
   fvars.add(fnorm1Profile->getVariables());
   fvars.add(fnorm2Profile->getVariables());
   fvars.add(fzero2Profile->getVariables());

   norm1 = fvars[fnorm1Name];
   norm2 = fvars[fnorm2Name];

}

void PlaneProfile2Mode::updateVariableValues ( const utl::Variables &vars ) 
{

   fnorm1Profile->updateVariableValues(vars);
   fnorm2Profile->updateVariableValues(vars);
   fzero2Profile->updateVariableValues(vars);

   fvars.setFromVariables(vars);

   norm1 = fvars[fnorm1Name];
   norm2 = fvars[fnorm2Name];

}

double PlaneProfile2Mode::operator () ( double R, double theta ) const {
   return norm1 * (*fnorm1Profile)(R) + norm2 * (*fnorm2Profile)(R) * cos( theta - (*fzero2Profile)(R) );
}


void PlaneProfile2Mode::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   xercesc::DOMDocument *doc = node->getOwnerDocument();

   xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("PlaneProfile").unicodeForm());
   node->appendChild(profileEl);
   profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("2Mode").unicodeForm());

   for ( auto it = attributes.begin(); it != attributes.end(); ++it )
      profileEl->setAttribute(utl::XStr(it->first).unicodeForm(), utl::XStr(it->second).unicodeForm());

   //Add the name
   {
      xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
      profileEl->appendChild(element);

      xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
      element->appendChild(text);
   }

   std::vector<std::string> varNames(2), varIds(2);

   varIds[0] = "norm1";
   varIds[1] = "norm2";

   varNames[0] = fnorm1Name;
   varNames[1] = fnorm2Name;

   fvars.addToDOM( profileEl, varNames, varIds, profileName );

   std::map<std::string,std::string> attributeMap;
   attributeMap["intent"] = "Norm1";
   fnorm1Profile->addToDOM( profileEl, attributeMap );
   attributeMap["intent"] = "Norm2";
   fnorm2Profile->addToDOM( profileEl, attributeMap );
   attributeMap["intent"] = "Zero2";
   fzero2Profile->addToDOM( profileEl, attributeMap );

}


static utl::Registry1<PlaneProfile, const utl::Branch&>::Registrar<PlaneProfile3Mode> registrar3Mode("3Mode");

PlaneProfile3Mode::PlaneProfile3Mode( const utl::Branch &b ) :
   PlaneProfile(b)
{

   assert(b.GetAttributes()["type"] == "3Mode");

   std::vector<std::string> varNames(3);
   varNames[0] = "norm1";
   varNames[1] = "norm2";
   varNames[2] = "norm3";

   fvars.add(varNames, profileName, b);
   fnorm1Name = varNames[0];
   fnorm2Name = varNames[1];
   fnorm3Name = varNames[2];

   //Loop over the children to find the radial profiles
   for ( utl::Branch rp = b.GetFirstChild();  rp; rp = rp.GetNextSibling() ) {

      if ( rp.GetBranchNameString() == "RadialProfile" ) {

         //Look for the intent attribute
         std::string upper = rp.GetAttributes()["intent"];
         // explicit cast needed to resolve ambiguity
         std::transform(upper.begin(), upper.end(), upper.begin(), (int(*)(int)) std::toupper);

         if ( upper == "NORM1" ) {
            fnorm1Profile = RadialProfile::createProfile( rp );
         } else if ( upper == "NORM2" ) {
            fnorm2Profile = RadialProfile::createProfile( rp );
         } else if ( upper == "ZERO2" ) {
            fzero2Profile = RadialProfile::createProfile( rp );
         } else if ( upper == "NORM3" ) {
            fnorm3Profile = RadialProfile::createProfile( rp );
         } else if ( upper == "ZERO3" ) {
            fzero3Profile = RadialProfile::createProfile( rp );
         }
      }
   }

   //Make sure they are all there
   if ( fnorm1Profile.get() == nullptr ) {
      FATAL("Missing normalization radial profile for mode 1 in mode3 profile");
      throw(std::runtime_error("Radial profile missing"));
   }
   if ( fnorm2Profile.get() == nullptr ) {
      FATAL("Missing normalization radial profile for mode 2 in mode3 profile");
      throw(std::runtime_error("Radial profile missing"));
   }
   if ( fzero2Profile.get() == nullptr ) {
      FATAL("Missing zero point radial profile for mode 2 in mode3 profile");
      throw(std::runtime_error("Radial profile missing"));
   }
   if ( fnorm3Profile.get() == nullptr ) {
      FATAL("Missing normalization radial profile for mode 3 in mode3 profile");
      throw(std::runtime_error("Radial profile missing"));
   }
   if ( fzero3Profile.get() == nullptr ) {
      FATAL("Missing zero point radial profile for mode 3 in mode3 profile");
      throw(std::runtime_error("Radial profile missing"));
   }

   fvars.add(fnorm1Profile->getVariables());
   fvars.add(fnorm2Profile->getVariables());
   fvars.add(fzero2Profile->getVariables());
   fvars.add(fnorm3Profile->getVariables());
   fvars.add(fzero3Profile->getVariables());

   norm1 = fvars[fnorm1Name];
   norm2 = fvars[fnorm2Name];
   norm3 = fvars[fnorm3Name];
}

PlaneProfile3Mode::PlaneProfile3Mode ( const utl::Parameters &pars, 
      const std::string &norm1Name, 
      const std::string &norm2Name, 
      const std::string &norm3Name, 
      std::unique_ptr<RadialProfile> norm1profile,
      std::unique_ptr<RadialProfile> norm2profile,
      std::unique_ptr<RadialProfile> zero2profile,
      std::unique_ptr<RadialProfile> norm3profile,
      std::unique_ptr<RadialProfile> zero3profile ):
   fnorm1Profile(std::move(norm1profile)),
   fnorm2Profile(std::move(norm2profile)),
   fzero2Profile(std::move(zero2profile)),
   fnorm3Profile(std::move(norm3profile)),
   fzero3Profile(std::move(zero3profile)),
   fnorm1Name(norm1Name),
   fnorm2Name(norm2Name),
   fnorm3Name(norm3Name)
{

   //Add the variable
   std::vector<std::string> varNames(3);
   varNames[0] = norm1Name;
   varNames[1] = norm2Name;
   varNames[2] = norm3Name;

   fvars.add(varNames, pars);
   fvars.add(fnorm1Profile->getVariables());
   fvars.add(fnorm2Profile->getVariables());
   fvars.add(fzero2Profile->getVariables());
   fvars.add(fnorm3Profile->getVariables());
   fvars.add(fzero3Profile->getVariables());

   norm1 = fvars[fnorm1Name];
   norm2 = fvars[fnorm2Name];
   norm3 = fvars[fnorm3Name];

}

void PlaneProfile3Mode::updateVariableValues ( const utl::Variables &vars ) 
{
   fnorm1Profile->updateVariableValues(vars);
   fnorm2Profile->updateVariableValues(vars);
   fzero2Profile->updateVariableValues(vars);
   fnorm3Profile->updateVariableValues(vars);
   fzero3Profile->updateVariableValues(vars);

   fvars.setFromVariables(vars);

   norm1 = fvars[fnorm1Name];
   norm2 = fvars[fnorm2Name];
   norm3 = fvars[fnorm3Name];

}

double PlaneProfile3Mode::operator () ( double R, double theta ) const {
   return norm1 * (*fnorm1Profile)(R) + 
      norm2 * (*fnorm2Profile)(R) * cos( theta - (*fzero2Profile)(R) ) +
      norm3 * (*fnorm3Profile)(R) * cos( 2*theta - (*fzero3Profile)(R) );
}

void PlaneProfile3Mode::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   xercesc::DOMDocument *doc = node->getOwnerDocument();

   xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("PlaneProfile").unicodeForm());
   node->appendChild(profileEl);
   profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("3Mode").unicodeForm());

   for ( auto it = attributes.begin(); it != attributes.end(); ++it )
      profileEl->setAttribute(utl::XStr(it->first).unicodeForm(), utl::XStr(it->second).unicodeForm());

   //Add the name
   {
      xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
      profileEl->appendChild(element);

      xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
      element->appendChild(text);
   }

   std::vector<std::string> varNames(3), varIds(3);

   varIds[0] = "norm1";
   varIds[1] = "norm2";
   varIds[2] = "norm2";

   varNames[0] = fnorm1Name;
   varNames[1] = fnorm2Name;
   varNames[2] = fnorm3Name;

   fvars.addToDOM( profileEl, varNames, varIds, profileName );

   std::map<std::string,std::string> attributeMap;
   attributeMap["intent"] = "Norm1";
   fnorm1Profile->addToDOM( profileEl, attributeMap );
   attributeMap["intent"] = "Norm2";
   fnorm2Profile->addToDOM( profileEl, attributeMap );
   attributeMap["intent"] = "Norm3";
   fnorm3Profile->addToDOM( profileEl, attributeMap );
   attributeMap["intent"] = "Zero2";
   fzero2Profile->addToDOM( profileEl, attributeMap );
   attributeMap["intent"] = "Zero3";
   fzero3Profile->addToDOM( profileEl, attributeMap );

}



static utl::Registry1<PlaneProfile, const utl::Branch&>::Registrar<PlaneProfileNMode> registrarNMode("NMode");

PlaneProfileNMode::PlaneProfileNMode( const utl::Branch &b ) :
   PlaneProfile(b)
{

   assert(b.GetAttributes()["type"] == "NMode");

   //Read the number of modes
   utl::Branch nmb = b.GetChild("numberOfModes");
   size_t nModes(0);
   nmb.GetData(nModes);

   if (nModes == 0) {
      FATAL("Number of modes set to 0 in NMODE profile");
      throw(std::runtime_error("Need at least 1 mode"));
   }

   fnormNames.resize(nModes);
   for (size_t i(0); i < nModes; ++i) {
      std::ostringstream os;
      os << "norm"<<(i+1);
      fnormNames[i] = os.str();
   }

   fzeroNames.resize(nModes-1);
   for (size_t i(0); i < fzeroNames.size(); ++i) {
      std::ostringstream os;
      os << "zero"<<(i+2);
      fzeroNames[i] = os.str();
   }

   fvars.add(fnormNames, profileName, b);
   fvars.add(fzeroNames, profileName, b);

   fnormProfiles.reserve(nModes);
   fzeroProfiles.reserve(nModes-1);
   //Find the correct radial profiles
   auto findRadialProfileIntent = []( const utl::Branch &tb, std::string intent )->utl::Branch {
      //Do a case insensitive comparison
      std::transform(intent.begin(), intent.end(), intent.begin(), (int(*)(int)) std::toupper);
      for ( utl::Branch rp = tb.GetFirstChild();  rp; rp = rp.GetNextSibling() ) {
         if ( rp.GetBranchNameString() == "RadialProfile" ) {
            std::string upper = rp.GetAttributes()["intent"];
            std::transform(upper.begin(), upper.end(), upper.begin(), (int(*)(int)) std::toupper);
            if (upper == intent)
               return rp;
         }
      }
      return utl::Branch();
   };

   std::ostringstream os;
   for (size_t nn(1); nn <= nModes; ++nn) { 

      os.str("");
      os<<"NORM"<<nn;
      utl::Branch rp = findRadialProfileIntent(b, os.str());

      if ( rp )
         fnormProfiles.push_back(RadialProfile::createProfile( rp ) );
      else {
         os.str("");
         os << "Missing normalization radial profile for mode "<<nn<<" in modeN profile";
         FATAL(os.str());
         throw(std::runtime_error("Radial profile missing"));
      }

      if ( nn > 1 ) {
         os.str("");
         os<<"ZERO"<<nn;
         rp = findRadialProfileIntent(b, os.str());

         if ( rp )
            fzeroProfiles.push_back(RadialProfile::createProfile( rp ) );
         else {
            os.str("");
            os << "Missing zero point radial profile for mode "<<nn<<" in modeN profile";
            FATAL(os.str());
            throw(std::runtime_error("Radial profile missing"));
         }
      }
   }

   for (size_t i(0); i < fnormProfiles.size(); ++i)
      fvars.add(fnormProfiles[i]->getVariables());
   for (size_t i(0); i < fzeroProfiles.size(); ++i)
      fvars.add(fzeroProfiles[i]->getVariables());

   norm.resize(fnormNames.size());
   for (size_t i(0); i < fnormNames.size(); ++i)
      norm[i] = fvars[fnormNames[i]];
   zero.resize(fzeroNames.size());
   for (size_t i(0); i < fzeroNames.size(); ++i)
      zero[i] = fvars[fzeroNames[i]];
}

PlaneProfileNMode::PlaneProfileNMode ( const utl::Parameters &pars, 
      const std::vector<std::string> &normNames,
      const std::vector<std::string> &zeroNames,
      std::vector<std::unique_ptr<RadialProfile> > &&normProfiles,
      std::vector<std::unique_ptr<RadialProfile> > &&zeroProfiles) :
   fnormProfiles(std::move(normProfiles)),
   fzeroProfiles(std::move(zeroProfiles)),
   fnormNames(normNames),
   fzeroNames(zeroNames)
{

   //Make sure length matches
   if (fnormNames.size() != fnormProfiles.size() ) {
      std::cerr<<"Number of normaliztion variable names in NMode PlaneProfile does not match profile number"<<std::endl;
      throw(std::runtime_error("The number of normProfiles should be equal to the norm variable names"));
   }
   if (fzeroNames.size() != fzeroProfiles.size() ) {
      std::cerr<<"Number of zero point variable names in NMode PlaneProfile does not match profile number"<<std::endl;
      throw(std::runtime_error("The number of zeroProfiles should be equal to the zero variable names"));
   }

   if (fnormProfiles.size() != fzeroProfiles.size() + 1) {
      std::cerr<<"Number of profiles in NMode PlaneProfile does not match"<<std::endl;
      throw(std::runtime_error("The number of normProfiles should be larger than zeroProfiles by exactly 1"));
   }

   //Add the variable
   fvars.add(fnormNames, pars);
   fvars.add(fzeroNames, pars);

   for (size_t i(0); i < fnormProfiles.size(); ++i)
      fvars.add(fnormProfiles[i]->getVariables());
   for (size_t i(0); i < fzeroProfiles.size(); ++i)
      fvars.add(fzeroProfiles[i]->getVariables());

   norm.resize(fnormNames.size());
   for (size_t i(0); i < fnormNames.size(); ++i)
      norm[i] = fvars[fnormNames[i]];
   zero.resize(fzeroNames.size());
   for (size_t i(0); i < fzeroNames.size(); ++i)
      zero[i] = fvars[fzeroNames[i]];
      
}

PlaneProfileNMode::PlaneProfileNMode ( const utl::Parameters &pars, 
      const std::string &prefix,
      std::vector<std::unique_ptr<RadialProfile> > &&normProfiles,
      std::vector<std::unique_ptr<RadialProfile> > &&zeroProfiles) :
   fnormProfiles(std::move(normProfiles)),
   fzeroProfiles(std::move(zeroProfiles))
{

   //Make sure length matches
   if (fnormProfiles.size() != fzeroProfiles.size() + 1) {
      std::cerr<<"Number of profiles in NMode PlaneProfile does not match"<<std::endl;
      throw(std::runtime_error("The number of normProfiles should be larger than zeroProfiles by exactly 1"));
   }

   //Add the variable
   fnormNames.resize(fnormProfiles.size());

   for (size_t i(0); i < fnormNames.size(); ++i) {
      std::ostringstream ii;
      ii << i+1;
      fnormNames[i] = prefix + "_norm_" + ii.str();
   }

   fvars.add(fnormNames, pars);

   fzeroNames.resize(fzeroProfiles.size());

   for (size_t i(0); i < fzeroNames.size(); ++i) {
      std::ostringstream ii;
      ii << i+2;
      fzeroNames[i] = prefix + "_zero_" + ii.str();
   }

   fvars.add(fzeroNames, pars);

   for (size_t i(0); i < fnormProfiles.size(); ++i)
      fvars.add(fnormProfiles[i]->getVariables());
   for (size_t i(0); i < fzeroProfiles.size(); ++i)
      fvars.add(fzeroProfiles[i]->getVariables());

   norm.resize(fnormNames.size());
   for (size_t i(0); i < fnormNames.size(); ++i)
      norm[i] = fvars[fnormNames[i]];
      
   zero.resize(fzeroNames.size());
   for (size_t i(0); i < fzeroNames.size(); ++i)
      zero[i] = fvars[fzeroNames[i]];
      
}

void PlaneProfileNMode::updateVariableValues ( const utl::Variables &vars ) 
{

   for (size_t i(0); i < fnormProfiles.size(); ++i)
      fnormProfiles[i]->updateVariableValues(vars);
   for (size_t i(0); i < fzeroProfiles.size(); ++i)
      fzeroProfiles[i]->updateVariableValues(vars);

   fvars.setFromVariables(vars);

   for (size_t i(0); i < fnormNames.size(); ++i)
      norm[i] = fvars[fnormNames[i]];
   for (size_t i(0); i < fzeroNames.size(); ++i)
      zero[i] = fvars[fzeroNames[i]];
      
}

double PlaneProfileNMode::operator () ( double R, double theta ) const {

   double output = norm[0] * (*fnormProfiles[0])(R);

   for (size_t i(1); i < fnormProfiles.size(); ++i) {
      const double it = i*theta;
      const double zp = zero[i-1]*(*fzeroProfiles[i-1])(R);
      output += norm[i] * (*fnormProfiles[i])(R) * (cos(it)*cos(zp)+sin(it)*sin(zp));
   }

   return output;
}

void PlaneProfileNMode::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   xercesc::DOMDocument *doc = node->getOwnerDocument();

   xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("PlaneProfile").unicodeForm());
   node->appendChild(profileEl);
   profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("NMode").unicodeForm());

   for ( auto it = attributes.begin(); it != attributes.end(); ++it )
      profileEl->setAttribute(utl::XStr(it->first).unicodeForm(), utl::XStr(it->second).unicodeForm());

   //Add the name
   {
      xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
      profileEl->appendChild(element);

      xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
      element->appendChild(text);
   }

   xercesc::DOMElement* element = doc->createElement(utl::XStr("numberOfModes").unicodeForm());
   profileEl->appendChild(element);

   std::ostringstream oss;
   oss << fnormProfiles.size();
   xercesc::DOMText* text = doc->createTextNode(utl::XStr(oss.str()).unicodeForm());
   element->appendChild(text);
   oss.str("");

   std::vector<std::string> varNames(fnormNames.size()+fzeroNames.size()), varIds(fnormNames.size()+fzeroNames.size());

   for (size_t i(0); i < fnormNames.size(); ++i) {
      oss << "norm"<<(i+1);
      varIds[i] = oss.str();
      oss.str("");
      varNames[i] = fnormNames[i];
   }

   for (size_t i(0); i < fzeroNames.size(); ++i) {
      oss << "zero"<<(i+2);
      varIds[fnormNames.size()+i] = oss.str();
      oss.str("");
      varNames[fnormNames.size()+i] = fzeroNames[i];
   }

   fvars.addToDOM( profileEl, varNames, varIds, profileName );

   std::map<std::string,std::string> attributeMap;
   for (size_t i = 0; i < fnormProfiles.size(); ++i) {
      oss << "Norm" << i+1;
      attributeMap["intent"] = oss.str();
      fnormProfiles[i]->addToDOM( profileEl, attributeMap );
      oss.str("");
   }
   for (size_t i = 0; i < fzeroProfiles.size(); ++i) {
      oss << "Zero" << i+2;
      attributeMap["intent"] = oss.str();
      fzeroProfiles[i]->addToDOM( profileEl, attributeMap );
      oss.str("");
   }

}


#ifdef HAVE_OPENCL
std::string PlaneProfileNMode::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   //Update the name
   name = "PlaneProfileNMode_"+profileName+"_"+name;

   //The name for the profiles
   std::vector<std::string> normNames(fnormProfiles.size()), zeroNames(fzeroProfiles.size());

   //Use an ostringstream object to store the function name, start with the external function
   //Also keep track of parIndex, store our parameters at the back
   std::ostringstream codeStr;

   for (size_t i(0); i < fnormProfiles.size(); ++i){
      std::ostringstream buf;
      buf<<name<<"norm"<<i;
      normNames[i] = buf.str();
      codeStr<<fnormProfiles[i]->getOpenCLFunction(normNames[i], parIndex);
      parIndex += fnormProfiles[i]->getOpenCLNPars();
   }

   for (size_t i(0); i < fzeroProfiles.size(); ++i){
      std::ostringstream buf;
      buf<<name<<"zero"<<i;
      zeroNames[i] = buf.str();
      codeStr<<fzeroProfiles[i]->getOpenCLFunction(zeroNames[i], parIndex);
      parIndex += fzeroProfiles[i]->getOpenCLNPars();
   }

   codeStr<<"inline float "<<name<<"( float R, float theta, __constant float *pars )\n"
      <<"{\n"
      <<"return pars["<<parIndex<<"] * "<<normNames[0]<<"(R,pars)";
   for (size_t i(1); i < fnormProfiles.size(); ++i)
      codeStr<<" + pars["<<parIndex+i<<"] * "<<normNames[i]<<"(R,pars) * (cos("<<i<<"*theta)*cos(pars["<<parIndex+fnormProfiles.size()+i-1<<"]*"<<zeroNames[i-1]<<"(R,pars)) + sin("<<i<<"*theta)*sin(pars["<<parIndex+fnormProfiles.size()+i-1<<"]*"<<zeroNames[i-1]<<"(R,pars)))";
   codeStr<<";\n"
      <<"}\n";

   return codeStr.str();
}

size_t PlaneProfileNMode::getOpenCLNPars() const
{
   size_t npars(0);
   for (const auto & p: fnormProfiles)
      npars += p->getOpenCLNPars();
   for (const auto & p: fzeroProfiles)
      npars += p->getOpenCLNPars();
   return npars + fnormProfiles.size() + fzeroProfiles.size();
}

std::vector<cl_float> PlaneProfileNMode::getOpenCLPars() const
{
   std::vector<cl_float> pars; 
   for (const auto & p : fnormProfiles) {
      std::vector<cl_float> tmpPars = p->getOpenCLPars();
      pars.insert(pars.end(), tmpPars.begin(), tmpPars.end());
   }
   for (const auto & p : fzeroProfiles) {
      std::vector<cl_float> tmpPars = p->getOpenCLPars();
      pars.insert(pars.end(), tmpPars.begin(), tmpPars.end());
   }
   for (const auto & n : norm)
      pars.push_back(n);
   for (const auto & z : zero)
      pars.push_back(z);
   return pars;
}

#endif



static utl::Registry1<PlaneProfile, const utl::Branch&>::Registrar<PlaneProfileArm> registrarArm("Arm");

PlaneProfileArm::PlaneProfileArm ( const utl::Branch &b ) :
   PlaneProfile(b),
   norm(0)
{
   assert(b.GetAttributes()["type"] == "Arm");

   //The variables for the arm
   fvarNames.resize(5);
   fvarNames[0] = "norm";
   fvarNames[1] = "a";
   fvarNames[2] = "rMin";
   fvarNames[3] = "phiMin";
   fvarNames[4] = "width";

   fvars.add(fvarNames, profileName, b);

   fvarNames.resize(6);
   use_rMax = false;
   try {
      std::vector<std::string> varNames(1);
      varNames[0] = "rMax";
      fvars.add(varNames, profileName, b);
      fvarNames[5] = varNames[0];
      use_rMax = true;
   } catch (utl::Variables::VariableError) {}

   try {
      std::vector<std::string> varNames(1);
      varNames[0] = "phiExtent";
      fvars.add(varNames, profileName, b);
      fvarNames[5] = varNames[0];
      if (use_rMax) {
         ERROR("Both phiExtent and rMax given, only one can be used.");
         throw(std::invalid_argument("Invalid variables in XML files"));
      }
   } catch (utl::Variables::VariableError &e) {
      if (! use_rMax)
         throw(e);
   }

   //The arm function
   utl::Branch tb = b.GetChild("armfunction");
   if ( ! tb ) {
      std::ostringstream os;
      os << "armfunction element missing in PlaneProfileArm \""<<profileName<<"\"";
      FATAL(os.str());
      throw(std::runtime_error("armfunction element is required"));
   }

   std::string armFName;
   tb.GetData(armFName);

   auto scale = utl::Registry0<ScaleFunction>::create(armFName);

   //Create an empty arm function, it will be updated with the real values later
   arm = std::unique_ptr<ArmFunction>(new ArmFunction(1.0, 1.0, 0.0, 0.0, 1.0, std::move(scale)));

   //Loop over the children to find the radial profiles
   for ( utl::Branch rp = b.GetFirstChild();  rp; rp = rp.GetNextSibling() ) {
      if ( rp.GetBranchNameString() == "RadialProfile" ) {
         fradProfile = RadialProfile::createProfile( rp );
         break;
      }
   }

   fvars.add(fradProfile->getVariables());

   updateVariableValues( fvars );
   
}

PlaneProfileArm::PlaneProfileArm ( const utl::Parameters& pars, 
      std::unique_ptr<ArmFunction> armFunction,
      std::unique_ptr<RadialProfile> radProfile, 
      const std::string &prefix,
      bool rMax) :
   norm(0),
   arm(std::move(armFunction)),
   fradProfile(std::move(radProfile)),
   use_rMax(rMax)
{

   //Create the variables vector
   fvarNames.resize(6);
   fvarNames[0] = prefix + "_norm";
   fvarNames[1] = prefix + "_a";
   fvarNames[2] = prefix + "_rMin";
   fvarNames[3] = prefix + "_phiMin";
   fvarNames[4] = prefix + "_width";
   if (rMax)
      fvarNames[5] = prefix + "_rMax";
   else
      fvarNames[5] = prefix + "_phiExtent";

   fvars.add(fvarNames, pars);

   fvars.add(fradProfile->getVariables());

   updateVariableValues( fvars );
}

PlaneProfileArm::PlaneProfileArm ( const utl::Parameters& pars, 
      std::unique_ptr<ArmFunction> armFunction,
      std::unique_ptr<RadialProfile> radProfile, 
      const std::vector<std::string> &varNames,
      bool rMax):
   norm(0),
   arm(std::move(armFunction)),
   fradProfile(std::move(radProfile)),
   fvarNames(varNames),
   use_rMax(rMax)
{

   if (fvarNames.size() != 6) {
      std::cerr<<"Should have exactly 6 variables for arm profile."<<std::endl;
      throw(std::runtime_error("Number of variables incorrect in model."));
   }

   fvars.add(fvarNames, pars);

   fvars.add(fradProfile->getVariables());

   updateVariableValues( fvars );

}

void PlaneProfileArm::updateVariableValues( const utl::Variables &vars ) {

   fradProfile->updateVariableValues(vars);

   fvars.setFromVariables(vars);

   norm = vars[fvarNames[0]];
   const double a = vars[fvarNames[1]];
   const double rMin = vars[fvarNames[2]];
   const double phiMin = vars[fvarNames[3]];
   const double width = vars[fvarNames[4]];

   arm->setA(a);
   arm->setRMin(rMin);
   arm->setPhiMin(phiMin);
   arm->setWidth(width);
   if (use_rMax) {
      const double rMax = vars[fvarNames[5]];
      arm->setPhiExtentFromRMax(rMax);
   } else {
      const double phiExtent = vars[fvarNames[5]];
      arm->setPhiExtent(phiExtent);
   }

}

double PlaneProfileArm::operator() ( double radius, double theta ) const {
   return norm * (*fradProfile)(radius) * (*arm)(radius,theta);
}

void PlaneProfileArm::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   xercesc::DOMDocument *doc = node->getOwnerDocument();

   xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("PlaneProfile").unicodeForm());
   node->appendChild(profileEl);
   profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("Arm").unicodeForm());

   for ( auto it = attributes.begin(); it != attributes.end(); ++it )
      profileEl->setAttribute(utl::XStr(it->first).unicodeForm(), utl::XStr(it->second).unicodeForm());

   //Add the name
   {
      xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
      profileEl->appendChild(element);

      xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
      element->appendChild(text);
   }

   xercesc::DOMElement* element = doc->createElement(utl::XStr("armfunction").unicodeForm());
   profileEl->appendChild(element);

   xercesc::DOMText* text = doc->createTextNode(utl::XStr(arm->scaleName()).unicodeForm());
   element->appendChild(text);

   std::vector<std::string> varNames(6), varIds(6);

   varIds[0] = "norm";
   varIds[1] = "a";
   varIds[2] = "rMin";
   varIds[3] = "phiMin";
   varIds[4] = "width";
   if (use_rMax)
      varIds[5] = "rMax";
   else
      varIds[5] = "phiExtent";

   fvars.addToDOM( profileEl, fvarNames, varIds, profileName );

   std::map<std::string,std::string> attributeMap;
   fradProfile->addToDOM( profileEl, attributeMap );
}


#ifdef HAVE_OPENCL
std::string PlaneProfileArm::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   //Update the name
   name = "PlaneProfileArm_"+profileName+"_"+name;

   //The name for the radial profile
   std::string fradProfileName = name;
   std::string armProfileName = name;

   //Use an ostringstream object to store the function name, start with the external function
   //Also keep track of parIndex, store our parameters at the back
   std::ostringstream codeStr;

   codeStr<<fradProfile->getOpenCLFunction(fradProfileName, parIndex);
   parIndex += fradProfile->getOpenCLNPars();
   codeStr<<arm->getOpenCLFunction(armProfileName, parIndex);
   parIndex += arm->getOpenCLNPars();

   codeStr<<"inline float "<<name<<"( float R, float theta, __constant float *pars )\n"
      <<"{\n"
      <<"return pars["<<parIndex<<"] * "<<fradProfileName<<"(R,pars) * "<<armProfileName<<"(R,theta,pars);\n"
      <<"}\n";

   return codeStr.str();
}

size_t PlaneProfileArm::getOpenCLNPars() const
{
   return fradProfile->getOpenCLNPars() + arm->getOpenCLNPars() + 1;
}

std::vector<cl_float> PlaneProfileArm::getOpenCLPars() const
{
   std::vector<cl_float> pars = fradProfile->getOpenCLPars();
   std::vector<cl_float> tmpPars = arm->getOpenCLPars();
   pars.insert(pars.end(), tmpPars.begin(), tmpPars.end());
   pars.push_back(norm);
   return pars;
}

#endif



static utl::Registry1<PlaneProfile, const utl::Branch&>::Registrar<PlaneProfileSymmetricArms> registrarSymmetricArms("SymmetricArms");

PlaneProfileSymmetricArms::PlaneProfileSymmetricArms ( const utl::Branch &b ) :
   PlaneProfile(b),
   norm(0)
{
   assert(b.GetAttributes()["type"] == "SymmetricArms");

   //The variables for the arm
   fvarNames.resize(5);
   fvarNames[0] = "norm";
   fvarNames[1] = "a";
   fvarNames[2] = "rMin";
   fvarNames[3] = "phiMin";
   fvarNames[4] = "width";

   fvars.add(fvarNames, profileName, b);

   fvarNames.resize(6);
   use_rMax = false;
   try {
      std::vector<std::string> varNames(1);
      varNames[0] = "rMax";
      fvars.add(varNames, profileName, b);
      fvarNames[5] = varNames[0];
      use_rMax = true;
   } catch (utl::Variables::VariableError) {}

   try {
      std::vector<std::string> varNames(1);
      varNames[0] = "phiExtent";
      fvars.add(varNames, profileName, b);
      fvarNames[5] = varNames[0];
      if (use_rMax) {
         ERROR("Both phiExtent and rMax given, only one can be used.");
         throw(std::invalid_argument("Invalid variables in XML files"));
      }
   } catch (utl::Variables::VariableError &e) {
      if (! use_rMax)
         throw(e);
   }

   //The arm function
   utl::Branch tb = b.GetChild("armfunction");
   if ( ! tb ) {
      std::ostringstream os;
      os << "armfunction element missing in PlaneProfileSymmetricArms \""<<profileName<<"\"";
      FATAL(os.str());
      throw(std::runtime_error("armfunction element is required"));
   }

   std::string armFName;
   tb.GetData(armFName);

   //Number of arms
   tb = b.GetChild("numberOfArms");
   if ( ! tb ) {
      std::ostringstream os;
      os << "numberOfArms element missing in PlaneProfileSymmetricArms \""<<profileName<<"\"";
      FATAL(os.str());
      throw(std::runtime_error("numberOfArms element is required"));
   }

   size_t nArms(0);
   tb.GetData(nArms);

   if (nArms < 2) {
      std::ostringstream os;
      os << "Need at least two arms in PlaneProfileSymmetricArms \""<<profileName<<"\"";
      FATAL(os.str());
      throw(std::runtime_error("Less than 2 arms specified"));
   }

   //The parameters of the arm functions are updated later
   arms.reserve(nArms);
   for (size_t i(0); i < nArms; ++i)
      arms.push_back(std::unique_ptr<ArmFunction>(new ArmFunction(1.0, 1.0, 0.0, 0.0, 1.0, utl::Registry0<ScaleFunction>::create(armFName))));

   //Loop over the children to find the radial profiles
   for ( utl::Branch rp = b.GetFirstChild();  rp; rp = rp.GetNextSibling() ) {
      if ( rp.GetBranchNameString() == "RadialProfile" ) {
         fradProfile = RadialProfile::createProfile( rp );
         break;
      }
   }

   fvars.add(fradProfile->getVariables());

   updateVariableValues( fvars );
   
}

PlaneProfileSymmetricArms::PlaneProfileSymmetricArms ( const utl::Parameters& pars, 
      std::vector<std::unique_ptr<ArmFunction> > &&armFunctions,
      std::unique_ptr<RadialProfile> radProfile, 
      const std::string &prefix, 
      bool rMax) :
   norm(0),
   arms(std::move(armFunctions)),
   fradProfile(std::move(radProfile)),
   use_rMax(rMax)
{

   if (arms.size() < 2) {
      std::cerr<<"should have at least two arms for symmetric arm profile."<<std::endl;
      throw(std::runtime_error("too few arms in model."));
   }

   //Create the variables vector
   fvarNames.resize(6);
   fvarNames[0] = prefix + "_norm";
   fvarNames[1] = prefix + "_a";
   fvarNames[2] = prefix + "_rMin";
   fvarNames[3] = prefix + "_phiMin";
   fvarNames[4] = prefix + "_width";
   if (use_rMax)
      fvarNames[5] = prefix + "_rMax";
   else
      fvarNames[5] = prefix + "_phiExtent";

   fvars.add(fvarNames, pars);

   fvars.add(fradProfile->getVariables());

   updateVariableValues( fvars );
}

PlaneProfileSymmetricArms::PlaneProfileSymmetricArms ( const utl::Parameters& pars, 
      std::vector<std::unique_ptr<ArmFunction> > &&armFunctions,
      std::unique_ptr<RadialProfile> radProfile, 
      const std::vector<std::string> &varNames,
      bool rMax):
   norm(0),
   arms(std::move(armFunctions)),
   fradProfile(std::move(radProfile)),
   fvarNames(varNames),
   use_rMax(rMax)
{

   if (arms.size() < 2) {
      std::cerr<<"Should have at least two arms for symmetric arm profile."<<std::endl;
      throw(std::runtime_error("Too few arms in model."));
   }

   if (fvarNames.size() != 6) {
      std::cerr<<"Should have exactly 6 variables for arm profile."<<std::endl;
      throw(std::runtime_error("Number of variables incorrect in model."));
   }

   fvars.add(fvarNames, pars);

   fvars.add(fradProfile->getVariables());

   updateVariableValues( fvars );

}


void PlaneProfileSymmetricArms::updateVariableValues( const utl::Variables &vars ) {

   fradProfile->updateVariableValues(vars);

   fvars.setFromVariables(vars);

   norm = vars[fvarNames[0]];
   const double a = vars[fvarNames[1]];
   const double rMin = vars[fvarNames[2]];
   const double phiMin = vars[fvarNames[3]];
   const double width = vars[fvarNames[4]];
   const double var5 = vars[fvarNames[5]];

   const double dphiMin = 2*M_PI/arms.size();

   for (size_t i=0; i < arms.size(); ++i) {
      arms[i]->setA(a);
      arms[i]->setRMin(rMin);
      arms[i]->setPhiMin(phiMin+i*dphiMin);
      arms[i]->setWidth(width);
      if (use_rMax)
         arms[i]->setPhiExtentFromRMax(var5);
      else
         arms[i]->setPhiExtent(var5);
   }

}

double PlaneProfileSymmetricArms::operator() ( double radius, double theta ) const {
   double output(0);
   for (size_t i=0; i < arms.size(); ++i)
      output += (*arms[i])(radius,theta);
   output *= norm * (*fradProfile)(radius);
   return output;
}


void PlaneProfileSymmetricArms::addToDOM( xercesc::DOMNode *node, 
      const std::map<std::string,std::string> &attributes ) const 
{
   xercesc::DOMDocument *doc = node->getOwnerDocument();

   xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("PlaneProfile").unicodeForm());
   node->appendChild(profileEl);
   profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("SymmetricArms").unicodeForm());

   for ( auto it = attributes.begin(); it != attributes.end(); ++it )
      profileEl->setAttribute(utl::XStr(it->first).unicodeForm(), utl::XStr(it->second).unicodeForm());

   //Add the name
   {
      xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
      profileEl->appendChild(element);

      xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
      element->appendChild(text);
   }

   xercesc::DOMElement* element1 = doc->createElement(utl::XStr("numberOfArms").unicodeForm());
   profileEl->appendChild(element1);

   std::ostringstream oss;
   oss << arms.size();
   xercesc::DOMText* text1 = doc->createTextNode(utl::XStr(oss.str()).unicodeForm());
   element1->appendChild(text1);

   xercesc::DOMElement* element = doc->createElement(utl::XStr("armfunction").unicodeForm());
   profileEl->appendChild(element);

   xercesc::DOMText* text = doc->createTextNode(utl::XStr(arms[0]->scaleName()).unicodeForm());
   element->appendChild(text);

   std::vector<std::string> varNames(6), varIds(6);

   varIds[0] = "norm";
   varIds[1] = "a";
   varIds[2] = "rMin";
   varIds[3] = "phiMin";
   varIds[4] = "width";
   if (use_rMax)
      varIds[5] = "rMax";
   else
      varIds[5] = "phiExtent";

   fvars.addToDOM( profileEl, fvarNames, varIds, profileName );

   std::map<std::string,std::string> attributeMap;
   fradProfile->addToDOM( profileEl, attributeMap );
}


#ifdef HAVE_OPENCL
std::string PlaneProfileSymmetricArms::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   //Update the name
   name = "PlaneProfileSymmetricArms_"+profileName+"_"+name;

   //The name for the radial profile
   std::string fradProfileName = name;

   //Use an ostringstream object to store the function name, start with the external function
   //Also keep track of parIndex, store our parameters at the back
   std::ostringstream codeStr;

   codeStr<<fradProfile->getOpenCLFunction(fradProfileName, parIndex);
   parIndex += fradProfile->getOpenCLNPars();

   std::vector<std::string> armProfileNames(arms.size());
   for (size_t i(0); i < arms.size(); ++i) {
      std::ostringstream armname;
      armname << name << i;
      armProfileNames[i] = armname.str();
      codeStr<<arms[i]->getOpenCLFunction(armProfileNames[i], parIndex);
      parIndex += arms[i]->getOpenCLNPars();
   }

   codeStr<<"float "<<name<<"( float R, float theta, __constant float *pars )\n"
      <<"{\n"
      <<"  return pars["<<parIndex<<"] * "<<fradProfileName<<"(R,pars) * (";
   for (size_t i(0); i < arms.size(); ++i) {
      codeStr<<armProfileNames[i]<<"(R,theta,pars)";
      if (i != arms.size() - 1)
         codeStr<<" + ";
   }
   codeStr<<");\n"
      <<"}\n";

   return codeStr.str();
}

size_t PlaneProfileSymmetricArms::getOpenCLNPars() const
{
   size_t npars = fradProfile->getOpenCLNPars() + 1;
   for (size_t i(0); i < arms.size(); ++i)
      npars += arms[i]->getOpenCLNPars();
   return npars;
}

std::vector<cl_float> PlaneProfileSymmetricArms::getOpenCLPars() const
{
   std::vector<cl_float> pars = fradProfile->getOpenCLPars();
   for (size_t i=0; i < arms.size(); ++i) {
      std::vector<cl_float> tmpPars = arms[i]->getOpenCLPars();
      pars.insert(pars.end(), tmpPars.begin(), tmpPars.end());
   }
   pars.push_back(norm);
   return pars;
}

#endif


}
