#include "cylindricalprofiles.h"
#include "Reader.h"
#include "spiralarms.h"
#include "Registry.h"
#include <ReaderErrorReporter.h>
#include <PhysicalConstants.h>

namespace GalacticStructure {

   CylindricalProfile::CylindricalProfile( const utl::Branch &b )
   {

      assert(b.GetBranchNameString() == "CylindricalProfile");

      utl::Branch nb = b.GetChild("name");
      if ( ! nb) {
         std::ostringstream os;
         os<<"A name is required for cylindrical profiles. Please specify the name element.";
         FATAL(os.str());
         throw(std::runtime_error("No name element found"));
      }
      nb.GetData(profileName);

   }


   std::unique_ptr<CylindricalProfile> CylindricalProfile::createProfile( const utl::Branch &b ) {
      assert(b.GetBranchNameString() == "CylindricalProfile");

      const std::string name = b.GetAttributes()["type"];

      return utl::Registry1<CylindricalProfile, const utl::Branch&>::create(name,b);
   }

#ifdef HAVE_OPENCL
   std::string CylindricalProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
   {
      ERROR("OpenCL implementation not ready for this class");
      throw(std::runtime_error("Not Implemented"));
   }

   size_t CylindricalProfile::getOpenCLNPars() const
   {
      ERROR("OpenCL implementation not ready for this class");
      throw(std::runtime_error("Not Implemented"));
   }

   std::vector<cl_float> CylindricalProfile::getOpenCLPars() const
   {
      ERROR("OpenCL implementation not ready for this class");
      throw(std::runtime_error("Not Implemented"));
   }
#endif




   static utl::Registry1<CylindricalProfile, const utl::Branch&>::Registrar<GenericDiskProfile> registrarGeneric("Generic");

   GenericDiskProfile::GenericDiskProfile( const utl::Parameters &pars,
         std::unique_ptr<PlaneProfile> planeDensity,
         std::unique_ptr<PlaneProfile> planeCenter,
         std::unique_ptr<PlaneProfile> northScaleHeight,
         std::unique_ptr<PlaneProfile> southScaleHeight,
         std::unique_ptr<ScaleFunction> northScaleFunction,
         std::unique_ptr<ScaleFunction> southScaleFunction ) :
      density(std::move(planeDensity)),
      center(std::move(planeCenter)),
      northHeight(std::move(northScaleHeight)),
      southHeight(std::move(southScaleHeight)),
      northScale(std::move(northScaleFunction)),
      southScale(std::move(southScaleFunction))
   {

      fvars.add(density->getVariables());
      fvars.add(center->getVariables());
      fvars.add(northHeight->getVariables());
      fvars.add(southHeight->getVariables());

   }


   GenericDiskProfile::GenericDiskProfile( const utl::Branch &b ) :
      CylindricalProfile(b)
   {

      assert(b.GetAttributes()["type"] == "Generic");

      //Loop over the children to find the plane profiles
      for ( utl::Branch rp = b.GetFirstChild();  rp; rp = rp.GetNextSibling() ) {

         if ( rp.GetBranchNameString() == "PlaneProfile" ) {

            //Look for the intent attribute
            std::string upper = rp.GetAttributes()["intent"];
            // explicit cast needed to resolve ambiguity
            std::transform(upper.begin(), upper.end(), upper.begin(), (int(*)(int)) std::toupper);

            if ( upper == "PLANEDENSITY" ) {
               density = PlaneProfile::createProfile( rp );
            } else if ( upper == "PLANECENTER" ) {
               center = PlaneProfile::createProfile( rp );
            } else if ( upper == "NORTHSCALEHEIGHT" ) {
               northHeight = PlaneProfile::createProfile( rp );
            } else if ( upper == "SOUTHSCALEHEIGHT" ) {
               southHeight = PlaneProfile::createProfile( rp );
            }
         }
      }

      //And the scale functions
      std::string scaleFunc;
      utl::Branch tb = b.GetChild("northScaleFunction");
      tb.GetData(scaleFunc);
      northScale = utl::Registry0<ScaleFunction>::create(scaleFunc);

      tb = b.GetChild("southScaleFunction");
      tb.GetData(scaleFunc);
      southScale = utl::Registry0<ScaleFunction>::create(scaleFunc);

      //Register all the variables
      fvars.add(density->getVariables());
      fvars.add(center->getVariables());
      fvars.add(northHeight->getVariables());
      fvars.add(southHeight->getVariables());
   }


   void GenericDiskProfile::updateVariableValues( const utl::Variables & vars ) {

      density->updateVariableValues(vars);
      center->updateVariableValues(vars);
      northHeight->updateVariableValues(vars);
      southHeight->updateVariableValues(vars);
   
   }


   double GenericDiskProfile::operator() ( double r, double theta, double z ) const {

      const double zm = z - (*center)(r,theta);

      const double pd = (*density)(r,theta);
      
      if (zm > 0)
         return pd * (*northScale)(zm*(*northHeight)(r,theta));
      else if (zm < 0)
         return pd * (*southScale)(-zm*(*southHeight)(r,theta));
      else
         return pd;

   }

   void GenericDiskProfile::addToDOM( xercesc::DOMNode *node ) const 
   {
      xercesc::DOMDocument *doc = node->getOwnerDocument();

      xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("CylindricalProfile").unicodeForm());
      node->appendChild(profileEl);
      profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("Generic").unicodeForm());

      //Add the name
      {
         xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
         profileEl->appendChild(element);

         xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
         element->appendChild(text);
      }

      std::map<std::string, std::string> attributeMap;
      attributeMap["intent"] = "PlaneDensity";
      density->addToDOM(profileEl, attributeMap);
      attributeMap["intent"] = "PlaneCenter";
      center->addToDOM(profileEl, attributeMap);
      attributeMap["intent"] = "NorthScaleHeight";
      northHeight->addToDOM(profileEl, attributeMap);
      attributeMap["intent"] = "SouthScaleHeight";
      southHeight->addToDOM(profileEl, attributeMap);

      xercesc::DOMElement* elementnorth = doc->createElement(utl::XStr("northScaleFunction").unicodeForm());
      profileEl->appendChild(elementnorth);
      xercesc::DOMText* textnorth = doc->createTextNode(utl::XStr(northScale->name()).unicodeForm());
      elementnorth->appendChild(textnorth);

      xercesc::DOMElement* elementsouth = doc->createElement(utl::XStr("southScaleFunction").unicodeForm());
      profileEl->appendChild(elementsouth);
      xercesc::DOMText* textsouth = doc->createTextNode(utl::XStr(southScale->name()).unicodeForm());
      elementsouth->appendChild(textsouth);

   }

#ifdef HAVE_OPENCL
   std::string GenericDiskProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
   {
      //Update the name, use this name for the other functions
      name = "GenericDiskProfile_"+profileName+"_"+name;

      //The function names
      std::string densityName = name+"density";
      std::string centerName = name+"center";
      std::string northHeightName = name+"northHeight";
      std::string southHeightName = name+"southHeight";
      std::string northScaleName = name+"northScale";
      std::string southScaleName = name+"southScale";

      //Use an ostringstream object to store the function name, start with the external functions
      //Keep track of the parIndex
      std::ostringstream codeStr;

      codeStr<<density->getOpenCLFunction(densityName, parIndex);
      parIndex += density->getOpenCLNPars();
      
      codeStr<<center->getOpenCLFunction(centerName, parIndex);
      parIndex += center->getOpenCLNPars();
      
      codeStr<<northHeight->getOpenCLFunction(northHeightName, parIndex);
      parIndex += northHeight->getOpenCLNPars();
      
      codeStr<<southHeight->getOpenCLFunction(southHeightName, parIndex);
      parIndex += southHeight->getOpenCLNPars();
      
      codeStr<<northScale->getOpenCLFunction(northScaleName);
      
      codeStr<<southScale->getOpenCLFunction(southScaleName);

      //Now it is the code for the disk profile
      codeStr<<"inline float "<<name<<"( float R, float theta, float z, __constant float *pars )\n"
         <<"{\n"
         <<"const float zm = z - "<<centerName<<"(R,theta,pars);\n"
         <<"if (zm > 0)\n"
         <<"return "<<densityName<<"(R,theta,pars) * "<<northScaleName<<"(zm * "<<northHeightName<<"(R,theta,pars) );\n"
         <<"else\n"
         <<"return "<<densityName<<"(R,theta,pars) * "<<southScaleName<<"(-zm * "<<southHeightName<<"(R,theta,pars) );\n"
         <<"}\n";

      return codeStr.str();

   }

   size_t GenericDiskProfile::getOpenCLNPars() const
   {
      size_t nPars(0);
      nPars += density->getOpenCLNPars();
      nPars += center->getOpenCLNPars();
      nPars += northHeight->getOpenCLNPars();
      nPars += southHeight->getOpenCLNPars();
      return nPars;
   }

   std::vector<cl_float> GenericDiskProfile::getOpenCLPars() const
   {
      std::vector<cl_float> pars, tmpPars;

      tmpPars = density->getOpenCLPars();
      pars.insert(pars.end(), tmpPars.begin(), tmpPars.end());

      tmpPars = center->getOpenCLPars();
      pars.insert(pars.end(), tmpPars.begin(), tmpPars.end());

      tmpPars = northHeight->getOpenCLPars();
      pars.insert(pars.end(), tmpPars.begin(), tmpPars.end());

      tmpPars = southHeight->getOpenCLPars();
      pars.insert(pars.end(), tmpPars.begin(), tmpPars.end());

      return pars;
   }

#endif





   static utl::Registry1<CylindricalProfile, const utl::Branch&>::Registrar<NE2001ArmProfile> registrarNE2001Arm("NE2001Arm");

   NE2001ArmProfile::NE2001ArmProfile( const utl::Branch &b ) :
      CylindricalProfile(b),
      arm1(4.25, 3.48, 0.0,   15.0, 0.6, utl::Registry0<ScaleFunction>::create("Gaussian")),
      arm2(4.25, 3.48, 3.141, 15.0, 0.6*1.5, utl::Registry0<ScaleFunction>::create("Gaussian")),
      arm3(4.89, 4.90, 2.525, 15.0, 0.6, utl::Registry0<ScaleFunction>::create("Gaussian")),
      arm4(4.89, 3.76, 4.24,  15.0, 0.6*0.8, utl::Registry0<ScaleFunction>::create("Gaussian")),
      arm5(4.57, 8.10, 5.847, 0.55, 0.6, utl::Registry0<ScaleFunction>::create("Gaussian"))
   {
      assert(b.GetAttributes()["type"] == "NE2001Arm");
      init();
   }

   NE2001ArmProfile::NE2001ArmProfile( ) :
      arm1(4.25, 3.48, 0.0,   15.0, 0.6, utl::Registry0<ScaleFunction>::create("Gaussian")),
      arm2(4.25, 3.48, 3.141, 15.0, 0.6*1.5, utl::Registry0<ScaleFunction>::create("Gaussian")),
      arm3(4.89, 4.90, 2.525, 15.0, 0.6, utl::Registry0<ScaleFunction>::create("Gaussian")),
      arm4(4.89, 3.76, 4.24,  15.0, 0.6*0.8, utl::Registry0<ScaleFunction>::create("Gaussian")),
      arm5(4.57, 8.10, 5.847, 0.55, 0.6, utl::Registry0<ScaleFunction>::create("Gaussian"))

   {
      init();
   }

   void NE2001ArmProfile::init() 
   {

      na = 0.030;
      wa = 0.6; 
      ha = 0.25;
      Aa = 11.0;

      //Set the arm variables
      arm1.setA(4.25);
      arm2.setA(4.25);
      arm3.setA(4.89);
      arm4.setA(4.89);
      arm5.setA(4.57);

      arm1.setRMin(3.48);
      arm2.setRMin(3.48);
      arm3.setRMin(4.90);
      arm4.setRMin(3.76);
      arm5.setRMin(8.10);

      arm1.setPhiMin(0.0);
      arm2.setPhiMin(3.141);
      arm3.setPhiMin(2.525);
      arm4.setPhiMin(4.24);
      arm5.setPhiMin(5.847);

      //Cut all arms (except 5) at 2*Aa to speed things up
      arm1.setPhiExtentFromRMax(2*Aa);
      arm2.setPhiExtentFromRMax(2*Aa);
      arm3.setPhiExtentFromRMax(2*Aa);
      arm4.setPhiExtentFromRMax(2*Aa);
      arm5.setPhiExtent(0.55);

      //Shape the arms according to NE2001 model
      //Start with 2nd arm which is arm 3 in TC model
      std::vector< std::pair< double, double > > armRP = arm2.RPArmData();

      const double ll1 = 370*M_PI/180.;
      const double ul1 = 410*M_PI/180.;
      const double c1 = (ul1+ll1)/2.;
      const double scale1 = M_PI/(ul1-ll1);
      const double a1 = 0.04;

      const double ll2 = 315*M_PI/180.;
      const double ul2 = 370*M_PI/180.;
      const double c2 = (ul2+ll2)/2.;
      const double scale2 = M_PI/(ul2-ll2);
      const double a2 = -0.07;

      const double ll3 = 180*M_PI/180.;
      const double ul3 = 315*M_PI/180.;
      const double c3 = (ul3+ll3)/2.;
      const double scale3 = M_PI/(ul3-ll3);
      const double a3 = 0.16;

      for ( size_t i(0); i < armRP.size(); ++i ) {

         // 370 deg < phi < 410 deg
         if ( armRP[i].second > ll1 && armRP[i].second < ul1 )
            armRP[i].first *= 1. + a1*cos((armRP[i].second-c1)*scale1);

         // 315 deg < phi < 370 deg 
         if ( armRP[i].second > ll2 && armRP[i].second < ul2 )
            armRP[i].first *= 1. + a2*cos((armRP[i].second-c2)*scale2);

         // 180 deg < phi < 315 deg
         if ( armRP[i].second > ll3 && armRP[i].second < ul3 )
            armRP[i].first *= 1. + a3*cos((armRP[i].second-c3)*scale3);

      }

      arm2.setRPArmData(armRP);

      //4th arm which is 2 in TC model
      armRP = arm4.RPArmData();

      const double ll4 = 290*M_PI/180.;
      const double ul4 = 395*M_PI/180.;
      const double c4 = (ul4+ll4)/2.;
      const double scale4 = M_PI/(ul4-ll4);
      const double a4 = -0.11;

      for ( size_t i(0); i < armRP.size(); ++i ) {

         // 290 < phi < 395
         if ( armRP[i].second > ll4 && armRP[i].second < ul4 ) 
            armRP[i].first *= 1. + a4*cos((armRP[i].second-c4)*scale4);

      }

      arm4.setRPArmData(armRP);

      updateVariableValues(fvars);
   }

   void NE2001ArmProfile::updateVariableValues( const utl::Variables &vars ) {

      wa1 = wa; 
      wa2 = wa*1.5; 
      wa3 = wa; 
      wa4 = wa*0.8; 
      wa5 = wa;

      arm1.setWidth(wa1);
      arm2.setWidth(wa2);
      arm3.setWidth(wa3);
      arm4.setWidth(wa4);
      arm5.setWidth(wa5);

      na1 = na*0.5;
      na2 = na*1.2;
      na3 = na*1.3;
      na4 = na;
      na5 = na*0.25;

      ha1 = ha; 
      ha2 = ha*0.8; 
      ha3 = ha*1.3;
      ha4 = ha*1.5;
      ha5 = ha;

   }

   double NE2001ArmProfile::operator() (double r, double theta, double z) const {

      double out = arm1(r, theta) * na1 * sech2( z / ha1 );
      out += arm2(r, theta) * na2 * sech2( z / ha2 );
      out += arm3(r, theta) * na3 * sech2( z / ha3 );
      out += arm4(r, theta) * na4 * sech2( z / ha4 );
      out += arm5(r, theta) * na5 * sech2( z / ha5 );

      if (r > Aa)
         return out*sech2((r-Aa)/2.);
      else
         return out;

   }

   void NE2001ArmProfile::addToDOM( xercesc::DOMNode *node ) const 
   {
      xercesc::DOMDocument *doc = node->getOwnerDocument();

      xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("CylindricalProfile").unicodeForm());
      node->appendChild(profileEl);
      profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("NE2001Arm").unicodeForm());

      //Add the name
      {
         xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
         profileEl->appendChild(element);

         xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
         element->appendChild(text);
      }

   }



   static utl::Registry1<CylindricalProfile, const utl::Branch&>::Registrar<WainscoatBulgeProfile> registrarWainscoatBulge("WainscoatBulge");

   WainscoatBulgeProfile::WainscoatBulgeProfile ( double Rmin_, double R1, double k1_, double norm_,
           const std::string &Rminname, const std::string &R1name, const std::string &k1name, const std::string &normname ) :
      fRminname(Rminname),
      fR1name(R1name),
      fk1name(k1name),
      fnormname(normname)
   {
      fvars.add(Rminname, Rmin_);
      fvars.add(R1name, R1);
      fvars.add(k1name, k1_);
      fvars.add(normname, norm_);

      updateVariableValues( fvars );
   }

   WainscoatBulgeProfile::WainscoatBulgeProfile ( const utl::Parameters &pars, 
           const std::string &Rminname, const std::string &R1name, const std::string &k1name, const std::string &normname ) :
      fRminname(Rminname),
      fR1name(R1name),
      fk1name(k1name),
      fnormname(normname)
   {
      std::vector<std::string> varNames(4);
      varNames[0] = Rminname;
      varNames[1] = R1name;
      varNames[2] = k1name;
      varNames[3] = normname;

      fvars.add(varNames, pars);

      updateVariableValues( fvars );
   }

   WainscoatBulgeProfile::WainscoatBulgeProfile ( const utl::Branch &b ) :
      CylindricalProfile(b)
   {
      assert(b.GetAttributes()["type"] == "WainscoatBulge");

      std::vector<std::string> varNames(4);
      varNames[0] = "Rmin";
      varNames[1] = "R1";
      varNames[2] = "k1";
      varNames[3] = "norm";

      fvars.add(varNames, profileName, b);

      fRminname = varNames[0];
      fR1name = varNames[1];
      fk1name = varNames[2];
      fnormname = varNames[3];

      updateVariableValues( fvars );
   }

   void WainscoatBulgeProfile::updateVariableValues( const utl::Variables &vars )
   {
      fvars.setFromVariables(vars);

      Rmin = fvars[fRminname];
      oneOver_R1 = 1./fvars[fR1name];
      k1 = fvars[fk1name];
      norm = fvars[fnormname];
   }

   double WainscoatBulgeProfile::operator() (double r, double theta, double z ) const
   {
      const double rs = (r < Rmin ? Rmin : r), zs = k1*z,
            x = sqrt(rs*rs + zs*zs)*oneOver_R1;

      return norm * exp(-x*x*x)/pow(x, 1.8);
   }

   void WainscoatBulgeProfile::addToDOM( xercesc::DOMNode *node ) const 
   {
      xercesc::DOMDocument *doc = node->getOwnerDocument();

      xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("CylindricalProfile").unicodeForm());
      node->appendChild(profileEl);
      profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("WainscoatBulge").unicodeForm());

      //Add the name
      {
         xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
         profileEl->appendChild(element);

         xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
         element->appendChild(text);
      }

      std::vector<std::string> varIds(4), varNames(4);

      varIds[0] = "Rmin";
      varIds[1] = "R1";
      varIds[2] = "k1";
      varIds[3] = "norm";

      varNames[0] = fRminname;
      varNames[1] = fR1name;
      varNames[2] = fk1name;
      varNames[3] = fnormname;

      fvars.addToDOM( profileEl, varNames, varIds, profileName );
   }

   void WainscoatBulgeProfile::setRmin(double r) {
      Rmin = r;
   }

   void WainscoatBulgeProfile::setR1(double r) {
      oneOver_R1 = 1./r;
   }

   void WainscoatBulgeProfile::setk1(double k) {
      k1 = k;
   }
      
   void WainscoatBulgeProfile::setnorm(double k) {
      norm = k;
   }
      



   static utl::Registry1<CylindricalProfile, const utl::Branch&>::Registrar<WainscoatHaloProfile> registrarWainscoatHalo("WainscoatHalo");

   WainscoatHaloProfile::WainscoatHaloProfile ( double rE, double zE_, double norm_, 
         const std::string &rEname, const std::string &zEname, const std::string &normname ) :
      frEname(rEname),
      fzEname(zEname),
      fnormname(normname)
   {
      fvars.add(rEname, rE);
      fvars.add(zEname, zE_);
      fvars.add(normname, norm_);

      updateVariableValues( fvars );
   }

   WainscoatHaloProfile::WainscoatHaloProfile ( const utl::Parameters &pars, 
         const std::string &rEname, const std::string &zEname, const std::string &normname ) :
      frEname(rEname),
      fzEname(zEname),
      fnormname(normname)
   {
      std::vector<std::string> varNames(3);
      varNames[0] = rEname;
      varNames[1] = zEname;
      varNames[2] = normname;

      fvars.add(varNames, pars);

      updateVariableValues( fvars );
   }

   WainscoatHaloProfile::WainscoatHaloProfile ( const utl::Branch &b ) :
      CylindricalProfile(b)
   {
      assert(b.GetAttributes()["type"] == "WainscoatHalo");

      std::vector<std::string> varNames(3);
      varNames[0] = "rE";
      varNames[1] = "zE";
      varNames[2] = "norm";

      fvars.add(varNames, profileName, b);

      frEname = varNames[0];
      fzEname = varNames[1];
      fnormname = varNames[2];

      updateVariableValues( fvars );
   }

   void WainscoatHaloProfile::updateVariableValues( const utl::Variables &vars )
   {
      fvars.setFromVariables(vars);

      oneOver_rE = 1./fvars[frEname];
      zE = fvars[fzEname];
      norm = fvars[fnormname];
   }

   double WainscoatHaloProfile::operator() (double r, double theta, double z ) const
   {
      const double zp = zE*z;

      const double alpha = sqrt(r*r + zp*zp)*oneOver_rE;

      return norm * pow(10., -3.3307*(pow(alpha, 0.25) - 1.));
   }


   void WainscoatHaloProfile::addToDOM( xercesc::DOMNode *node ) const 
   {
      xercesc::DOMDocument *doc = node->getOwnerDocument();

      xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("CylindricalProfile").unicodeForm());
      node->appendChild(profileEl);
      profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("WainscoatHalo").unicodeForm());

      //Add the name
      {
         xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
         profileEl->appendChild(element);

         xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
         element->appendChild(text);
      }

      std::vector<std::string> varIds(3), varNames(3);

      varIds[0] = "rE";
      varIds[1] = "zE";
      varIds[2] = "norm";

      varNames[0] = frEname;
      varNames[1] = fzEname;
      varNames[2] = fnormname;

      fvars.addToDOM( profileEl, varNames, varIds, profileName );
   }

   void WainscoatHaloProfile::setrE(double r) {
      oneOver_rE = 1./r;
   }

   void WainscoatHaloProfile::setzE(double z) {
      zE = z;
   }
      
   void WainscoatHaloProfile::setnorm(double z) {
      norm = z;
   }
      



   static utl::Registry1<CylindricalProfile, const utl::Branch&>::Registrar<SDSSHaloProfile> registrarSDSSHalo("SDSSHalo");

   SDSSHaloProfile::SDSSHaloProfile ( double rS_, double qH, double alpha_, double norm_,
         const std::string &rSname, const std::string &qHname, const std::string &alphaname, const std::string &normname ) :
      frSname(rSname),
      fqHname(qHname),
      falphaname(alphaname),
      fnormname(normname)
   {
      fvars.add(rSname, rS_);
      fvars.add(qHname, qH);
      fvars.add(alphaname, alpha);
      fvars.add(normname, norm);

      updateVariableValues( fvars );
   }

   SDSSHaloProfile::SDSSHaloProfile ( const utl::Parameters &pars, 
         const std::string &rSname, const std::string &qHname, const std::string &alphaname, const std::string &normname ) :
      frSname(rSname),
      fqHname(qHname),
      falphaname(alphaname),
      fnormname(normname)
   {
      std::vector<std::string> varNames(4);
      varNames[0] = rSname;
      varNames[1] = qHname;
      varNames[2] = alphaname;
      varNames[3] = normname;

      fvars.add(varNames, pars);

      updateVariableValues( fvars );
   }

   SDSSHaloProfile::SDSSHaloProfile ( const utl::Branch &b ) :
      CylindricalProfile(b)
   {
      assert(b.GetAttributes()["type"] == "SDSSHalo");

      std::vector<std::string> varNames(4);
      varNames[0] = "rS";
      varNames[1] = "qH";
      varNames[2] = "alpha";
      varNames[3] = "norm";

      fvars.add(varNames, profileName, b);

      frSname = varNames[0];
      fqHname = varNames[1];
      falphaname = varNames[2];
      fnormname = varNames[3];

      updateVariableValues( fvars );
   }

   void SDSSHaloProfile::updateVariableValues( const utl::Variables &vars )
   {
      fvars.setFromVariables(vars);

      rS = fvars[frSname];
      oneOver_qH = 1./fvars[fqHname];
      alpha = fvars[falphaname];
      norm = fvars[fnormname];
   }

   double SDSSHaloProfile::operator() (double r, double theta, double z ) const
   {
      const double zOnQH = z*oneOver_qH, rp = (r < 0.1 ? 0.1 : r);

      return norm * pow(rS/sqrt(rp*rp + zOnQH*zOnQH), alpha);
   }

   void SDSSHaloProfile::addToDOM( xercesc::DOMNode *node ) const 
   {
      xercesc::DOMDocument *doc = node->getOwnerDocument();

      xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("CylindricalProfile").unicodeForm());
      node->appendChild(profileEl);
      profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("SDSSHalo").unicodeForm());

      //Add the name
      {
         xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
         profileEl->appendChild(element);

         xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
         element->appendChild(text);
      }

      std::vector<std::string> varIds(4), varNames(4);

      varIds[0] = "rS";
      varIds[1] = "qH";
      varIds[2] = "alpha";
      varIds[3] = "norm";

      varNames[0] = frSname;
      varNames[1] = fqHname;
      varNames[2] = falphaname;
      varNames[3] = fnormname;

      fvars.addToDOM( profileEl, varNames, varIds, profileName );
   }

   void SDSSHaloProfile::setrS(double r) {
      rS = r;
   }

   void SDSSHaloProfile::setqH(double q) {
      oneOver_qH = 1./q;
   }

   void SDSSHaloProfile::setalpha(double a) {
      alpha = a;
   }

   void SDSSHaloProfile::setnorm(double a) {
      norm = a;
   }




   static utl::Registry1<CylindricalProfile, const utl::Branch&>::Registrar<LopezCorredoiraBulgeProfile> registrarLopezCorredoiraBulge("LopezCorredoiraBulge");

   LopezCorredoiraBulgeProfile::LopezCorredoiraBulgeProfile ( double phiOffset, double ya, double za, double index, double scaleLength, double norm,
         const std::string &phiOffsetname, const std::string &yaname, const std::string &zaname, 
         const std::string &indexname, const std::string &scaleLengthname, const std::string &normname ) :
      fphiOffsetname(phiOffsetname),
      fyaname(yaname),
      fzaname(zaname),
      findexname(indexname),
      fscaleLengthname(scaleLengthname),
      fnormname(normname)
   {
      fvars.add(phiOffsetname, phiOffset);
      fvars.add(yaname, ya);
      fvars.add(zaname, za);
      fvars.add(indexname, index);
      fvars.add(scaleLengthname, scaleLength);
      fvars.add(normname, norm);

      updateVariableValues( fvars );
   }

   LopezCorredoiraBulgeProfile::LopezCorredoiraBulgeProfile ( const utl::Parameters &pars, 
         const std::string &phiOffsetname, const std::string &yaname, const std::string &zaname, 
         const std::string &indexname, const std::string &scaleLengthname, const std::string &normname ) :
      fphiOffsetname(phiOffsetname),
      fyaname(yaname),
      fzaname(zaname),
      findexname(indexname),
      fscaleLengthname(scaleLengthname),
      fnormname(normname)
   {
      std::vector<std::string> varNames(6);
      varNames[0] = phiOffsetname;
      varNames[1] = yaname;
      varNames[2] = zaname;
      varNames[3] = indexname;
      varNames[4] = scaleLengthname;
      varNames[5] = normname;

      fvars.add(varNames, pars);

      updateVariableValues( fvars );
   }

   LopezCorredoiraBulgeProfile::LopezCorredoiraBulgeProfile ( const utl::Branch &b ) :
      CylindricalProfile(b)
   {
      assert(b.GetAttributes()["type"] == "LopezCorredoiraBulge");

      std::vector<std::string> varNames(6);
      varNames[0] = "phiOffset";
      varNames[1] = "ya";
      varNames[2] = "za";
      varNames[3] = "index";
      varNames[4] = "scaleLength";
      varNames[5] = "norm";

      fvars.add(varNames, profileName, b);

      fphiOffsetname = varNames[0];
      fyaname = varNames[1];
      fzaname = varNames[2];
      findexname = varNames[3];
      fscaleLengthname = varNames[4];
      fnormname = varNames[5];

      updateVariableValues( fvars );
   }

   void LopezCorredoiraBulgeProfile::updateVariableValues( const utl::Variables &vars )
   {
      fvars.setFromVariables(vars);

      const double phiOffset = fvars[fphiOffsetname];
      cosPhiOffset = cos(phiOffset);
      sinPhiOffset = sin(phiOffset);
      oneOver_ya = 1./fvars[fyaname];
      oneOver_za = 1./fvars[fzaname];
      index = fvars[findexname];
      oneOver_index = 1./index;
      oneOver_scaleLength = 1./fvars[fscaleLengthname];
      norm = fvars[fnormname];
   }

   double LopezCorredoiraBulgeProfile::operator() (double r, double theta, double z ) const
   {
      const double x = r*cos(theta), y = r*sin(theta),
            xp = cosPhiOffset*x + sinPhiOffset*y,
            yp = -sinPhiOffset*x + cosPhiOffset*y;

      const double xpIndex = pow(fabs(xp), index), ypIndex = pow(fabs(yp*oneOver_ya), index),
            zpIndex = pow(fabs(z*oneOver_za), index);

      const double t = pow(xpIndex + ypIndex + zpIndex, oneOver_index);

      return norm*exp(-t*oneOver_scaleLength);

   }

   void LopezCorredoiraBulgeProfile::addToDOM( xercesc::DOMNode *node ) const 
   {
      xercesc::DOMDocument *doc = node->getOwnerDocument();

      xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("CylindricalProfile").unicodeForm());
      node->appendChild(profileEl);
      profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("LopezCorredoiraBulge").unicodeForm());

      //Add the name
      {
         xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
         profileEl->appendChild(element);

         xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
         element->appendChild(text);
      }

      std::vector<std::string> varIds(6), varNames(6);

      varIds[0] = "phiOffset";
      varIds[1] = "ya";
      varIds[2] = "za";
      varIds[3] = "index";
      varIds[4] = "scaleLength";
      varIds[5] = "norm";

      varNames[0] = fphiOffsetname;
      varNames[1] = fyaname;
      varNames[2] = fzaname;
      varNames[3] = findexname;
      varNames[4] = fscaleLengthname;
      varNames[5] = fnormname;

      fvars.addToDOM( profileEl, varNames, varIds, profileName );
   }

   void LopezCorredoiraBulgeProfile::setphiOffset(double phi) {
      cosPhiOffset = cos(phi);
      sinPhiOffset = sin(phi);
   }

   void LopezCorredoiraBulgeProfile::setya(double q) {
      oneOver_ya = 1./q;
   }

   void LopezCorredoiraBulgeProfile::setza(double q) {
      oneOver_za = 1./q;
   }

   void LopezCorredoiraBulgeProfile::setindex(double i) {
      index = i;
      oneOver_index = 1./i;
   }

   void LopezCorredoiraBulgeProfile::setscaleLength(double q) {
      oneOver_scaleLength = 1./q;
   }

   void LopezCorredoiraBulgeProfile::setnorm(double a) {
      norm = a;
   }




   static utl::Registry1<CylindricalProfile, const utl::Branch&>::Registrar<FreudenreichBarProfile> registrarFreudenreichBar("FreudenreichBar");

  FreudenreichBarProfile::FreudenreichBarProfile ( double phiOffset, double pitchAngle, double barX, double barY, double barZ, double barPerp, double barPara, double barREnd, double barHEnd, double norm, const std::string &phiOffsetname, const std::string& pitchAnglename, const std::string &barXname, const std::string &barYname, const std::string &barZname, const std::string &barPerpname, const std::string &barParaname, const std::string &barREndname, const std::string &barHEndname, const std::string &normname ) :
    fphiOffsetname(phiOffsetname),
    fpitchAnglename(pitchAnglename),
    fbarXname(barXname),
    fbarYname(barYname),
    fbarZname(barZname),
    fbarPerpname(barPerpname),
    fbarParaname(barParaname),
    fbarREndname(barREndname),
    fbarHEndname(barHEndname),
    fnormname(normname)
  {
    fvars.add(phiOffsetname, phiOffset);
    fvars.add(pitchAnglename, pitchAngle);
    fvars.add(barXname, barX);
    fvars.add(barYname, barY);
    fvars.add(barZname, barZ);
    fvars.add(barPerpname, barPerp);
    fvars.add(barParaname, barPara);
    fvars.add(barREndname, barREnd);
    fvars.add(barHEndname, barHEnd);
    fvars.add(normname, norm);
    
    updateVariableValues( fvars );
  }
  
  FreudenreichBarProfile::FreudenreichBarProfile ( const utl::Parameters &pars, 
						   const std::string &phiOffsetname, const std::string& pitchAnglename, const std::string &barXname, const std::string &barYname, const std::string &barZname, const std::string &barPerpname, const std::string &barParaname, const std::string &barREndname, const std::string &barHEndname, const std::string &normname ) :
    fphiOffsetname(phiOffsetname),
    fpitchAnglename(pitchAnglename),
    fbarXname(barXname),
    fbarYname(barYname),
    fbarZname(barZname),
    fbarPerpname(barPerpname),
    fbarParaname(barParaname),
    fbarREndname(barREndname),
    fbarHEndname(barHEndname),
    fnormname(normname)
   {
      std::vector<std::string> varNames(10);

      varNames[0] = phiOffsetname;
      varNames[1] = pitchAnglename;
      varNames[2] = barXname;
      varNames[3] = barYname;
      varNames[4] = barZname;
      varNames[5] = barPerpname;
      varNames[6] = barParaname;
      varNames[7] = barREndname;
      varNames[8] = barHEndname;
      varNames[9] = normname;

      fvars.add(varNames, pars);

      updateVariableValues( fvars );
   }

   FreudenreichBarProfile::FreudenreichBarProfile ( const utl::Branch &b ) :
      CylindricalProfile(b)
   {
      assert(b.GetAttributes()["type"] == "FreudenreichBar");

      std::vector<std::string> varNames(10);
      varNames[0] = "phiOffset";
      varNames[1] = "pitchAngle";
      varNames[2] = "barX";
      varNames[3] = "barY";
      varNames[4] = "barZ";
      varNames[5] = "barPerp";
      varNames[6] = "barPara";
      varNames[7] = "barREnd";
      varNames[8] = "barHEnd";
      varNames[9] = "norm";

      fvars.add(varNames, profileName, b);

      fphiOffsetname = varNames[0];
      fpitchAnglename = varNames[1];
      fbarXname = varNames[2];
      fbarYname = varNames[3];
      fbarZname = varNames[4];
      fbarPerpname = varNames[5];
      fbarParaname = varNames[6];
      fbarREndname = varNames[7];
      fbarHEndname = varNames[8];
      fnormname = varNames[9];

      updateVariableValues( fvars );
   }

   void FreudenreichBarProfile::updateVariableValues( const utl::Variables &vars )
   {
      fvars.setFromVariables(vars);

      setphiOffset(fvars[fphiOffsetname]);
      setpitchAngle(fvars[fpitchAnglename]);
      setbarX(fvars[fbarXname]);
      setbarY(fvars[fbarYname]);
      setbarZ(fvars[fbarZname]);
      setbarPerp(fvars[fbarPerpname]);
      setbarPara(fvars[fbarParaname]);
      setbarREnd(fvars[fbarREndname]);
      setbarHEnd(fvars[fbarHEndname]);
      setnorm(fvars[fnormname]);

   }

   void FreudenreichBarProfile::addToDOM( xercesc::DOMNode *node ) const 
   {
      xercesc::DOMDocument *doc = node->getOwnerDocument();

      xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("CylindricalProfile").unicodeForm());
      node->appendChild(profileEl);
      profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("FreudenreichBar").unicodeForm());

      //Add the name
      {
         xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
         profileEl->appendChild(element);

         xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
         element->appendChild(text);
      }

      std::vector<std::string> varIds(10), varNames(10);

      varIds[0] = "phiOffset";
      varIds[1] = "pitchAngle";
      varIds[2] = "barX";
      varIds[3] = "barY";
      varIds[4] = "barZ";
      varIds[5] = "barPerp";
      varIds[6] = "barPara";
      varIds[7] = "barREnd";
      varIds[8] = "barHEnd";
      varIds[9] = "norm";

      varNames[0] = fphiOffsetname;
      varNames[1] = fpitchAnglename;
      varNames[2] = fbarXname;
      varNames[3] = fbarYname;
      varNames[4] = fbarZname;
      varNames[5] = fbarPerpname;
      varNames[6] = fbarParaname;
      varNames[7] = fbarREndname;
      varNames[8] = fbarHEndname;
      varNames[9] = fnormname;

      fvars.addToDOM( profileEl, varNames, varIds, profileName );
   }

#ifdef HAVE_OPENCL
   std::string FreudenreichBarProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
   {
      //Update the name, use this name for the other functions
      name = "FreudenreichBar_"+profileName+"_"+name;

      //Use an ostringstream object to store the function 
      std::ostringstream codeStr;
      
      //To reduce the number of registers, put calculation of rs into one big line.
      std::ostringstream xp,yp,zp,rperp;
      xp <<"( pars["<<parIndex+0<<"]*pars["<<parIndex+2<<"]*R*cos(theta) "
         <<"+ pars["<<parIndex+1<<"]*R*sin(theta) "
         <<"- pars["<<parIndex+0<<"]*pars["<<parIndex+3<<"]*z )";
      yp <<"(-pars["<<parIndex+1<<"]*pars["<<parIndex+2<<"]*R*cos(theta) "
         <<"+ pars["<<parIndex+0<<"]*R*sin(theta) "
         <<"+ pars["<<parIndex+1<<"]*pars["<<parIndex+3<<"]*z )";
      zp <<"( pars["<<parIndex+3<<"]*R*cos(theta) + pars["<<parIndex+2<<"]*z )";
      rperp <<"powr( powr(fabs"<<xp.str()<<"*pars["<<parIndex+4<<"], pars["<<parIndex+7<<"]) "
            <<"+ powr(fabs"<<yp.str()<<"*pars["<<parIndex+5<<"], pars["<<parIndex+7<<"]), (float)1./pars["<<parIndex+7<<"] )";


      codeStr<<"inline float "<<name<<"( float R, float theta, float z, __constant float *pars )\n"
         <<"{\n"
         <<"const float rs = powr( powr( "<<rperp.str()<<", pars["<<parIndex+8<<"]) +"
         <<" powr( fabs"<<zp.str()<<"*pars["<<parIndex+6<<"], pars["<<parIndex+8<<"]), (float)1./pars["<<parIndex+8<<"] );\n"
         <<"return pars["<<parIndex+11<<"] * (rs <= pars["<<parIndex+9<<"] ? 1 : exp( -pown( (rs - pars["<<parIndex+9<<"])*pars["<<parIndex+10<<"], 2)))"
         <<" * pown( (float)2.0/(exp(rs) + exp(-rs)), 2);\n"
         <<"}\n";

      return codeStr.str();

   }

   size_t FreudenreichBarProfile::getOpenCLNPars() const
   {
      return 12;
   }

   std::vector<cl_float> FreudenreichBarProfile::getOpenCLPars() const
   {
      std::vector<cl_float> pars(12);

      pars[0] = cosPhiOffset;
      pars[1] = sinPhiOffset;
      pars[2] = cosPitchAngle;
      pars[3] = sinPitchAngle;
      pars[4] = oneOver_barX;
      pars[5] = oneOver_barY;
      pars[6] = oneOver_barZ;
      pars[7] = barPerp;
      pars[8] = barPara;
      pars[9] = barREnd;
      pars[10] = oneOver_barHEnd;
      pars[11] = norm;

      return pars;
   }

#endif

   double FreudenreichBarProfile::operator() (double r, double theta, double z ) const
   {
     const double x = r*std::cos(theta), y = r*std::sin(theta);
     // Coordinates rotated into bar coordinate frame -- need rotation through
     // the phiOffset and bar pitch angle. Rotate first through y-axis, then 
     // about z-axis. The matrix elements for the z-rotation are changed for
     // the off-diagonal elements because of how phiOffset is defined.
     const double xp = cosPhiOffset*cosPitchAngle*x + sinPhiOffset*y - cosPhiOffset*sinPitchAngle*z;
     const double yp = -sinPhiOffset*cosPitchAngle*x + cosPhiOffset*y + sinPhiOffset*sinPitchAngle*z;
     const double zp = sinPitchAngle*x + cosPitchAngle*z;
 
     const double rPerp = pow( pow(fabs(xp)*oneOver_barX, barPerp) + 
			       pow(fabs(yp)*oneOver_barY, barPerp), oneOver_barPerp);
     
     const double rs = pow( pow(rPerp, barPara) + pow(fabs(zp)*oneOver_barZ, barPara), oneOver_barPara);
     
     return norm * (rs <= barREnd ? 1 : exp(-pow((rs - barREnd)*oneOver_barHEnd, 2.))) * sech2(rs);
     
   }

  void FreudenreichBarProfile::setphiOffset(double phi) {
    cosPhiOffset = std::cos(phi);
    sinPhiOffset = std::sin(phi);
  }

  void FreudenreichBarProfile::setpitchAngle(double pa) {
    sinPitchAngle = std::sin(pa);
    cosPitchAngle = std::cos(pa);
  }
  
  void FreudenreichBarProfile::setbarX(double q) {
    oneOver_barX = 1./q;
  }
  
  void FreudenreichBarProfile::setbarY(double q) {
    oneOver_barY = 1./q;
  }
  
  void FreudenreichBarProfile::setbarZ(double q) {
    oneOver_barZ = 1./q;
  }
  
  void FreudenreichBarProfile::setbarPerp(double q) {
    barPerp = q;
    oneOver_barPerp = 1./q;
  }
  
  void FreudenreichBarProfile::setbarPara(double q) {
    barPara = q;
    oneOver_barPara = 1./q;
  }
  
  void FreudenreichBarProfile::setbarREnd(double i) {
    barREnd = i;
  }
  
  void FreudenreichBarProfile::setbarHEnd(double i) {
    oneOver_barHEnd = 1./i;
  }
  
  void FreudenreichBarProfile::setnorm(double a) {
    norm = a;
  }
  
  /*static utl::Registry1<CylindricalProfile, const utl::Branch&>::Registrar<ExpDiscWithEccentricHoleProfile> registrarExpDiscWithEccentricHole("ExpDiscWithEccentricHole");

  ExpDiscWithEccentricHoleProfile::ExpDiscWithEccentricHoleProfile ( double R0, double Rs_, double Rh, double hi_, double phiOffset, double holeEccentricity, double Zs, double Zoffset, double norm, const std::string &R0name, const std::string &Rsname, const std::string &Rhname, const std::string &hiname, const std::string &phiOffsetname, const std::string &holeEccentricityname, const std::string &Zsname, const std::string &Zoffsetname, const std::string &normname ) :
    fR0name(R0name),
    fRsname(Rsname),
    fRhname(Rhname),
    fhiname(hiname),
    fphiOffsetname(phiOffsetname),
    feccentricityname(holeEccentricityname),
    fZsname(Zsname),
    fZoffsetname(Zoffsetname),
    fnormname(normname)
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
     fvars.add(holeEccentricityname, holeEccentricity);
     fvars.add(Zsname, Zs);
     oneOver_Zs = 1./Zs;
     fvars.add(Zoffsetname, Zoffset);
     fvars.add(normname, norm);
     
     updateVariableValues( fvars );
   }

  ExpDiscWithEccentricHoleProfile::ExpDiscWithEccentricHoleProfile ( const utl::Parameters &pars, const std::string &R0name, const std::string &Rsname, const std::string &Rhname, const std::string &hiname, const std::string &phiOffsetname, const std::string &holeEccentricityname, const std::string &Zsname, const std::string &Zoffsetname, const std::string &normname) :
    fR0name(R0name),
    fRsname(Rsname),
    fRhname(Rhname),
    fhiname(hiname),
    fphiOffsetname(phiOffsetname),
    feccentricityname(holeEccentricityname),
    fZsname(Zsname),
    fZoffsetname(Zoffsetname),
    fnormname(normname)
   {
      std::vector<std::string> varNames(9);
      varNames[0] = R0name;
      varNames[1] = Rsname;
      varNames[2] = Rhname;
      varNames[3] = hiname;
      varNames[4] = phiOffsetname;
      varNames[5] = holeEccentricityname;
      varNames[6] = Zsname;
      varNames[7] = Zoffsetname;
      varNames[8] = normname;

      fvars.add(varNames, pars);

      updateVariableValues( fvars );
   }

   ExpDiscWithEccentricHoleProfile::ExpDiscWithEccentricHoleProfile ( const utl::Branch &b ) :
      CylindricalProfile(b)
   {
      assert(b.GetAttributes()["type"] == "ExpDiscWithEccentricHole");

      std::vector<std::string> varNames(9);
      varNames[0] = "R0";
      varNames[1] = "Rs";
      varNames[2] = "Rh";
      varNames[3] = "hi";
      varNames[4] = "phiOffset";
      varNames[5] = "eccentricity";
      varNames[6] = "Zs";
      varNames[7] = "Zoffset";
      varNames[8] = "norm";

      fvars.add(varNames, profileName, b);

      fR0name = varNames[0];
      fRsname = varNames[1];
      fRhname = varNames[2];
      fhiname = varNames[3];
      fphiOffsetname = varNames[4];
      feccentricityname = varNames[5];
      fZsname = varNames[6];
      fZoffsetname = varNames[7];
      fnormname = varNames[8];

      updateVariableValues( fvars );
   }

   void ExpDiscWithEccentricHoleProfile::updateVariableValues( const utl::Variables &vars )
   {
      fvars.setFromVariables(vars);

      setR0(fvars[fR0name]);
      setRs(fvars[fRsname]);
      setRh(fvars[fRhname]);
      sethi(fvars[fhiname]);
      setphiOffset(fvars[fphiOffsetname]);
      setholeEccentricity(fvars[feccentricityname]);
      setZs(fvars[fZsname]);
      setZoffset(fvars[fZoffsetname]);
      setnorm(fvars[fnormname]);

   }

   double ExpDiscWithEccentricHoleProfile::operator() (double r, double theta, double z ) const
   {
     const auto x = r*std::cos(theta), y = r*std::sin(theta),
       xp = cosPhiOffset*x + sinPhiOffset*y,
       yp = -sinPhiOffset*x + cosPhiOffset*y;

     const auto rhole = std::sqrt(xp*xp + eccentricity*eccentricity*yp*yp);

     return norm*std::exp(-(r-Rs)*oneOver_R0) * (1 - std::exp(-std::pow(rhole*oneOver_Rh, hi)))*std::exp(-std::fabs(z - Zoffset)*oneOver_Zs);

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
  
  void ExpDiscWithEccentricHoleProfile::setholeEccentricity(double q) {
    eccentricity = q;
  }
  
  void ExpDiscWithEccentricHoleProfile::setZs(double q) {
    oneOver_Zs = 1./q;
  }
  
  void ExpDiscWithEccentricHoleProfile::setZoffset(double q) {
    Zoffset = q;
  }
    
  void ExpDiscWithEccentricHoleProfile::setnorm(double a) {
    norm = a;
  }
  */  


   static utl::Registry1<CylindricalProfile, const utl::Branch&>::Registrar<FerriereCMZProfile> registrarFerriereCMZ("FerriereCMZ");

   FerriereCMZProfile::FerriereCMZProfile ( double phiOffset, double Xmax, double Hc, double x0, double y0, double norm, 
            const std::string &phiOffsetname, const std::string &Xmaxname, const std::string &Hcname,
            const std::string &x0name, const std::string &y0name, const std::string &normname ) :
      Yc(2.5),
      fphiOffsetname(phiOffsetname),
      fXmaxname(Xmaxname),
      fHcname(Hcname),
      fx0name(x0name),
      fy0name(y0name),
      fnormname(normname)
   {
      fvars.add(phiOffsetname, phiOffset);
      fvars.add(Xmaxname, Xmax);
      fvars.add(Hcname, Hc);
      fvars.add(x0name, x0);
      fvars.add(y0name, y0);
      fvars.add(normname, norm);

      updateVariableValues( fvars );
   }

   FerriereCMZProfile::FerriereCMZProfile ( double phiOffset, double Yc, double Xc, double Lc, double Hc, double x0, double y0, double norm, 
         const std::string &phiOffsetname, const std::string &Ycname, const std::string &Xcname, const std::string &Lcname,
         const std::string &Hcname, const std::string &x0name, const std::string &y0name, const std::string &normname ) :
      fphiOffsetname(phiOffsetname),
      fXmaxname(""),
      fYcname(Ycname),
      fXcname(Xcname),
      fLcname(Lcname),
      fHcname(Hcname),
      fx0name(x0name),
      fy0name(y0name),
      fnormname(normname)
   {
      fvars.add(phiOffsetname, phiOffset);
      fvars.add(Ycname, Yc);
      fvars.add(Xcname, Xc);
      fvars.add(Lcname, Lc);
      fvars.add(Hcname, Hc);
      fvars.add(x0name, x0);
      fvars.add(y0name, y0);
      fvars.add(normname, norm);

      updateVariableValues( fvars );
   }

   FerriereCMZProfile::FerriereCMZProfile ( const utl::Parameters &pars, 
            const std::string &phiOffsetname, const std::string &Xmaxname, const std::string &Hcname,
            const std::string &x0name, const std::string &y0name, const std::string &normname ) :
      Yc(2.5),
      fphiOffsetname(phiOffsetname),
      fXmaxname(Xmaxname),
      fHcname(Hcname),
      fx0name(x0name),
      fy0name(y0name),
      fnormname(normname)
   {
      std::vector<std::string> varNames(6);
      varNames[0] = phiOffsetname;
      varNames[1] = Xmaxname;
      varNames[2] = Hcname;
      varNames[3] = x0name;
      varNames[4] = y0name;
      varNames[5] = normname;

      fvars.add(varNames, pars);

      updateVariableValues( fvars );
   }

   FerriereCMZProfile::FerriereCMZProfile ( const utl::Parameters &pars, 
         const std::string &phiOffsetname, const std::string &Ycname, const std::string &Xcname, const std::string &Lcname,
         const std::string &Hcname, const std::string &x0name, const std::string &y0name, const std::string &normname ) :
      fphiOffsetname(phiOffsetname),
      fXmaxname(""),
      fYcname(Ycname),
      fXcname(Xcname),
      fLcname(Lcname),
      fHcname(Hcname),
      fx0name(x0name),
      fy0name(y0name),
      fnormname(normname)
   {
      std::vector<std::string> varNames(8);
      varNames[0] = phiOffsetname;
      varNames[1] = Ycname;
      varNames[2] = Xcname;
      varNames[3] = Lcname;
      varNames[4] = Hcname;
      varNames[5] = x0name;
      varNames[6] = y0name;
      varNames[7] = normname;

      fvars.add(varNames, pars);

      updateVariableValues( fvars );
   }

   FerriereCMZProfile::FerriereCMZProfile ( const utl::Branch &b ) :
      CylindricalProfile(b),
      fXmaxname("")
   {
      assert(b.GetAttributes()["type"] == "FerriereCMZ");

      std::vector<std::string> varNames(8);
      varNames[0] = "phiOffset";
      varNames[1] = "Yc";
      varNames[2] = "Xc";
      varNames[3] = "Lc";
      varNames[4] = "Hc";
      varNames[5] = "x0";
      varNames[6] = "y0";
      varNames[7] = "norm";

      fvars.add(varNames, profileName, b);

      fphiOffsetname = varNames[0];
      fYcname = varNames[1];
      fXcname = varNames[2];
      fLcname = varNames[3];
      fHcname = varNames[4];
      fx0name = varNames[5];
      fy0name = varNames[6];
      fnormname = varNames[7];

      updateVariableValues( fvars );
   }

   void FerriereCMZProfile::updateVariableValues( const utl::Variables &vars )
   {
      fvars.setFromVariables(vars);

      setphiOffset(fvars[fphiOffsetname]);
      setHc(fvars[fHcname]);
      setx0(fvars[fx0name]);
      sety0(fvars[fy0name]);
      setnorm(fvars[fnormname]);

      if (fXmaxname != "") {
         setXmax(fvars[fXmaxname]);
      } else {
         setYc(fvars[fYcname]);
         setXc(fvars[fXcname]);
         setLc(fvars[fLcname]);
      }
   }

   void FerriereCMZProfile::addToDOM( xercesc::DOMNode *node ) const 
   {
      xercesc::DOMDocument *doc = node->getOwnerDocument();

      xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("CylindricalProfile").unicodeForm());
      node->appendChild(profileEl);
      profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("FerriereCMZ").unicodeForm());

      //Add the name
      {
         xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
         profileEl->appendChild(element);

         xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
         element->appendChild(text);
      }

      std::vector<std::string> varIds(8), varNames(8);

      varIds[0] = "phiOffset";
      varIds[1] = "Yc";
      varIds[2] = "Xc";
      varIds[3] = "Lc";
      varIds[4] = "Hc";
      varIds[5] = "x0";
      varIds[6] = "y0";
      varIds[7] = "norm";

      varNames[0] = fphiOffsetname;
      varNames[1] = fYcname;
      varNames[2] = fXcname;
      varNames[3] = fLcname;
      varNames[4] = fHcname;
      varNames[5] = fx0name;
      varNames[6] = fy0name;
      varNames[7] = fnormname;

      fvars.addToDOM( profileEl, varNames, varIds, profileName );
   }

#ifdef HAVE_OPENCL
   std::string FerriereCMZProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
   {
      //Update the name, use this name for the other functions
      name = "FerriereCMZ_"+profileName+"_"+name;

      //Use an ostringstream object to store the function 
      std::ostringstream codeStr;
      
      //To reduce the number of registers, put calculation of rs into one big line.
      std::ostringstream xp,yp;
      xp <<"( pars["<<parIndex+0<<"]*(R*cos(theta)-pars["<<parIndex+2<<"]) "
         <<"- pars["<<parIndex+1<<"]*(R*sin(theta)+pars["<<parIndex+3<<"]) )";
      yp <<"( - pars["<<parIndex+1<<"]*(R*cos(theta)-pars["<<parIndex+2<<"]) "
         <<"- pars["<<parIndex+0<<"]*(R*sin(theta)+pars["<<parIndex+3<<"]) )";

      codeStr<<"inline float "<<name<<"( float R, float theta, float z, __constant float *pars )\n"
         <<"{\n"
         <<"return pars["<<parIndex+8<<"] * exp( -pown( (sqrt( pown("<<xp.str()<<",2) + pown(pars["<<parIndex+4<<"]*"<<yp.str()<<",2) )"
         <<" - pars["<<parIndex+5<<"]) * pars["<<parIndex+6<<"], 4) - pown(z,2)*pars["<<parIndex+7<<"]);\n"
         <<"}\n";

      return codeStr.str();

   }

   size_t FerriereCMZProfile::getOpenCLNPars() const
   {
      return 9;
   }

   std::vector<cl_float> FerriereCMZProfile::getOpenCLPars() const
   {
      std::vector<cl_float> pars(9);

      pars[0] = cosPhiOffset;
      pars[1] = sinPhiOffset;
      pars[2] = x0;
      pars[3] = y0;
      pars[4] = Yc;
      pars[5] = Xc;
      pars[6] = oneOver_Lc;
      pars[7] = oneOver_Hc2;
      pars[8] = norm;

      return pars;
   }

#endif

   double FerriereCMZProfile::operator() (double r, double theta, double z ) const
   {
      //Note that Ferriere uses a left handed system, so y => -y
      const double x = r*cos(theta), y = r*sin(theta),
            xp = cosPhiOffset*(x-x0) - sinPhiOffset*(y+y0),
            yp = -sinPhiOffset*(x-x0) - cosPhiOffset*(y+y0);

      return norm * exp(-pow( (sqrt(xp*xp + Yc*Yc*yp*yp) - Xc)*oneOver_Lc, 4 ) - z*z*oneOver_Hc2);
   }

   void FerriereCMZProfile::setphiOffset(double phi) {
      cosPhiOffset = cos(phi);
      sinPhiOffset = sin(phi);
   }

   void FerriereCMZProfile::setXmax(double q) {
      Xc = q/2.;
      oneOver_Lc = 2*pow(log(2.), 0.25)/q;
   }

   void FerriereCMZProfile::setYc(double q) {
      Yc = q;
   }

   void FerriereCMZProfile::setXc(double q) {
      Xc = q;
   }

   void FerriereCMZProfile::setLc(double q) {
      oneOver_Lc = 1./q;
   }

   void FerriereCMZProfile::setHc(double q) {
      oneOver_Hc2 = 1./(q*q);
   }

   void FerriereCMZProfile::setx0(double q) {
      x0 = q;
   }

   void FerriereCMZProfile::sety0(double q) {
      y0 = q;
   }

   void FerriereCMZProfile::setnorm(double a) {
      norm = a;
   }




   static utl::Registry1<CylindricalProfile, const utl::Branch&>::Registrar<FerriereDiskProfile> registrarFerriereDisk("FerriereDisk");

   FerriereDiskProfile::FerriereDiskProfile ( double alpha, double beta, double phid, double Xmax, double Xmin, double Hd, double norm, 
            const std::string &alphaname, const std::string &betaname, const std::string &phidname, 
            const std::string &Xmaxname, const std::string &Xminname, const std::string &Hdname, const std::string &normname ) :
      Xmax(0),
      Xmin(0),
      falphaname(alphaname),
      fbetaname(betaname),
      fphidname(phidname),
      fXmaxname(Xmaxname),
      fXminname(Xminname),
      fYdname("FerriereDisk_Yd"),
      fHdname(Hdname),
      fnormname(normname)
   {
      fvars.add(alphaname, alpha);
      fvars.add(betaname, beta);
      fvars.add(phidname, phid);
      fvars.add(Xmaxname, Xmax);
      fvars.add(Xminname, Xmin);
      fvars.add(fYdname, 3.1);
      fvars.add(Hdname, Hd);
      fvars.add(normname, norm);

      updateVariableValues( fvars );
   }

   FerriereDiskProfile::FerriereDiskProfile ( double alpha, double beta, double phid, double Xmax, double Xmin, double Yd, double Hd, double norm, 
            const std::string &alphaname, const std::string &betaname, const std::string &phidname, const std::string &Xmaxname, 
            const std::string &Xminname, const std::string &Ydname, const std::string &Hdname, const std::string &normname ) :
      Xmax(0),
      Xmin(0),
      falphaname(alphaname),
      fbetaname(betaname),
      fphidname(phidname),
      fXmaxname(Xmaxname),
      fXminname(Xminname),
      fYdname(Ydname),
      fHdname(Hdname),
      fnormname(normname)
   {
      fvars.add(alphaname, alpha);
      fvars.add(betaname, beta);
      fvars.add(phidname, phid);
      fvars.add(Xmaxname, Xmax);
      fvars.add(Xminname, Xmin);
      fvars.add(Ydname, Yd);
      fvars.add(Hdname, Hd);
      fvars.add(normname, norm);

      updateVariableValues( fvars );
   }

   FerriereDiskProfile::FerriereDiskProfile ( const utl::Parameters &pars, 
            const std::string &alphaname, const std::string &betaname, const std::string &phidname, 
            const std::string &Xmaxname, const std::string &Xminname, const std::string &Hdname, const std::string &normname ) :
      Xmax(0),
      Xmin(0),
      falphaname(alphaname),
      fbetaname(betaname),
      fphidname(phidname),
      fXmaxname(Xmaxname),
      fXminname(Xminname),
      fYdname("FerriereDisk_Yd"),
      fHdname(Hdname),
      fnormname(normname)
   {
      std::vector<std::string> varNames(7);
      varNames[0] = alphaname;
      varNames[1] = betaname;
      varNames[2] = phidname;
      varNames[3] = Xmaxname;
      varNames[4] = Xminname;
      varNames[5] = Hdname;
      varNames[6] = normname;

      fvars.add(varNames, pars);

      fvars.add(fYdname, 3.1);

      updateVariableValues( fvars );
   }

   FerriereDiskProfile::FerriereDiskProfile ( const utl::Parameters &pars, 
            const std::string &alphaname, const std::string &betaname, const std::string &phidname, const std::string &Xmaxname, 
            const std::string &Xminname, const std::string &Ydname, const std::string &Hdname, const std::string &normname ) :
      Xmax(0),
      Xmin(0),
      falphaname(alphaname),
      fbetaname(betaname),
      fphidname(phidname),
      fXmaxname(Xmaxname),
      fXminname(Xminname),
      fYdname(Ydname),
      fHdname(Hdname),
      fnormname(normname)
   {
      std::vector<std::string> varNames(8);
      varNames[0] = alphaname;
      varNames[1] = betaname;
      varNames[2] = phidname;
      varNames[3] = Xmaxname;
      varNames[4] = Xminname;
      varNames[5] = Ydname;
      varNames[6] = Hdname;
      varNames[7] = normname;

      fvars.add(varNames, pars);

      updateVariableValues( fvars );
   }

   FerriereDiskProfile::FerriereDiskProfile ( const utl::Branch &b ) :
      CylindricalProfile(b),
      Xmax(0),
      Xmin(0)
   {
      assert(b.GetAttributes()["type"] == "FerriereDisk");

      std::vector<std::string> varNames(8);
      varNames[0] = "alpha";
      varNames[1] = "beta";
      varNames[2] = "phid";
      varNames[3] = "Xmax";
      varNames[4] = "Xmin";
      varNames[5] = "Yd";
      varNames[6] = "Hd";
      varNames[7] = "norm";

      fvars.add(varNames, profileName, b);

      falphaname = varNames[0];
      fbetaname = varNames[1];
      fphidname = varNames[2];
      fXmaxname = varNames[3];
      fXminname = varNames[4];
      fYdname = varNames[5];
      fHdname = varNames[6];
      fnormname = varNames[7];

      updateVariableValues( fvars );
   }

   void FerriereDiskProfile::updateVariableValues( const utl::Variables &vars )
   {
      fvars.setFromVariables(vars);

      setalpha(fvars[falphaname]);
      setbeta(fvars[fbetaname]);
      setphid(fvars[fphidname]);
      setXmax(fvars[fXmaxname]);
      setXmin(fvars[fXminname]);
      setYd(fvars[fYdname]);
      setHd(fvars[fHdname]);
      setnorm(fvars[fnormname]);

   }

   void FerriereDiskProfile::addToDOM( xercesc::DOMNode *node ) const 
   {
      xercesc::DOMDocument *doc = node->getOwnerDocument();

      xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("CylindricalProfile").unicodeForm());
      node->appendChild(profileEl);
      profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("FerriereDisk").unicodeForm());

      //Add the name
      {
         xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
         profileEl->appendChild(element);

         xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
         element->appendChild(text);
      }

      std::vector<std::string> varIds(8), varNames(8);

      varIds[0] = "alpha";
      varIds[1] = "beta";
      varIds[2] = "phid";
      varIds[3] = "Xmax";
      varIds[4] = "Xmin";
      varIds[5] = "Yd";
      varIds[6] = "Hd";
      varIds[7] = "norm";

      varNames[0] = falphaname;
      varNames[1] = fbetaname;
      varNames[2] = fphidname;
      varNames[3] = fXmaxname;
      varNames[4] = fXminname;
      varNames[5] = fYdname;
      varNames[6] = fHdname;
      varNames[7] = fnormname;

      fvars.addToDOM( profileEl, varNames, varIds, profileName );
   }

#ifdef HAVE_OPENCL
   std::string FerriereDiskProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
   {
      //Update the name, use this name for the other functions
      name = "FerriereDisk_"+profileName+"_"+name;

      //Use an ostringstream object to store the function 
      std::ostringstream codeStr;
      
      //To reduce the number of registers, put calculation of rs into one big line.
      std::ostringstream xp,yp,zp;
      xp <<"( pars["<<parIndex+0<<"]*pars["<<parIndex+2<<"]*R*cos(theta) "
         <<"+ (pars["<<parIndex+1<<"]*pars["<<parIndex+2<<"]*pars["<<parIndex+5<<"]-pars["<<parIndex+3<<"]*pars["<<parIndex+4<<"])*R*sin(theta)"
         <<"- (pars["<<parIndex+1<<"]*pars["<<parIndex+2<<"]*pars["<<parIndex+4<<"]+pars["<<parIndex+3<<"]*pars["<<parIndex+5<<"])*z )";
      yp <<"( -pars["<<parIndex+0<<"]*pars["<<parIndex+3<<"]*R*cos(theta) "
         <<"- (pars["<<parIndex+1<<"]*pars["<<parIndex+3<<"]*pars["<<parIndex+5<<"]+pars["<<parIndex+2<<"]*pars["<<parIndex+4<<"])*R*sin(theta)"
         <<"+ (pars["<<parIndex+1<<"]*pars["<<parIndex+3<<"]*pars["<<parIndex+4<<"]-pars["<<parIndex+2<<"]*pars["<<parIndex+5<<"])*z )";
      zp <<"( pars["<<parIndex+1<<"]*R*cos(theta) "
         <<"- pars["<<parIndex+0<<"]*pars["<<parIndex+5<<"]*R*sin(theta)"
         <<"+ pars["<<parIndex+0<<"]*pars["<<parIndex+4<<"]*z )";

      codeStr<<"inline float "<<name<<"( float R, float theta, float z, __constant float *pars )\n"
         <<"{\n"
         <<"return pars["<<parIndex+10<<"] * exp( -pown( (sqrt( pown("<<xp.str()<<",2) + pown(pars["<<parIndex+6<<"]*"<<yp.str()<<",2) )"
         <<" - pars["<<parIndex+7<<"]) * pars["<<parIndex+8<<"], 4) - pown("<<zp.str()<<",2)*pars["<<parIndex+9<<"]);\n"
         <<"}\n";

      return codeStr.str();

   }

   size_t FerriereDiskProfile::getOpenCLNPars() const
   {
      return 11;
   }

   std::vector<cl_float> FerriereDiskProfile::getOpenCLPars() const
   {
      std::vector<cl_float> pars(11);

      pars[0] = cosBeta;
      pars[1] = sinBeta;
      pars[2] = cosPhid;
      pars[3] = sinPhid;
      pars[4] = cosAlpha;
      pars[5] = sinAlpha;
      pars[6] = Yd;
      pars[7] = Xd;
      pars[8] = oneOver_Ld;
      pars[9] = oneOver_Hd2;
      pars[10] = norm;

      return pars;
   }

#endif

   double FerriereDiskProfile::operator() (double r, double theta, double z ) const
   {
      //Ferriere uses a left handed system so y => -y
      const double x = r*cos(theta), y = r*sin(theta),
            xp =  x*cosBeta*cosPhid + y*(sinAlpha*sinBeta*cosPhid - cosAlpha*sinPhid) - z*(cosAlpha*sinBeta*cosPhid + sinAlpha*sinPhid),
            yp = -x*cosBeta*sinPhid - y*(sinAlpha*sinBeta*sinPhid + cosAlpha*cosPhid) + z*(cosAlpha*sinBeta*sinPhid - sinAlpha*cosPhid),
            zp =  x*sinBeta         - y*(sinAlpha*cosBeta                           ) + z*(cosAlpha*cosBeta                           );
            
      return norm * exp(-pow( (sqrt(xp*xp + Yd*Yd*yp*yp) - Xd)*oneOver_Ld, 4 ) - zp*zp*oneOver_Hd2);
   }

   void FerriereDiskProfile::setalpha(double phi) {
      cosAlpha = cos(phi);
      sinAlpha = sin(phi);
   }

   void FerriereDiskProfile::setbeta(double phi) {
      cosBeta = cos(phi);
      sinBeta = sin(phi);
   }

   void FerriereDiskProfile::setphid(double phi) {
      cosPhid = cos(phi);
      sinPhid = sin(phi);
   }

   void FerriereDiskProfile::setXmax(double q) {
      Xmax = q;
      setXdLd();
   }

   void FerriereDiskProfile::setXmin(double q) {
      Xmin = q;
      setXdLd();
   }

   void FerriereDiskProfile::setHd(double q) {
      oneOver_Hd2 = 1./(q*q);
   }

   void FerriereDiskProfile::setYd(double q) {
      Yd = q;
   }

   void FerriereDiskProfile::setnorm(double a) {
      norm = a;
   }

   void FerriereDiskProfile::setXdLd() {
      Xd = (Xmax+Xmin)/2.;
      oneOver_Ld = 2*pow(log(2.), 0.25)/(Xmax-Xmin);
   }





   static utl::Registry1<CylindricalProfile, const utl::Branch&>::Registrar<GenericBarProfile> registrarGenericBar("GenericBar");

   GenericBarProfile::GenericBarProfile ( double phi, double alpha, 
         double ySc, double rSc, double zSc, 
         double ri, double zi, double x0, double z0,
         double ei, double pi, double norm, 
         const std::string &phiname, const std::string& alphaname, 
         const std::string &yScname, const std::string &rScname, const std::string &zScname, 
         const std::string &riname, const std::string &ziname, const std::string &x0name, const std::string &z0name,
         const std::string &einame, const std::string &piname, const std::string &normname ) :
      fphiname(phiname),
      falphaname(alphaname),
      fyScname(yScname),
      frScname(rScname),
      fzScname(zScname),
      friname(riname),
      fziname(ziname),
      fx0name(x0name),
      fz0name(z0name),
      feiname(einame),
      fpiname(piname),
      fnormname(normname)
   {
      fvars.add(phiname, phi);
      fvars.add(alphaname, alpha);
      fvars.add(yScname, ySc);
      fvars.add(rScname, rSc);
      fvars.add(zScname, zSc);
      fvars.add(riname, ri);
      fvars.add(ziname, zi);
      fvars.add(x0name, x0);
      fvars.add(z0name, z0);
      fvars.add(einame, ei);
      fvars.add(piname, pi);
      fvars.add(normname, norm);
    
      updateVariableValues( fvars );
   }

   GenericBarProfile::GenericBarProfile ( const utl::Parameters &pars, 
         const std::string &phiname, const std::string& alphaname, 
         const std::string &yScname, const std::string &rScname, const std::string &zScname, 
         const std::string &riname, const std::string &ziname, const std::string &x0name, const std::string &z0name,
         const std::string &einame, const std::string &piname, const std::string &normname ) :
      fphiname(phiname),
      falphaname(alphaname),
      fyScname(yScname),
      frScname(rScname),
      fzScname(zScname),
      friname(riname),
      fziname(ziname),
      fx0name(x0name),
      fz0name(z0name),
      feiname(einame),
      fpiname(piname),
      fnormname(normname)
   {
      std::vector<std::string> varNames(12);

      varNames[0] = phiname;
      varNames[1] = alphaname;
      varNames[2] = yScname;
      varNames[3] = rScname;
      varNames[4] = zScname;
      varNames[5] = riname;
      varNames[6] = ziname;
      varNames[7] = x0name;
      varNames[8] = z0name;
      varNames[9] = einame;
      varNames[10] = piname;
      varNames[11] = normname;

      fvars.add(varNames, pars);

      updateVariableValues( fvars );
   }

   GenericBarProfile::GenericBarProfile ( const utl::Branch &b ) :
      CylindricalProfile(b)
   {
      assert(b.GetAttributes()["type"] == "GenericBar");

      std::vector<std::string> varNames(12);
      varNames[0] = "phi";
      varNames[1] = "alpha";
      varNames[2] = "ySc";
      varNames[3] = "rSc";
      varNames[4] = "zSc";
      varNames[5] = "ri";
      varNames[6] = "zi";
      varNames[7] = "x0";
      varNames[8] = "z0";
      varNames[9] = "ei";
      varNames[10] = "pi";
      varNames[11] = "norm";

      fvars.add(varNames, profileName, b);

      fphiname = varNames[0];
      falphaname = varNames[1];
      fyScname = varNames[2];
      frScname = varNames[3];
      fzScname = varNames[4];
      friname = varNames[5];
      fziname = varNames[6];
      fx0name = varNames[7];
      fz0name = varNames[8];
      feiname = varNames[9];
      fpiname = varNames[10];
      fnormname = varNames[11];

      updateVariableValues( fvars );
   }

   void GenericBarProfile::updateVariableValues( const utl::Variables &vars )
   {
      fvars.setFromVariables(vars);

      setphi(fvars[fphiname]);
      setalpha(fvars[falphaname]);
      setySc(fvars[fyScname]);
      setrSc(fvars[frScname]);
      setzSc(fvars[fzScname]);
      setri(fvars[friname]);
      setzi(fvars[fziname]);
      setx0(fvars[fx0name]);
      setz0(fvars[fz0name]);
      setei(fvars[feiname]);
      setpi(fvars[fpiname]);
      setnorm(fvars[fnormname]);

   }

   void GenericBarProfile::addToDOM( xercesc::DOMNode *node ) const 
   {
      xercesc::DOMDocument *doc = node->getOwnerDocument();

      xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("CylindricalProfile").unicodeForm());
      node->appendChild(profileEl);
      profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("GenericBar").unicodeForm());

      //Add the name
      {
         xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
         profileEl->appendChild(element);

         xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
         element->appendChild(text);
      }

      std::vector<std::string> varIds(12), varNames(12);

      varIds[0] = "phi";
      varIds[1] = "alpha";
      varIds[2] = "ySc";
      varIds[3] = "rSc";
      varIds[4] = "zSc";
      varIds[5] = "ri";
      varIds[6] = "zi";
      varIds[7] = "x0";
      varIds[8] = "z0";
      varIds[9] = "ei";
      varIds[10] = "pi";
      varIds[11] = "norm";

      varNames[0] = fphiname;
      varNames[1] = falphaname;
      varNames[2] = fyScname;
      varNames[3] = frScname;
      varNames[4] = fzScname;
      varNames[5] = friname;
      varNames[6] = fziname;
      varNames[7] = fx0name;
      varNames[8] = fz0name;
      varNames[9] = feiname;
      varNames[10] = fpiname;
      varNames[11] = fnormname;

      fvars.addToDOM( profileEl, varNames, varIds, profileName );
   }

#ifdef HAVE_OPENCL
   std::string GenericBarProfile::getOpenCLFunction(std::string &name, size_t parIndex) const
   {
      //Update the name, use this name for the other functions
      name = "GenericBar_"+profileName+"_"+name;

      //Use an ostringstream object to store the function 
      std::ostringstream codeStr;

      //To reduce the number of registers, put calculation of rr into one big line.
      std::ostringstream xpp,xp,yp,zp,rp;
      xpp<<"( R*cos(theta-pars["<<parIndex+0<<"]) )";
      yp <<"( R*sin(theta-pars["<<parIndex+0<<"]) )";
      xp <<"( pars["<<parIndex+1<<"]*"<<xpp.str()<<" + pars["<<parIndex+2<<"]*z - pars["<<parIndex+10<<"] )";
      zp <<"(-pars["<<parIndex+2<<"]*"<<xpp.str()<<" + pars["<<parIndex+1<<"]*z - pars["<<parIndex+11<<"] )";
      rp <<"powr( powr(fabs"<<yp.str()<<"*pars["<<parIndex+3<<"], pars["<<parIndex+6<<"]) "
         <<"+ powr(fabs"<<xp.str()<<", pars["<<parIndex+6<<"]), pars["<<parIndex+8<<"] )";


      codeStr<<"inline float "<<name<<"( float R, float theta, float z, __constant float *pars )\n"
         <<"{\n"
         <<"const float rr = powr( powr( "<<rp.str()<<"*pars["<<parIndex+4<<"], pars["<<parIndex+7<<"]) +"
         <<" powr( fabs"<<zp.str()<<"*pars["<<parIndex+5<<"], pars["<<parIndex+7<<"]), pars["<<parIndex+9<<"] );\n"
         <<"return pars["<<parIndex+14<<"] * exp( -powr( rr, pars["<<parIndex+12<<"]) ) * powr( rr, pars["<<parIndex+13<<"]);\n"
         <<"}\n";

      return codeStr.str();

   }

   size_t GenericBarProfile::getOpenCLNPars() const
   {
      return 15;
   }

   std::vector<cl_float> GenericBarProfile::getOpenCLPars() const
   {
      std::vector<cl_float> pars(15);

      pars[0] = fphi;
      pars[1] = fcosAlpha;
      pars[2] = fsinAlpha;
      pars[3] = fOO_ySc;
      pars[4] = fOO_rSc;
      pars[5] = fOO_zSc;
      pars[6] = fri;
      pars[7] = fzi;
      pars[8] = fOO_ri;
      pars[9] = fOO_zi;
      pars[10] = fx0;
      pars[11] = fz0;
      pars[12] = fei;
      pars[13] = fpi;
      pars[14] = fnorm;

      return pars;
   }

#endif

   double GenericBarProfile::operator() (double r, double theta, double z ) const
   {
      // Rotate coordinates to bar axis.  Rotate first around z axis and then y axis.
      // Allow for an offset in xp and zp to shift the bar around on the projection.
      // The rotation around the z axiz is performed in cylindrical coordinates
      const double yp =  r*std::sin(theta-fphi);
      const double xpp = r*std::cos(theta-fphi);
      const double xp =  xpp*fcosAlpha + z*fsinAlpha - fx0;
      const double zp = -xpp*fsinAlpha + z*fcosAlpha - fz0;

      const double rp = std::pow( std::pow( std::fabs(yp)*fOO_ySc, fri ) + std::pow( std::fabs(xp), fri), fOO_ri);
      const double rr = std::pow( std::pow( rp*fOO_rSc, fzi ) + std::pow( std::fabs(zp)*fOO_zSc, fzi), fOO_zi);

      return fnorm * exp( -pow( rr, fei ) ) * pow( rr, fpi );

   }

   void GenericBarProfile::setphi(double phi) {
      fphi = phi;
   }

   void GenericBarProfile::setalpha(double pa) {
      fsinAlpha = std::sin(pa);
      fcosAlpha = std::cos(pa);
   }

   void GenericBarProfile::setySc(double q) {
      fOO_ySc = 1./q;
   }

   void GenericBarProfile::setrSc(double q) {
      fOO_rSc = 1./q;
   }

   void GenericBarProfile::setzSc(double q) {
      fOO_zSc = 1./q;
   }

   void GenericBarProfile::setri(double i) {
      fri = i;
      fOO_ri = 1./i;
   }

   void GenericBarProfile::setzi(double i) {
      fzi = i;
      fOO_zi = 1./i;
   }

   void GenericBarProfile::setx0(double l) {
      fx0 = l;
   }

   void GenericBarProfile::setz0(double l) {
      fz0 = l;
   }

   void GenericBarProfile::setei(double i) {
      fei = i;
   }

   void GenericBarProfile::setpi(double i) {
      fpi = i;
   }

   void GenericBarProfile::setnorm(double a) {
      fnorm = a;
   }
  


   static utl::Registry1<CylindricalProfile, const utl::Branch&>::Registrar<GALPROPnHIProfile> registrarGALPROPnHI("GALPROPnHI");

   GALPROPnHIProfile::GALPROPnHIProfile ( const utl::Branch &b ) :
      CylindricalProfile(b),
      Rsun(10)
   {
      assert(b.GetAttributes()["type"] == "GALPROPnHI");

      auto rb = b.GetChild("Rsun");
      if ( rb ) 
         rb.GetData(Rsun);
      
   }

   double GALPROPnHIProfile::operator() ( double Rkpc, double theta, double Zkpc ) const
   {
      int i;                                                             // Table 1 [GB76]
      static const double R[30] ={ 0.0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0,  // kpc, col.1
         6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,10.0,10.5,11.0,
         11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0},
                   Y[30] ={ .10, .13, .14, .16, .19, .25, .30, .33, .32, .31,  // nHI, cm^-3
                      .30, .37, .38, .36, .32, .29, .38, .40, .25, .23,  // (col.3)
                      .32, .36, .32, .25, .16, .10, .09, .08, .06, .00};
      double fR, fZ,fZ1=0.,fZ2=0., R1,R2=R[29], Y1,Y2=Y[29];
      const double nGB =0.33, nDL =0.57;       // cm^-3, disk density @ 4-8 kpc; [GB76], [DL90]
      const double A1=0.395,     z1=0.212/2.,  // cm^-3, kpc; Z-distribution parameters from [DL90]
            A2=0.107,     z2=0.530/2.,
            B =0.064,     zh=0.403;

      const double HIzero=1e-20;

      //Fix the scaling according to Rsun
      Rkpc *= 10./Rsun;
      Zkpc *= 8.5/Rsun;

      for (i=0; i<28; i++)  if(R[i] <= Rkpc && Rkpc <= R[i+1])  break;  //Gulli20070810  i=28 if condition never met 

      R1 = (R[i]+R[i+1])/2;   Y1 = Y[i];
      if(Rkpc < R1)
      {  
         if(i> 0)    { R2 = (R[i-1]+R[i])/2;   Y2 = Y[i-1]; }
         else        { R2 = R[0];              Y2 = Y[0];   }
      }
      else  if(i<28) { R2 = (R[i+1]+R[i+2])/2; Y2 = Y[i+1]; }

      fR = Y1 +(Y2 -Y1)/(R2 -R1)*(Rkpc -R1);                             // interpolation in R

      R2 = (R[28] +R[29]) /2;
      if(Rkpc > R2) fR = Y[28]*exp(-(Rkpc-R2)/3);                        // extrapolation in R

      // calculation of Z-dependence
      if(Rkpc <10.)                                                      // [DL90]
         fZ1 =A1*exp(-utl::kLogTwo*pow(Zkpc/z1,2))+A2*exp(-utl::kLogTwo*pow(Zkpc/z2,2))+B*exp(-fabs(Zkpc)/zh);
      if(Rkpc > 8.) fZ2=nDL*exp(-pow(Zkpc /(0.0523*exp(0.11*Rkpc)), 2)); // [C86] IMOS20010220

      if(Rkpc <= 8.) fZ = fZ1;
      else
      {   if(Rkpc >=10.) fZ = fZ2;
         else fZ = fZ1 +(fZ2 -fZ1)/2.*(Rkpc -8.);                       // interp. [DL90] & [C86]
      }
      double nHI = fZ *fR/nGB;
      nHI *= 8.5/Rsun;
      return nHI < HIzero ? HIzero : nHI;
   }

   void GALPROPnHIProfile::addToDOM( xercesc::DOMNode *node ) const
   {
      xercesc::DOMDocument *doc = node->getOwnerDocument();

      xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("CylindricalProfile").unicodeForm());
      node->appendChild(profileEl);
      profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("GALPROPnHI").unicodeForm());

      //Add the name
      {
         xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
         profileEl->appendChild(element);

         xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
         element->appendChild(text);
      }
   }


   static utl::Registry1<CylindricalProfile, const utl::Branch&>::Registrar<GALPROPnH2Profile> registrarGALPROPnH2("GALPROPnH2");

   GALPROPnH2Profile::GALPROPnH2Profile ( const utl::Branch &b ) :
      CylindricalProfile(b),
      Rsun(10)
   {
      assert(b.GetAttributes()["type"] == "GALPROPnH2");

      auto rb = b.GetChild("Rsun");
      if ( rb ) 
         rb.GetData(Rsun);
      
   }

   double GALPROPnH2Profile::operator() ( double Rkpc, double theta, double Zkpc ) const
   {
      int i;
      const double COzero(1e-40);
      double nH2_(COzero), fR,fZ0,fZh;                                              // [B88]/Table 3
      static const double R[45] ={ 0.00, 0.12, 0.24, 0.35, 0.47, 0.59, 0.71, 0.82, 0.94,
         1.06, 1.18, 1.29, 1.41, 1.53, 1.65, 1.76,             // F07 model
         2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75,       // kpc, col.1
         6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75,
         10.6 ,11.8 ,12.9 ,14.1 ,15.3 ,16.5 ,17.6 ,18.8 ,
         20.0 ,21.2 ,22.4 ,23.5 ,24.7 },                                 //Add data from [W90] Gulli 20100218
                   Y[45] ={ 9378,  5261, 2290,  95.2, 83.7,  91.8, 65.1,  54.7, 50.1,
                      43.5,  36.0,  29.2,  22.9,  16.3,  9.4,   3.7,             // F07 model converted into average surface density and the CO using average scale height
                      1.5,  3.3,  5.8,  5.5,  8.4,  9.0,  9.6,  8.6,       // CO, K km s^-1 
                      9.1,  7.9,  9.2,  7.7,  5.0,  3.6,  4.8,  1.7,
                      1.7,  0.3,  0.9,  0.7,  0.5, 0.21, 0.11, 0.07,
                      0.05,0.04,0.003,0.006,0.002},// (col.4)
                   Z0[45]={0.000,    0,    0,    0,    0,    0,    0,    0,    0,       //Ferriere's model is symmetric around z (on average at least)
                      0,    0,    0,    0,    0,    0,    0,
                      0.039,0.036,0.000,-.008,0.001,-.010,-.001,-.004,       // kpc, col.7
                      -.019,-.022,-.014,-.009,-.004,0.013,-.004,-.020,
                      0    ,    0,    0,    0,    0,    0,    0,    0,
                      0    ,    0,    0,    0,    0},
                   Zh[45]={0.015,0.015,0.015,0.068,0.082,0.082,0.082,0.082,0.082,
                      0.082,0.082,0.082,0.082,0.082,0.082,0.082,             // Using the Ferriere surface density, but assuming the Gaussian functional form for the height, similar to Levine et al. 2006, ApJ
                      0.077,0.080,0.061,0.065,0.071,0.072,0.082,0.083,       // kpc, col.10
                      0.073,0.063,0.058,0.072,0.080,0.066,0.023,0.147,
                      0.091,0.131,0.160,0.173,0.188,0.262,0.302,0.259,
                      0.235,0.221,0.235,0.235,0.235};

      Rkpc *= 10./Rsun;
      Zkpc *= 10./Rsun;

      if(Rkpc > R[43]) return nH2_;
      for (i=0; i<43; i++)  if(R[i] <= Rkpc && Rkpc <= R[i+1])  break;

      fR =  Y[i] + ( Y[i+1] - Y[i])/(R[i+1] - R[i])*(Rkpc - R[i]); 
      fZ0= Z0[i] + (Z0[i+1] -Z0[i])/(R[i+1] - R[i])*(Rkpc - R[i]);
      fZh= Zh[i] + (Zh[i+1] -Zh[i])/(R[i+1] - R[i])*(Rkpc - R[i]);
      nH2_ =  fR * exp( -log(2.)*pow( (Zkpc-fZ0)/fZh, 2 ) );
      //We clip CO rather than H2 so the average X_CO factor doesn't get biased.
      nH2_ = nH2_< COzero ? COzero: nH2_;
      nH2_ *= 10./Rsun;
      return nH2_;
   }

   void GALPROPnH2Profile::addToDOM( xercesc::DOMNode *node ) const
   {
      xercesc::DOMDocument *doc = node->getOwnerDocument();

      xercesc::DOMElement* profileEl = doc->createElement(utl::XStr("CylindricalProfile").unicodeForm());
      node->appendChild(profileEl);
      profileEl->setAttribute(utl::XStr("type").unicodeForm(), utl::XStr("GALPROPnH2").unicodeForm());

      //Add the name
      {
         xercesc::DOMElement* element = doc->createElement(utl::XStr("name").unicodeForm());
         profileEl->appendChild(element);

         xercesc::DOMText* text = doc->createTextNode(utl::XStr(profileName).unicodeForm());
         element->appendChild(text);
      }
   }

}
