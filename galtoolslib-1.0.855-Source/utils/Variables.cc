#include "Variables.h"
#include <ErrorLogger.h>
#include <ReaderErrorReporter.h>
#include <Parameters.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <cmath>

namespace utl {
//Adding variables from parameters 
void Variables::add(const std::vector<std::string> & names, const utl::Parameters &pars){

   //Set up the variables from the parameters
   for (size_t i = 0; i < names.size(); ++i){
      
      std::vector<double> parInput;
      
      pars.getParameter(names[i], parInput);
      
      switch ( parInput.size() ) {
         // A fixed parameter
         case 1:
            add(names[i], parInput[0]);
            break;
      
         // A free parameter with no bounds
         case 2:
            add(names[i], parInput[0], parInput[1]);
            break;
      
         // A free parameter with lower bound set
         case 3:
            add(names[i], parInput[0], parInput[1]);
            setLowerBound(names[i], parInput[2]);
            break;

         // A free parameter with both bounds or only upper bound set
         case 4:
            if (parInput[2] < parInput[3]) {
               // Both bound because lower bound is smaller than upper bound
               add(names[i], parInput[0], parInput[1], parInput[2], parInput[3]);
            } else {
               // Only upper bound because lower bound is larger than upper bound
               add(names[i], parInput[0], parInput[1]);
               setUpperBound(names[i], parInput[3]);
            }
            break;
      
         default:
#ifdef DEBUG
            std::cerr<<"Number of values for variable: "<<parInput.size()<<std::endl;
            for (size_t j = 0; j < parInput.size(); ++j){
               std::cerr<<parInput[j]<<" ";
            }
            std::cerr<<std::endl;
#endif
            throw(VariableError("Variable \""+names[i]+"\" has incorrect number of values.  There should be 1 to 4 values"));
      
      }

      //Try to read in the prior, set a default one if that is not possible
      try {

         std::string priorVals;
         pars.getParameter(names[i] + "_prior", priorVals, true);

         //Need to use our own interpretor of the values because they are a mix of string and numbers
         std::string priorString;
         double par1, par2;
         std::istringstream iss(priorVals);

         iss>>priorString>>par1>>par2;

         //Throw an exception if there was an error
         if ( iss.fail() || ! iss.eof() ) {
#ifdef DEBUG
            std::cerr<<priorString<<", "<<par1<<", "<<par2<<std::endl;
#endif
            throw(VariableError("Error reading prior for variable \""+names[i]+"\"."));
         }

         setPrior(names[i], stringToPrior(priorString), par1, par2);

      } catch (utl::Parameters::ParameterError) {

         variable var = getIterator(names[i])->second;
         
         //Set the prior depending on the bounds
         switch (type(names[i])) {

            case FIXED:
            case UPPBND:
            case FREE:
            case LOWBND:
               var.prior = NONE;
               break;

            case BOTHBND:
               var.prior = UNIFORM;
               var.priorPar1 = var.lBound;
               var.priorPar2 = var.uBound;
               break;
         }

      }
   }
}

void Variables::add(std::vector<std::string> &varNames, const std::string &prefix, const utl::Branch &b ) {

   for (size_t i = 0; i < varNames.size(); ++i) {

      //Throw an exception if the named variable is not found
      const utl::Branch varBr = b.GetChild("variable", varNames[i]);
      if (! varBr ) {
         std::ostringstream os;
         os << "Variable named \""<<varNames[i]<<"\" not found in branch.";
         //INFO(os.str());
         throw(VariableError(os.str()));
      }

      //Set the variable name if found
      utl::Branch br = varBr.GetChild("name");
      if ( br ) 
         br.GetData(varNames[i]);
      else
         varNames[i] = prefix + "_" + varNames[i];

      br = varBr.GetChild("value");
      if ( !br ) {
         std::ostringstream os;
         os<<"Variable \""<<varNames[i]<<"\" has no value element";
         FATAL(os.str());
         throw(VariableError("Variables should have value element in XML"));
      }

      double value(0);
      br.GetData(value);

      add(varNames[i], value);

      //Ignore min and max if step is not given as it makes no sense.
      br = varBr.GetChild("step");
      if ( br ) {
         double error(0);
         br.GetData(error);
         setError(varNames[i], error);

         br = varBr.GetChild("min");
         if ( br ) {
            double min(0);
            br.GetData(min);
            setLowerBound(varNames[i], min);
         }

         br = varBr.GetChild("max");
         if ( br ) {
            double max(0);
            br.GetData(max);
            setUpperBound(varNames[i], max);
         }
      }

      //Set prior from the attributes, need all of them
      br = varBr.GetChild("prior");
      if ( br ) {

         std::string prior;
         br.GetData(prior);

         //First check for pval1
         br = varBr.GetChild("pval1");
         if ( ! br ) {
            std::ostringstream os;
            os<<"Prior attribute of variable \""<<varNames[i]<<"\" given but not pval1.";
            WARNING(os.str());
         } else {

            double pval1(0);
            br.GetData(pval1);

            br = varBr.GetChild("pval2");
            if ( ! br ) {
               std::ostringstream os;
               os<<"Prior attribute of variable \""<<varNames[i]<<"\" given but not pval2.";
               WARNING(os.str());
            } else {

               double pval2(0);
               br.GetData(pval2);

               //Set the prior
               setPrior( varNames[i], stringToPrior(prior), pval1, pval2);
            }
         }
      }

   }
         
}


void Variables::addToDOM( xercesc::DOMNode *node, 
      const std::vector<std::string> &names, 
      const std::vector<std::string> &ids, 
      const std::string &prefix ) const
{
   //First make sure all the variables are available.
   for (const auto& name: names) {
      if (fvars.find(name) == fvars.end()) {
         throw(VariableNotFound("Variable \""+name+"\" does not exist;"));
      }
   }

   //names and ids have to have the same size
   if (names.size() != ids.size())
      throw(VariableError("Names and IDs must have the same size"));
   
   //Now we can start creating the variable elements
   xercesc::DOMDocument *doc = node->getOwnerDocument();

   //To convert values to string
   std::ostringstream oss;

   for (size_t i = 0; i < names.size(); ++i) {
      auto it = getIterator(names[i]);

      xercesc::DOMElement* varEl = doc->createElement(utl::XStr("variable").unicodeForm());
      node->appendChild(varEl);
      varEl->setAttribute(utl::XStr("id").unicodeForm(), utl::XStr(ids[i]).unicodeForm());

      xercesc::DOMElement* valueEl = doc->createElement(utl::XStr("value").unicodeForm());
      varEl->appendChild(valueEl);
      oss << it->second.value;
      xercesc::DOMText* valueT = doc->createTextNode(utl::XStr(oss.str()).unicodeForm());
      oss.str("");
      valueEl->appendChild(valueT);

      //step
      if (it->second.type != FIXED) {
         xercesc::DOMElement* errorEl = doc->createElement(utl::XStr("step").unicodeForm());
         varEl->appendChild(errorEl);
         oss << it->second.error;
         xercesc::DOMText* errorT = doc->createTextNode(utl::XStr(oss.str()).unicodeForm());
         oss.str("");
         errorEl->appendChild(errorT);
      }

      //min
      if (it->second.type == LOWBND || it->second.type == BOTHBND) {
         xercesc::DOMElement* lBoundEl = doc->createElement(utl::XStr("min").unicodeForm());
         varEl->appendChild(lBoundEl);
         oss << it->second.lBound;
         xercesc::DOMText* lBoundT = doc->createTextNode(utl::XStr(oss.str()).unicodeForm());
         oss.str("");
         lBoundEl->appendChild(lBoundT);
      }

      //max
      if (it->second.type == UPPBND || it->second.type == BOTHBND) {
         xercesc::DOMElement* uBoundEl = doc->createElement(utl::XStr("max").unicodeForm());
         varEl->appendChild(uBoundEl);
         oss << it->second.uBound;
         xercesc::DOMText* uBoundT = doc->createTextNode(utl::XStr(oss.str()).unicodeForm());
         oss.str("");
         uBoundEl->appendChild(uBoundT);
      }

      //prior
      if (it->second.prior != NONE) {
         xercesc::DOMElement* priorEl = doc->createElement(utl::XStr("prior").unicodeForm());
         varEl->appendChild(priorEl);
         xercesc::DOMText* priorT = doc->createTextNode(utl::XStr(priorToString(it->second.prior)).unicodeForm());
         priorEl->appendChild(priorT);

         //pval1
         xercesc::DOMElement* pval1El = doc->createElement(utl::XStr("pval1").unicodeForm());
         varEl->appendChild(pval1El);
         oss << it->second.priorPar1;
         xercesc::DOMText* pval1T = doc->createTextNode(utl::XStr(oss.str()).unicodeForm());
         oss.str("");
         pval1El->appendChild(pval1T);

         //pval2
         xercesc::DOMElement* pval2El = doc->createElement(utl::XStr("pval2").unicodeForm());
         varEl->appendChild(pval2El);
         oss << it->second.priorPar2;
         xercesc::DOMText* pval2T = doc->createTextNode(utl::XStr(oss.str()).unicodeForm());
         oss.str("");
         pval2El->appendChild(pval2T);
      }

      //name
      if (names[i] != prefix+"_"+ids[i]) {
         xercesc::DOMElement* priorEl = doc->createElement(utl::XStr("name").unicodeForm());
         varEl->appendChild(priorEl);
         xercesc::DOMText* priorT = doc->createTextNode(utl::XStr(names[i]).unicodeForm());
         priorEl->appendChild(priorT);
      }
      
   }

}



//Adding variables from another variables object
void Variables::add(const Variables & vars) {
#pragma omp critical(Variables)
   {
      //Use the insert operator of the map to insert all the element of vars
      fvars.insert(vars.fvars.begin(), vars.fvars.end());

      //Find the maximum length
      fmaxLength = std::max(fmaxLength,vars.fmaxLength);
      if (fmaxLength > 500)
         std::cerr<<"Variable name is too long: "<<std::endl;
   }
}

void Variables::add(const std::string &name, double value){
   //Check for existance of the variable and throw an error if it already exists
#pragma omp critical(Variables)
   {
      if (fvars.find(name) != fvars.end()){

         throw(VariableAlreadyCreated("Variable \""+name+"\" already exists; cannot add a variable with the same name"));

      } else {

         //Create the variable struct 
         variable var;
         var.value = value;
         var.error = 0;
         var.uError = 0;
         var.lError = 0;
         var.upLim = 0;
         var.uSet = false;
         var.lSet = false;
         var.lBound = std::numeric_limits<double>::min();
         var.uBound = std::numeric_limits<double>::max();
         var.type = FIXED;

         //Prior is by default NONE
         var.prior = NONE;
         var.priorPar1 = 0;
         var.priorPar2 = 0;

         fvars[name] = var;

         fmaxLength = std::max(fmaxLength,int(name.size()));
         if (int(name.size()) > 500)
            std::cerr<<"Variable name is too long: "<<name<<std::endl;
         if (fmaxLength > 500)
            std::cerr<<"Variable name is too long: "<<name<<std::endl;

      }
   }
}

//Adding a variable to the ensemble
void Variables::add(const std::string &name, double value, double error){
   //Check for existance of the variable and throw an error if it already exists
#pragma omp critical(Variables)
   {
      if (fvars.find(name) != fvars.end()){

         throw(VariableAlreadyCreated("Variable \""+name+"\" already exists; cannot add a variable with the same name"));
      
      }else{

         //Create the variable struct 
         variable var;
         var.value = value;

         //Make sure the error is strictly positive.  Default to 1 percent error if not fulfilled or 0.01 if val is 0
         if ( error > 0 ) {
            var.error = error;
         } else {
            if (var.value == 0)
               var.error = 0.01;
            else
               var.error = 0.01*fabs(var.value);
         }
         var.uError = 0;
         var.lError = 0;
         var.upLim = 0;         
         var.uSet = false;
         var.lSet = false;
         var.lBound = std::numeric_limits<double>::min();
         var.uBound = std::numeric_limits<double>::max();
         var.type = FREE;

         //Use no prior as the default
         var.prior = NONE;
         var.priorPar1 = 0;
         var.priorPar2 = 0;

         fvars[name] = var;

         fmaxLength = std::max(fmaxLength,int(name.size()));
         if (int(name.size()) > 500)
            std::cerr<<"Variable name is too long: "<<name<<std::endl;
         if (fmaxLength > 500)
            std::cerr<<"Variable name is too long: "<<name<<std::endl;

      }
   }
}

void Variables::add(const std::string &name, double value, double error, double lowerBound, double upperBound){
   //Check for existance of the variable and throw an error if it already exists
#pragma omp critical(Variables)
   {
      if (fvars.find(name) != fvars.end()){

         throw(VariableAlreadyCreated("Variable \""+name+"\" already exists; cannot add a variable with the same name"));
      
      }else{
      
         //Create the variable struct
         variable var;
         var.value = value;
         //Make sure the error is strictly positive.  Default to 1 percent error if not fulfilled or 0.01 if val is 0
         if ( error > 0 ) {
            var.error = error;
         } else {
            if (var.value == 0)
               var.error = 0.01;
            else
               var.error = 0.01*fabs(var.value);
         }
         var.uError = 0;
         var.lError = 0;
         var.upLim = 0;
         var.uSet = true;
         var.lSet = true;
         var.uBound = upperBound;
         var.lBound = lowerBound;
         var.type = BOTHBND;

         //Set the prior, in this case the default is uniform
         var.prior = UNIFORM;
         var.priorPar1 = lowerBound;
         var.priorPar2 = upperBound;

         fvars[name] = var;

         fmaxLength = std::max(fmaxLength,int(name.size()));
         if (int(name.size()) > 500)
            std::cerr<<"Variable name is too long: "<<name<<std::endl;
         if (fmaxLength > 500)
            std::cerr<<"Variable name is too long: "<<name<<std::endl;

      }
   }
}


void Variables::setFromVariables(const Variables &vars) {

   // Loop over all the names in the current object and set the values from
   // the given object
   std::vector<std::string> varNames = getNames();

   for (size_t i(0); i < varNames.size(); ++i) {

      variable &var = getIterator(varNames[i])->second;
      const variable &newvar = vars.getIterator(varNames[i])->second;

      var = newvar;

   }

}

void Variables::removeUpperBound(const std::string &name) {
   //Get the variable
   variable &var = getIterator(name)->second;

   //Nothing to do if bound was not set
   if ( var.uSet ) {

      //Remove the upper bound
      var.uSet = false;

      //Fix the type and priors
      switch (var.type) {
         case BOTHBND:
            var.type = LOWBND;

            if ( var.prior == UNIFORM && var.priorPar1 == var.lBound && var.priorPar2 == var.uBound ) {
               var.prior = NONE;
               var.priorPar1 = 0;
               var.priorPar2 = 0;
            }
            break;

         case UPPBND:
            var.type = FREE;
            break;

         case LOWBND:
         case FREE:
         case FIXED:
         default:
            //These should not occur
            throw(VariableError("Variable type is wrong when removing upper bound.  There is a bug in the code."));
            break;

      }
   }

}

void Variables::removeLowerBound(const std::string &name) {
   //Get the variable
   variable &var = getIterator(name)->second;

   //Nothing to do if bound was not set
   if ( var.lSet ) {

      //Remove the upper bound
      var.lSet = false;

      //Fix the type and priors
      switch (var.type) {
         case BOTHBND:
            var.type = UPPBND;

            if ( var.prior == UNIFORM && var.priorPar1 == var.lBound && var.priorPar2 == var.uBound ) {
               var.prior = NONE;
               var.priorPar1 = 0;
               var.priorPar2 = 0;
            }
            break;

         case LOWBND:
            var.type = FREE;
            break;

         case UPPBND:
         case FREE:
         case FIXED:
         default:
            //These should not occur
            throw(VariableError("Variable type is wrong when removing lower bound for variable \"" + name + "\".  There is a bug in the code."));
            break;

      }
   }

}

//Set the variable bounds
void Variables::setUpperBound(const std::string &name, double bound) {

   // Get the variable
   variable &var = getIterator(name)->second;

   // Throw an error if trying to set an upper bound for a fixed variable
   if (var.type == FIXED)
      throw(VariableError("Cannot set an upper bound for the fixed variable \"" + name + "\".  Set the error first to free it"));

   var.uSet = true;
   var.uBound = bound;

   // Fix the variable type and prior if needed
   if ( var.lSet ) {

      var.type = BOTHBND;

      if ( var.prior == NONE ) {
         var.prior = UNIFORM;
         var.priorPar1 = var.lBound;
         var.priorPar2 = var.uBound;
      }

   } else {

      var.type = UPPBND;

   }
}

void Variables::setLowerBound(const std::string &name, double bound){

   // Get the variable
   variable &var = getIterator(name)->second;

   // Throw an error if trying to set an upper bound for a fixed variable
   if (var.type == FIXED)
      throw(VariableError("Cannot set a lower bound for the fixed variable \"" + name + "\".  Set the error first to free it"));

   var.lSet = true;
   var.lBound = bound;

   // Fix the variable type and prior if needed
   if ( var.uSet ) {

      var.type = BOTHBND;

      if ( var.prior == NONE ) {
         var.prior = UNIFORM;
         var.priorPar1 = var.lBound;
         var.priorPar2 = var.uBound;
      }

   } else {

      var.type = LOWBND;

   }
}

//Get the values of the variables as an array
std::vector<double> Variables::getValueVector() const{
   std::vector<double> values;

   //Loop the fvars map
   for (cvarIt it=fvars.begin(); it != fvars.end(); ++it){
      values.push_back((*it).second.value);
   }

   return values;
}

//Set the values of the variables from an array
void Variables::setValueVector(const std::vector<double> &values){

   //Check the sizes of the vector and fvars
   if (values.size() != fvars.size()){

      throw(VariableError("The size of the value array is not the same as the number of variables"));
   
   } else {

      //Loop the array and set the values
      varIt vit = fvars.begin();
      for (std::vector<double>::const_iterator it = values.begin(); it != values.end(); ++it, ++vit){
         (*vit).second.value = *it;
      }

   }
}

//Return the type of the variable
Variables::varType Variables::type(const std::string &name) const{
   const variable & var = getIterator(name)->second;

   return var.type;
}

//Return the value of the variable
double & Variables::value(const std::string &name){
   variable & var = getIterator(name)->second;

   return var.value;
}

const double & Variables::value(const std::string &name) const{
   const variable & var = getIterator(name)->second;

   return var.value;
}

double & Variables::operator [] (const std::string &name){
   return value(name);
}
const double & Variables::operator [] (const std::string &name) const{
   return value(name);
}

//Return the error of the variable
double Variables::error(const std::string &name) const {
   const variable & var = getIterator(name)->second;

   return var.error;
}

//Return asymetric errors.
std::pair<double,double> Variables::errors(const std::string &name) const {
   const variable & var = getIterator(name)->second;

   return std::pair<double,double>(var.lError,var.uError);
}

//Return upper limit
double Variables::upperLimit(const std::string &name) const {
   const variable & var = getIterator(name)->second;
   return var.upLim;
}

//Set the error of a variable
void Variables::setError( const std::string &name, double error ) {
   variable & var = getIterator(name)->second;

   //Make sure the error is strictly positive.  Default to 1 percent error if not fulfilled or 0.01 if val is 0
   if ( error > 0 ) {
      var.error = error;
   } else {
      if (var.value == 0)
         var.error = 0.01;
      else
         var.error = 0.01*fabs(var.value);
   }

   //Fix the type if needed
   if ( var.type == FIXED ) {
      var.type = FREE;

   }

}

//Set assymetric errors of a variable
void Variables::setErrors( const std::string &name, std::pair<double,double> errors) {
   variable & var = getIterator(name)->second;
   var.lError = errors.first;
   var.uError = errors.second;
}

//Set upper limit of a variable
void Variables::setUpperLimit( const std::string &name, double uplim) {
   variable & var = getIterator(name)->second;
   var.upLim = uplim;
}

void Variables::setPrior(const std::string &name, PRIOR prior, double par1, double par2) {
   variable & var = getIterator(name)->second;

   var.prior = prior;
   var.priorPar1 = par1;
   var.priorPar2 = par2;
}

Variables::PRIOR Variables::prior (const std::string & name) const {
   const variable & var = getIterator(name)->second;

   return var.prior;
}

double Variables::priorPar1 (const std::string & name) const {
   const variable & var = getIterator(name)->second;

   return var.priorPar1;
}

double Variables::priorPar2 (const std::string & name) const {
   const variable & var = getIterator(name)->second;

   return var.priorPar2;
}


//Return the names of the variables in a vector
std::vector<std::string> Variables::getNames() const{

   std::vector<std::string> names;

   for (cvarIt it=fvars.begin(); it != fvars.end(); ++it){
      names.push_back((*it).first);
   }

   return names;
}

//Return the upper and lower bounds
double  Variables::upperBound(const std::string &name) const{
   const variable &var = getIterator(name)->second;

   if (! var.uSet){

      throw(BoundNotSet("Upper bound not set for variable \""+name+"\""));
   
   }else{
      
      return var.uBound;
   
   }
}

double Variables::lowerBound(const std::string &name) const{
   const variable &var = getIterator(name)->second;

   if (! var.lSet){

      throw(BoundNotSet("Lower bound not set for variable \""+name+"\""));
   
   }else{
      
      return var.lBound;
   
   }
}

void Variables::fix(const std::string &name) {
   variable &var = getIterator(name)->second;

   var.type = FIXED;

}

//Output operator
std::ostream & operator << (std::ostream &os, const Variables &vars){

   //Print the header
   int width = std::min(500, vars.fmaxLength);
   os << std::setw(width) << "Variable Name";
   os << " |     Value+High-Low (Error) Upper Limit     | Lower Bound | Upper Bound |   Prior    | Parameter 1 | Parameter 2 |\n";

   //Loop the variables
   for (Variables::cvarIt it = vars.fvars.begin(); it != vars.fvars.end(); ++it){
      
      //The name
      os << std::left << std::setw(width) << (*it).first << std::setw(3) << " | "<<std::right;
      
      //The value and errors
      os << std::setw(9) << std::setprecision(5) << (*it).second.value;
      if ((*it).second.type != Variables::FIXED) {
        if ((*it).second.lError!=0 || (*it).second.uError!=0) {
          //os <<"+"<<(*it).second.uError/(*it).second.error<<" -"<<(*it).second.lError/(*it).second.error<<std::setw(2);
          os <<"+"<<std::setw(7)<<std::setprecision(5)<<(*it).second.uError<<"-"<<std::setw(7)<<std::setprecision(5)<<(*it).second.lError;
        } else {
           os << std::setw(16)<<" ";
        }
        os <<"("<<std::setw(7)<<std::setprecision(5)<<(*it).second.error<<")";
        if((*it).second.upLim!=0.){
          os << std::setw(7) << std::setprecision(5)<<(*it).second.upLim<<" | ";
        } else {
          os << std::setw(7) <<" "<<" | ";
        }
      } else {
         os << std::setw(41) <<" "<<" | ";
      }
      
      //Check for the lower bound and print
      if ((*it).second.lSet){
         os << std::setw(11) << (*it).second.lBound << std::setw(3) << " | ";
      }else{
         os << std::setw(11) << "N/A" << std::setw(3) << " | ";
      }
      
      //Check for the upper bound and print
      if ((*it).second.uSet){
         os << std::setw(11) << (*it).second.uBound  << std::setw(3) << " | ";
      }else{
         os << std::setw(11) << "N/A"  << std::setw(3) << " | ";
      }

      //Check for prior
      if ((*it).second.prior == Variables::NONE) {
         os << std::setw(11) << "N/A" << std::setw(3) << " | ";
         os << std::setw(11) << "N/A" << std::setw(3) << " | ";
         os << std::setw(11) << "N/A" << std::setw(1) << "\n";
      } else {
         os << std::setw(11) << Variables::priorToString((*it).second.prior) << std::setw(3) << " | ";
         os << std::setw(11) << (*it).second.priorPar1 << std::setw(3) << " | ";
         os << std::setw(11) << (*it).second.priorPar2 << std::setw(1) << "\n";
      }
   }
   return os;
}

void Variables::clear() {
   fvars.clear();
   fmaxLength = 0;
}

Variables::PRIOR Variables::stringToPrior ( const std::string &priorString ) {

   // Convert to upper case
   std::string upper(priorString);
   // explicit cast needed to resolve ambiguity
   std::transform(upper.begin(), upper.end(), upper.begin(), (int(*)(int)) std::toupper);

   if ( upper == "NONE" ) {
      return NONE;
   } else if ( upper == "UNIFORM" ) {
      return UNIFORM;
   } else if ( upper == "LOGUNIFORM" ) {
      return LOGUNIFORM;
   } else if ( upper == "CAUCHY" ) {
      return CAUCHY;
   } else if ( upper == "GAUSSIAN" ) {
      return GAUSSIAN;
   } else if ( upper == "LOGNORMAL" ) {
      return LOGNORMAL;
   }

   throw (VariableError("Prior type \"" + upper + "\" is not known"));
   
}

std::string Variables::priorToString(Variables::PRIOR prior) {

   switch ( prior ) {
      case NONE:
         return "NONE";
      case UNIFORM:
         return "UNIFORM";
      case LOGUNIFORM:
         return "LOGUNIFORM";
      case GAUSSIAN:
         return "GAUSSIAN";
      case CAUCHY:
         return "CAUCHY";
      case LOGNORMAL:
         return "LOGNORMAL";
      default:
         throw (VariableError("Prior is not known"));
   }
}

Variables::varIt Variables::getIterator(const std::string &name) {
   //Find the variable and return an error if not found
   varIt it = fvars.find(name);
   if (it == fvars.end()){

      throw(VariableNotFound("Variable \""+name+"\" does not exist;"));

   }

   return it;
}
Variables::cvarIt Variables::getIterator(const std::string &name) const {
   //Find the variable and return an error if not found
   cvarIt it = fvars.find(name);
   if (it == fvars.end()){

      throw(VariableNotFound("Variable \""+name+"\" does not exist;"));

   }

   return it;
}
}
