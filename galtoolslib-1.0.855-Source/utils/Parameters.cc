#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include "Parameters.h"

namespace utl {

std::string & Parameters::StripComments(std::string & str, const std::string & comment) {
	//Throws an out of range exception if no comment
	try{
		str.erase(str.find(comment), std::string::npos);
	}catch (std::out_of_range & exception){}
	return (str);
}

std::string & Parameters::TrimWhitespace(std::string & str) {
	//Throws an out of range exception if no whitespace at end
	try{
		str.erase(str.find_last_not_of(" \t\n\v\f\r")+1, std::string::npos);
	}catch (std::out_of_range & exception){}
	try{
		str.erase(0,str.find_first_not_of(" \t\n\v\f\r"));
	}catch (std::out_of_range & exception){}
	return (str);
}

std::string & Parameters::ReplaceAll(std::string & str, const std::string & search, const std::string & replace) {
	//Loop with find and replace until end of string is reached.
	size_t pos = str.find(search);

	while (pos != std::string::npos) {
		str.replace(pos, search.size(), replace);
		//In case replace contains search, we start the search after the inserted string
		pos = str.find(search, pos+replace.size());
	}

	return(str);
}

void Parameters::ParseIstream(std::istream & is){
	//At most one parameter per line so we read in one line at a time
	std::string buffer;

	while (is.good()){ 
		std::getline(is, buffer);

		//Erease comments and trim whitespace off the end
		StripComments(buffer, fCommentString);
		TrimWhitespace(buffer);

		//Nothing to do for empty lines
		if (buffer.empty()) continue;

		//Find the = sign to split the string in 2
		size_t pos = buffer.find("=");
		if (pos == std::string::npos){
			std::cerr<<"Failed parsing line:\n"<<buffer<<"\n";
			continue;
		}

		std::string key = buffer.substr(0,pos);
		TrimWhitespace(key);

		std::string value = buffer.substr(pos+1);
		TrimWhitespace(value);

		//Insert the parameter
		fParameters[key] = value;
	}
}

void Parameters::PrintOstream(std::ostream & os) const {
	os<<"#Listing all parameters: "<<std::endl;
	for ( cmapiterator it = fParameters.begin();it != fParameters.end(); ++it) {
		os<<it->first<<"="<<it->second<<std::endl;
	}
}

void Parameters::setParameters(const Parameters &pars) {

   for (auto it = pars.fParameters.begin(); it != pars.fParameters.end(); ++it)
      fParameters[it->first] = it->second;

}


std::set<std::string> Parameters::getUnusedParameters() const {

   std::set<std::string> unusedParameters;

#ifdef USE_NATIVE_REGEX
   //Need to account for parameters using the name%index convention
   std::regex rex("(.*)%([0-9]+)");
   std::smatch sm;
#endif

   for (const auto &pair : fParameters) {

      std::string parameter = pair.first;

#ifdef USE_NATIVE_REGEX
      //Look for name%index parameters
      if ( std::regex_match(parameter, sm, rex) ) {
         parameter = sm.str(1);
      }
#endif

      if ( faccessedParameters.find(parameter) == faccessedParameters.end() )
         unusedParameters.insert(parameter);
   }

   return unusedParameters;
}

}
