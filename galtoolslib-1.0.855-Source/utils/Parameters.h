#ifndef _utl_Parameters_h
#define _utl_Parameters_h

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <exception>
#include <sstream>
#ifdef USE_NATIVE_REGEX
#include <regex>
#endif

namespace utl {

  /** \brief Class to handle parsing of simple parameter = value pairs.
   *
   * The parameter = value string is split into name and value pair on the first = occurance.
   *
   * The parameter value is parsed into variables using the >> operator.
   *
   * It is possible to create vectors or sets with either comma or white space separated values.
   * \note All commas will be removed from strings when parsing to a vector.
   *
   * Vectors can also be created with the notation
   * name%index = value
   * where name is the name of the parameter and index is integer representing the index.
   * The class first reads in the values given in the 
   * name = value0, value1, value2
   * and then looks for name%index kind.  
   *
   */
  class Parameters {
  private:
    typedef std::map<std::string, std::string> maptype; //!< Alias typedef for a map of string pairs, used to store the parameter name and value
    typedef maptype::const_iterator cmapiterator; //!< constant interator for the map
    maptype fParameters; //!< Storage for the parmeter name/value pairs
    std::string fCommentString; //!< The comment string is defined in the constructor.
    mutable std::set<std::string> faccessedParameters; //!< Keep track of parameters that have been accessed.
    
  public:
    /** \brief Error class for the Parameters class. */
    class ParameterError : public std::exception {
       private:
          std::string eMessage; //!< store the error message
       public:
          /** \brief Construct an error with the default message "Error in
           * Parmaeters class"
           */
          ParameterError() : eMessage("Error in Parameters class"){}
          /** \brief Construct an error with a customized error message.
           *
           * \param message is the error message to display.
           */
          ParameterError(const std::string & message) : eMessage(message){}
          ~ParameterError() throw(){}
          /** \brief Return the error message of this error */
          virtual const char* what () const throw(){
             return(eMessage.c_str());
          }
    };

  /** \brief Construct an empty Parameters object.
   *
   * \param commentString is the comment string; all characters on a line, including and after this string
   * in the parsed input streams will be discarded.  Defaults to "#".
   */
  Parameters(std::string commentString="#") : fCommentString(commentString) {}

  /** \brief Construct a Parameters object from an input stream.
   *
   * \param is is the input stream to be parsed.  It will be parsed line by
   * line and all comments stripped off.
   * \param commentString is the comment string; all characters on a line, including and after this string
   * in the parsed input streams will be discarded.  Defaults to "#".
   */
  Parameters(std::istream & is, std::string commentString="#") : fCommentString(commentString) {ParseIstream(is);}
  
  /** \brief Parse a parameter stream into the Parameters object.
   *
   * \param is is the input stream to be parsed.  It will be parsed line by
   * line and all comments stripped off.
   */
  void ParseIstream(std::istream & is);
  
  /** \brief Print all the parameters to a output stream
   *
   * \param os is the output stream to which the parameters will be printed.
   */
  void PrintOstream(std::ostream & os) const;
  
  /** \brief Return a parameter value given its name.
   *
   * This is a templated function and the results depend on the type of the
   * input parameter.
   * \param parameter is the name of the parameter
   * \param value is a reference to a variable which will contain the
   * parameter value on return.
   *
   * The basic input mechanism of c++ is used to parse the value string, so
   * only the first non-whitespace containing string will be parsed for the
   * value.
   */
  template<typename T>
     void getParameter(const std::string & parameter, T & value) const{
        const cmapiterator it = fParameters.find(parameter);
        if (it == fParameters.end()) throw(ParameterError("Parameter \""+parameter+"\" not found"));
        //Create a stream to read the parameter
        std::istringstream iss((*it).second);
        iss >> value;
        if (iss.fail() || iss.bad()) throw(ParameterError("Error reading parameter \""+parameter+"\""));
        faccessedParameters.insert(parameter);
     }
  
  /** \brief Return a parameter value as a string given its name.
   *
   * \param parameter is the name of the parameter
   * \param value is a reference to a string which will contain the
   * parameter value on return.
   *
   * The unparsed value string is returned.
   */
  void getParameter(const std::string & parameter, std::string & value, bool FullString) const {
    if (FullString) {
      const cmapiterator it = fParameters.find(parameter);
      if (it == fParameters.end()) throw(ParameterError("Parameter \""+parameter+"\" not found"));
      value = (*it).second;
      faccessedParameters.insert(parameter);
    } else {
      getParameter(parameter, value);
    }
  }
  
  /** \brief Return a list of parameters that have not been accessed
   */
  std::set<std::string> getUnusedParameters() const;
  
  /** \brief Return a parameter value as a vector given its name.
   *
   * This is a templated function and the results depend on the type of the
   * input parameter.
   * \param parameter is the name of the parameter
   * \param vec is a reference to a vector which will contain the
   * parameter values on return. 
   *
   * If the name%index = value notation is used, the vector values will not be erased. 
   * The values for the given indices will be updated and the vector size increased
   * to accomodate them.  No initialization is done for unsigned values.
   * If the usual name = value1 value2 ... is used the vector is cleared.
   *
   * The basic input mechanism of c++ is used to parse the value string, so
   * the value string will be parsed to type T for each white space seperated
   * object in the value string.
   *
   * All commas in the string will be replaced with a whitespace before parsing.
   */
  template<typename T>
    void getParameter(const std::string & parameter, std::vector<T> & vec) const{
    auto it = fParameters.find(parameter);
    bool found(false);
    if (it != fParameters.end()) {

       found = true;
       //Clear the vector
       vec.clear();

       //Replace commas with space
       std::string tempstr = (*it).second;
       ReplaceAll(tempstr, ",", " ");

       //Trim again to remove space left by comma at front or end of string
       TrimWhitespace(tempstr);

       //Create a stream to read the parameter
       std::istringstream iss(tempstr);
       while (iss.good()){
          T value;
          iss >> value;
          if (iss.fail() || iss.bad()) throw(ParameterError("Error reading parameter \""+parameter+"\""));
          vec.push_back(value);
       }
    } 

#ifdef USE_NATIVE_REGEX
    //Look for name%index parameters
    std::regex rex(parameter + "%([0-9]+)");
    std::smatch sm;
    for (it = fParameters.begin(); it != fParameters.end(); ++it) {

       if ( std::regex_match(it->first, sm, rex) ) {

          found = true;

          long indx = std::stol(sm.str(1));
          if (indx < 0)
             throw(ParameterError("Index < 0 when reading \"" + it->first +"\""));

          if (size_t(indx)+1 > vec.size()) 
             vec.resize(indx+1);

          std::istringstream iss((*it).second);
          iss >> vec[indx];
          if (iss.fail() || iss.bad()) throw(ParameterError("Error reading parameter \""+it->first+"\""));
       }
    }
#endif
    if (!found) throw(ParameterError("Parameter \""+parameter+"\" not found"));
    faccessedParameters.insert(parameter);
  }
  
  /** \brief Return a parameter value as a vector given its name.
   *
   * This is a templated function and the results depend on the type of the
   * input parameter.
   * \param parameter is the name of the parameter
   * \param s is a reference to a set which will contain the
   * parameter values on return.
   *
   * The basic input mechanism of c++ is used to parse the value string, so
   * the value string will be parsed to type T for each white space seperated
   * object in the value string.
   *
   * All commas in the string will be replaced with a whitespace before parsing.
   */
  template<typename T>
    void getParameter(const std::string & parameter, std::set<T> & s) const{
    const cmapiterator it = fParameters.find(parameter);
    if (it == fParameters.end()) throw(ParameterError("Parameter \""+parameter+"\" not found"));
    //Clear the vector
    s.clear();
    //Replace commas with space
    std::string tempstr = (*it).second;
    ReplaceAll(tempstr, ",", " ");
    //Trim again to remove space left by comma at front or end of string
    TrimWhitespace(tempstr);
    //Create a stream to read the parameter
    std::istringstream iss(tempstr);
    while (iss.good()){
      T value;
      iss >> value;
      if (iss.fail() || iss.bad()) throw(ParameterError("Error reading parameter \""+parameter+"\""));
      s.insert(value);
    }
  }
  
  /** \brief Set the value of a parameter given its name.
   *
   * This is a templated function and the results depend on the type of the
   * input parameter.
   * \param parameter is the name of the parameter
   * \param value is a reference to a variable which holds the parameter
   * value
   *
   * The basic output mechanism of c++ is used to convert the value to string.
   */
  template<typename T>
    void setParameter(const std::string & parameter, const T & value){
    std::ostringstream oss;
    oss << value;
    fParameters[parameter] = oss.str();
  }
  /** \brief Set the value of a parameter given its name.
   *
   * \param parameter is the name of the parameter
   * \param value is a reference to a string which holds the parameter
   * value
   *
   * The value of the parameter is set to that of the string.
   */
  void setParameter(const std::string & parameter, const std::string & value){
    fParameters[parameter] = value;
  }

  /** \brief Add the parameters from pars and override identically named parameters
   *
   * \param pars is a parameter object.
   */
  void setParameters(const Parameters &pars);
  
  
  /** \brief Strip comment string and everything after it from the input str.
   *
   * \param str is the string from which to strip comments.
   * \return a copy of the original string after comments have been removed.
   *
   * This utilises the comment string set in the constructor
   */
  static std::string & StripComments(std::string & str, const std::string &comment);
  
  /** \brief Trim white space from front and end of string */
  static std::string & TrimWhitespace(std::string & str);
  
  //! Replace all occurances of string with another string
  static std::string & ReplaceAll(std::string & str, const std::string & search, const std::string & replace);
  };
  
} //End namespace utl

#endif
