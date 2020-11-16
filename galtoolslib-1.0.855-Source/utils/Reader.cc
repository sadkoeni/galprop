// $Id: Reader.cc 6228 2007-08-14 19:40:36Z darko $

#include <Reader.h>
#include <ErrorLogger.h>
#include <sstream>

#include <CLHEP/Evaluator/Evaluator.h>
#include <Units.h>

#include <ReaderErrorReporter.h>
#include <xercesc/framework/MemBufInputSource.hpp>

using namespace utl;
using namespace xercesc;
using namespace std;

namespace utl {

  // Units definition to set the units system in expression evaluator
  const double r_length             = 1.0;                    // meter          
  const double r_mass               = 1.0 / (e_SI * 1.0e-18); // eV/c**2        
  const double r_time               = 1.0e+9;                 // nanosecond     
  const double r_current            = 1.0 / (e_SI / 1.0e-9);  // e/nanosecond   
  const double r_temperature        = 1.0;                    // Kelvin         
  const double r_amount             = 1.0;                    // mole           
  const double r_luminous_intensity = 1.0;                    // candela        

  HepTool::Evaluator gEval;

}

    
//! Old constructor with older validation options.
/*!
  Preserved to keep old code from breaking. 
  \param name string: name of the XML file
  \param b validation option
*/
Reader::Reader(const string& name, const bool v) :
  //fXMLFileName,
  //fMemInput,
  //fParser(0),
  //fErrReporter(0),
  //fCardname,
  //fDataCardBranch,
  fDocumentNode(0),
  fValidationMode(Reader::eNONE)
  //fTopBranch,
  //fSpaces
{
  fXMLFileName = name;
  if (v)
    SetUp(Reader::eDTD);
}


//! Constructor with validation options.
/*!
  \param name of the XML file
  \param v validation option
*/
Reader::Reader(const string& name, const Validation v) :
  //fXMLFileName,
  //fMemInput,
  //fParser(0),
  //fErrReporter(0),
  //fCardname,
  //fDataCardBranch,
  fDocumentNode(0)
  //fValidationMode(Reader::eNONE)
  //fTopBranch,
  //fSpaces
{
  fXMLFileName = name;  
  SetUp(v);
}


//! Constructor from an input string (as opposed to file). Validation not (yet) used.
/*!
  \param is input XML string.
*/
Reader::Reader(const ReaderStringInput& is, const Validation v) :
  //fXMLFileName,
  //fMemInput,
  fParser(0)
  //fErrReporter(0),
  //fCardname,
  //fDataCardBranch,
  //fDocumentNode(0),
  //fValidationMode(Reader::eNONE)
  //fTopBranch,
  //fSpaces
{
  fMemInput = is;
  SetUp(v);
}


//! Assign all the variables to the Reader
/*!
  \param v validation option
*/
void
Reader::SetUp(const Validation v)
{
  fValidationMode = v;
  Initialize();  // do not move Initialize() away from the beginning!! other DOM stuff needs it. 
  // Set up the evaluator and devine additional units
  gEval.setStdMath();
  gEval.setSystemOfUnits(r_length,r_mass,r_time,r_current,
			 r_temperature,r_amount,r_luminous_intensity);
  gEval.setVariable("EeV", "1e18*eV");

  Parse();
}


void
Reader::Initialize()
{
  XMLPlatformUtils::Initialize();
}


//! Parse the file 
void
Reader::Parse()
{
  fParser = new XercesDOMParser;

  fParser->setDoSchema(false);
  fParser->setValidationScheme(XercesDOMParser::Val_Never);
  fParser->setCreateEntityReferenceNodes(false); // MAY WANT TO CHANGE THIS AT SOME POINT
  
  switch (fValidationMode) {
  case Reader::eDTD:
    fParser->setValidationScheme(XercesDOMParser::Val_Auto);
    break;
  case Reader::eSCHEMA:
    fParser->setDoSchema(true);
    fParser->setValidationScheme(XercesDOMParser::Val_Auto);
    fParser->setValidationSchemaFullChecking(true);	
    fParser->setDoNamespaces(true);
    break;
  case Reader::eNONE:
    fParser->setValidationScheme(XercesDOMParser::Val_Never);
    break;
  default:
    ERROR("You have selected an invalid validation mode");
    exit(EXIT_FAILURE);
    break;
  }
  
  fErrReporter = (ErrorHandler*) new ReaderErrorReporter();
  fParser->setErrorHandler(fErrReporter.Get());
  
  if (fMemInput.GetInputString() != "") {
    // request was made to read directly from input (not from file)

    MemBufInputSource memSource((const XMLByte*)fMemInput.GetInputString().c_str(),
				fMemInput.GetInputString().size(), "MemoryXML", false);
    try {
      fParser->parse(memSource);
    } catch (const XMLException& e) {
      ostringstream msg;
      const StrX message(e.getMessage());
      msg << "An errror occurred during parsing, message: "
          << message;
      ERROR(msg);
      throw std::runtime_error(msg.str());
    }
  } else try {
    fParser->parse(fXMLFileName.c_str());
  } catch (const XMLException& e) {
    ostringstream msg;
    const StrX message(e.getMessage());
    msg << "An errror occurred during parsing, message: "
        << message;	  	
    ERROR(msg);
    throw std::runtime_error(msg.str());
  }
  
  fDocumentNode = fParser->getDocument();   
  
  // initialize top branch to the document element.
  fTopBranch.fDomNode = fDocumentNode->getDocumentElement();  
}


//! Dump the tree of nodes
void
Reader::DumpTree()
{
  cout << "\nDumping information in XML file : " << fXMLFileName << "\n"
          "-----------------------------------------------------------" << endl;
  DumpMe(GetTopBranch());
}


/*!
 * This utility function takes a string in the standard XML Schema dateTime format
 *
 * \code
 * CCYY-MM-DDThh:mm:ss
 * \endcode
 *
 * with an optional time zone appended in the form
 *
 * \code
 * +/-hh:mm
 * \endcode
 *
 * If the time zone is omitted, or if an optional 'Z' follows the :ss
 * field, then the time is in UTC. For example:
 *
 * \code
 * 2001-04-09T06:00:00         (UTC time)
 * 2001-04-09T06:00:00Z        (UTC time)
 * 2001-04-09T06:00:00-05:00   (-5 hours UTC)
 * 2001-04-09T06:00:00+01:00   (+1 hour UTC)
 * \endcode
 *
 * This function returns separate integers representing year, month,
 * day, hour, minute, second and the offset from UTC.  Note that <b>no
 * checking is done</b> to guarantee the formatting of XMLTime adheres
 * to the Schema standard - this kind of checking is what Schema is
 * for!
 * 
 * For more details on the dataTime format, see, for example: Aaron
 * Skonnard, Martin Gudgin, "Essential XML Quick Reference"
*/
/*void
Reader::BreakdownXMLTime(const string& xmlTime,
                         unsigned int& year,
                         unsigned int& month,
                         unsigned int& day,
                         unsigned int& hour,
                         unsigned int& minute,
                         double& second,
                         int& utmOffsetHour,
                         unsigned int& utmOffsetMinute)
{
  istringstream is(xmlTime);
  char c;
  is >> year >> c >> month >> c >> day >> c
     >> hour >> c >> minute >> c >> second;
  c = '\0';
  is >> c;
  utmOffsetHour = utmOffsetMinute = 0;
  if (c == '-' || c == '+') {
    char cc;
    is >> utmOffsetHour >> cc >> utmOffsetMinute;
    utmOffsetHour *= (c == '+' ? 1 : -1);
  }
}


void
Reader::GetXMLTime(const string& str, TimeStamp& time)
{
  unsigned int theYear;
  unsigned int theMonth;
  unsigned int theDay;
  unsigned int theHour;
  unsigned int theMinute;
  double theSec;
  int offHour;
  unsigned int offMinute;
  Reader::BreakdownXMLTime(str, theYear, theMonth, theDay,
                           theHour, theMinute, theSec,
                           offHour, offMinute);
  const unsigned int theSecond = int(theSec);
  const unsigned int theNanosecond = int((theSec - theSecond)*1e9 + 0.5);
  time = TimeStamp(theYear, theMonth, theDay, theHour, theMinute, theSecond, theNanosecond);
  if (offHour || offMinute) {
    const int sign = offHour < 0 ? 1 : -1;
    time += TimeInterval(sign * (std::abs(offHour)*hour + offMinute*minute));
  }
}
*/

//! Dump information in the XML tree (recursive)
/*!
  \param b Branch to be dumped
*/
void
Reader::DumpMe(const Branch& b)
{
  const string additionalSpace = "   ";
  
  vector<string> bData; // put data in a string vector to avoid spaces & CR in printout
  b.GetData(bData);

  map<string, string> atts = b.GetAttributes();

  cout << fSpaces << b.GetBranchNameString() << " : ";
  if (!atts.empty())
    cout << "[";
  for (map<string, string>::iterator it = atts.begin();
       it != atts.end() ; ++it) {
    cout << it->first  << "=";
    cout << it->second << " ";
  }
  if (!atts.empty())
    cout << "] " ;
  for (vector<string>::iterator it = bData.begin();
       it != bData.end(); ++it)
    cout << *it << " ";
  cout << endl;
  fSpaces += additionalSpace;

  for (Branch cb = b.GetFirstChild(); cb; cb = cb.GetNextSibling())
    DumpMe(cb);

  const int len = fSpaces.length();
  fSpaces.erase(len - additionalSpace.length());
}


namespace utl {

  string Branch::fgWarning;

}


//! Get the first branch child
/*!
  \return First child branch
*/
Branch
Branch::GetFirstChild()
  const
{
  if (!fDomNode) {
    WARNING(fgWarning);
    ERROR("Getting first child in a null-Branch");
    exit(EXIT_FAILURE);
  }

  DOMNode* childNode = fDomNode->getFirstChild();
  while (childNode && childNode->getNodeType() != DOMNode::ELEMENT_NODE)
    childNode = childNode->getNextSibling();

  if (!childNode) {
    ostringstream warn;
    warn << "First child in branch '" << GetBranchNameString() << "' not found";
    fgWarning = warn.str();
  }
  return Branch(childNode);
}


//! Get branch child
/*!
  \param requestedName Name requested in the branch
  \param requestedAttributeMap Map of the attributes attached to the branch
  \return Child branch
*/
Branch
Branch::GetChild(const string& requestedName,
                 map<string, string> requestedAttributeMap)
  const
{
  if (!fDomNode) {
    WARNING(fgWarning);
    ostringstream error;
    error << "Getting child '" << requestedName << "' on a null-Branch";
    ERROR(error);
    exit(EXIT_FAILURE);
  }

  requestedAttributeMap.erase("unit");
  requestedAttributeMap.erase("UNIT");

  for (DOMNode* childNode = fDomNode->getFirstChild(); 
       childNode; childNode = childNode->getNextSibling()) {

    if (childNode->getNodeType() == DOMNode::ELEMENT_NODE) {

      const string foundName = StrX(childNode->getNodeName()).localForm();

      if (foundName == requestedName) {

	// Get all the attributes into a map.
	map<string, string> foundAttributeMap;
	const DOMNamedNodeMap* const attributes = childNode->getAttributes();

	for (unsigned int j = 0; j < attributes->getLength(); ++j) {

	  const DOMNode* const attribute = attributes->item(j);

	  foundAttributeMap.
	    insert(pair<string, string>(StrX(attribute->getNodeName()).localForm(),
                                        StrX(attribute->getNodeValue()).localForm()));

	} // filled up foundAttributeMap
	  
	// First, discard any unit attributes in either the found or requested maps.
	// Units are treated as a special case.
	foundAttributeMap.erase("unit");
	foundAttributeMap.erase("UNIT");
	  
	// Check for exact match between remaining attributes and attribute
	// values in the requested and found maps
	if (foundAttributeMap.size() == requestedAttributeMap.size()) {

	  bool foundIt = true;

	  // try to disqualify equality
	  for (map<string, string>::iterator foundIter = foundAttributeMap.begin();
	       foundIter != foundAttributeMap.end(); ++foundIter) {

	    map<string, string>::iterator requestedIter =
	      requestedAttributeMap.find(foundIter->first);
	    if (requestedIter == requestedAttributeMap.end() ||
		foundIter->second != requestedIter->second)
	      foundIt = false;
	  }

	  if (foundIt)
	    return Branch(childNode);
	}
      } // same name
    } // condition ELEMENT_NODE
  } // loop child nodes

  ostringstream warn;
  warn << "Child '" << requestedName << "' in branch '"
       << GetBranchNameString() << "' not found";
  fgWarning = warn.str();
  return Branch(); // null-branch
}


//! Get child from the branch
/*!
  \param requestedName Name requested in the branch
  \return Child branch
*/
Branch
Branch::GetChild(const string& requestedName)
  const
{
  map<string, string> dummy;
  return GetChild(requestedName, dummy);
}


//! Get child from the branch
Branch
Branch::GetChild(const string& requestedName, const string& id)
  const
{
  map<string, string> idMap;
  idMap["id"] = id;
  return GetChild(requestedName, idMap);
}


//! Get the map of the branch attibutes
map<string, string>
Branch::GetAttributes()
  const
{
  if (!fDomNode) {
    WARNING(fgWarning);
    ERROR("Getting attributes of a null-Branch.");
    exit(EXIT_FAILURE);
  }

  map<string, string> attMap;

  const DOMNamedNodeMap* const attributes = fDomNode->getAttributes();
  for (unsigned int j = 0; j < attributes->getLength(); ++j) {
    const DOMNode* const attribute = attributes->item(j);

    attMap[StrX(attribute->getNodeName()).localForm()] =
      StrX(attribute->getNodeValue()).localForm();
  }
  return attMap;
}


string
Branch::GetBranchNameString()
  const
{
  if (!fDomNode) {
    WARNING(fgWarning);
    ERROR("Getting the name of a null-Branch.");
    exit(EXIT_FAILURE);
  }
  
  return StrX(fDomNode->getNodeName()).localForm();    
}


Branch
Branch::GetNextSibling()
  const
{
  if (!fDomNode) {
    WARNING(fgWarning);
    ERROR("Getting next sibling in a null-Branch");
    exit(EXIT_FAILURE);
  }

  DOMNode* siblingNode = fDomNode->getNextSibling();   
    
  while (siblingNode && siblingNode->getNodeType() != DOMNode::ELEMENT_NODE)
    siblingNode = siblingNode->getNextSibling();

  if (!siblingNode) {
    ostringstream warn;
    warn << "Next sibling of branch '" << GetBranchNameString() << "' not found";
    fgWarning = warn.str();
  }
  return Branch(siblingNode);
}


Branch
Branch::GetSibling(const string& requestedName, map<string, string>& attributeMap)
  const
{
  // check that current node is valid.
  if (!fDomNode) {
    WARNING(fgWarning);
    ERROR("Getting sibling of a null-Branch");
    exit(EXIT_FAILURE);
  }

  // back up to the parent.
  DOMNode* const parentNode = fDomNode->getParentNode();

  if (!parentNode) {
    ostringstream warn;
    warn << "Parent of branch '" << GetBranchNameString() << "' not found";
    fgWarning = warn.str();
  }

  Branch parentBranch(parentNode);
    
  return parentBranch.GetChild(requestedName, attributeMap);
}


Branch
Branch::GetSibling(const string& requestedName)
  const
{
  map<string, string> dummy;
  return GetSibling(requestedName, dummy);
}


//! Get sibling
/*!
  \param requestedName Requested name
  \param id Attribute name
  \return Sibling branch
*/
Branch
Branch::GetSibling(const string& requestedName, const string& id)
  const
{
  map<string, string> idMap;
  idMap["id"] = id;
  return GetSibling(requestedName, idMap);
}


// -------------------- GetData methods ------------------------------


//! Get Data string
/*!
  \return String with data attached
*/
string
Branch::GetDataString()
  const
{
  if (!fDomNode) {
    WARNING(fgWarning);
    ERROR("Getting data from a null-Branch");
    exit(EXIT_FAILURE);
  }

  // Find the data residing beneath this tag
  const DOMNodeList* const dataNodes = fDomNode->getChildNodes();

  if (dataNodes) {

    string dataString;

    for (unsigned int k = 0; k < dataNodes->getLength(); ++k) {
 
      const DOMNode* const currentNode = dataNodes->item(k);
	
      // For TEXT nodes, postpend the node onto the dataString. NB if
      // there are multiple TEXT nodes beneath this element, and they
      // are, for example, separated by one or more sub- elements,
      // these text nodes will be concatenated into a single data
      // string.
      if (currentNode->getNodeType() == DOMNode::TEXT_NODE)
	dataString += StrX(currentNode->getNodeValue()).localForm();
    }

    // Remove any leading and/or trailing whitespace (and CR's) from the string

    const string& dStr = dataString;
    const string::const_iterator end = dStr.end();
    string::const_iterator blackStart = dStr.begin();
    string::const_iterator whiteStart;

    // find last pair of starting black and white text
    do {
      whiteStart = find_if(blackStart, end, IsSpace());
      blackStart = find_if(whiteStart, end, IsNotSpace());
    } while (blackStart != end);

    // find starting black text at the beginning

    blackStart = find_if(dStr.begin(), end, IsNotSpace());
    string trimmedString;

    if (blackStart != dataString.end()) {  // case of empty element
      try {
	trimmedString.assign(blackStart, whiteStart);
      } catch (...) {
	ostringstream warn;
	warn << "String trimming failure. string = '"
	     << dataString << "'. "
	        "(string will be treated as empty!) ";
	WARNING(warn);
      }
    }

    return trimmedString;
  }

  return string();
}


//! Get the unit of the token
/*!
  \return Unit of the token
*/
double
Branch::GetUnit()
  const
{
  if (!fDomNode) {
    WARNING(fgWarning);
    ERROR("Getting unit from null-Branch");
    exit(EXIT_FAILURE);
  }

  double unit = 1.;
    
  // Check for a unit attribute
  const DOMNamedNodeMap* const attributes = fDomNode->getAttributes();

  if (attributes) {

    string unitString;
    for (unsigned int j = 0; j < attributes->getLength(); ++j) {

      const DOMNode* const attribute = attributes->item(j);
      const string name = StrX(attribute->getNodeName()).localForm();
      if (name == "UNIT" || name == "unit") {
	unitString = StrX(attribute->getNodeValue()).localForm();
	break;
      }
    }

    // If found, convert the unit to the appropriate scale factor,
    // otherwise set scale factor to 1.
    if (!unitString.empty()) {
      // HepTool::Evaluator eval;
      unit = gEval.evaluate(unitString.c_str());

      if (gEval.status() != HepTool::Evaluator::OK) {
	gEval.print_error();
	ostringstream msg;
	msg << "The unit '" << unitString << "' "
	       "was not understood by the expression evaluator of CLHEP" << endl;
	ERROR(msg);
	exit(EXIT_FAILURE);	    
      }
    }
  }

  return unit;
}


Branch&
Branch::operator=(const Branch& b)
{
  fDomNode = b.fDomNode;
  return *this;
}


//! Getting a std::string
void
Branch::GetData(string& s)
  const
{
  // do not use castData, as in this case
  // we want to return the ENTIRE string, including
  // possible white space.
  s = GetDataString();
}


// DV considered dangeres since user has to call free()
//! Getting a char*
void
Branch::GetData(char*& c)
  const
{
  // do not use castData, as in this case
  // we want to return the ENTIRE string, including
  // possible white space. Do not forget to call free.

  /*const string s = GetDataString();
  c = strdup(s.c_str());*/
  ERROR("This method requires user to call free() after "
        "done with the char*! Disabled for the moment.");
  exit(1);
}


//! Getting an utl::TimeStamp
/*void
Branch::GetData(TimeStamp& time)
  const
{
  Reader::GetXMLTime(GetDataString(), time);
}


void
Branch::GetData(vector<TimeStamp>& vt)
  const
{
  const string dataString = GetDataString();
  istringstream is(dataString);
  string str;
  TimeStamp t;
  while (!is.eof())
    if (is >> str) {
      Reader::GetXMLTime(str, t);
      vt.push_back(t);
    }
}
*/

#define UTL_BRANCH_GETDATA_WITH_CAST_ONLY(_Type_...) \
void                                                 \
Branch::GetData(_Type_& t)                           \
  const                                              \
{                                                    \
  CastData(t);                                       \
}


UTL_BRANCH_GETDATA_WITH_CAST_ONLY(bool)
UTL_BRANCH_GETDATA_WITH_CAST_ONLY(vector<bool>)
UTL_BRANCH_GETDATA_WITH_CAST_ONLY(list<bool>)
UTL_BRANCH_GETDATA_WITH_CAST_ONLY(vector<string>)
UTL_BRANCH_GETDATA_WITH_CAST_ONLY(list<string>)


#undef UTL_BRANCH_GETDATA



namespace utl {

  ostream&
  operator<<(ostream& os, const Branch& b)
  {
    // For now, just prints the branch name.
    // Could have more detailed information later.
  
    os << b.GetBranchNameString() << endl;
    return os;
  }

}


// Configure (x)emacs for this file ...
// Local Variables:
// mode: c++
// End:
