#ifndef _Reader_h_
#define _Reader_h_

/*!
  \file Reader.h

  \brief Utility for reading data from XML files
  
  \author T. Paul
  \author P. Cattaneo
  \author D. Veberic
  \version $Id: Reader.h 6230 2007-08-14 20:39:20Z darko $
*/

// XML includes
#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMImplementation.hpp> 
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/sax/HandlerBase.hpp>

#include <cstdio>
#include <sstream>

// stuff we want to be able to GetData()
#include <vector>
#include <list>
#include <map>
#include <utility>
#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>

#include <ShadowPtr.h>

namespace utl {

  //class TimeStamp;


  //! Predicate used in STL for searching for whitespace
  struct IsNotSpace {

    bool operator()(const char alpha)
    { return alpha != ' ' && alpha != '\015' && alpha != '\012' && alpha != '\011'; }

  };

  //! Predicate used in STL for searching for whitespace
  struct IsSpace {

    bool operator()(const char alpha)
    { return !IsNotSpace()(alpha); }

  };


  /**
     \class Branch

     \brief Class representing a document branch

     For an introduction to using Branches, see the documentation for the Reader.

     A branch is essentially a DOM element node (a subset of what the
     DOM considers to be a node). That is, DOM nodes such as comment
     or text nodes which may not have children do not qualify to be
     branches. The idea of the Branch is to provide a simple tool for
     navigating trees of data that might be more palatable for Auger
     applications (ie those which do not require especially
     sophisticated document traversal.)

     In any case, if you do not want to use the Branch class, you are
     free to retrieve the DOM document element via
     Reader::GetDocument() and apply standard DOM traversal tools.
     Thus the reader allows the full power of DOM2 should the user
     need it.

     \version $Id: Reader.h 6230 2007-08-14 20:39:20Z darko $
     \ingroup parsers
  */

  class Branch {

  public:
    Branch() : fDomNode(0) { }

    Branch(const Branch& branch) : fDomNode(branch.fDomNode) { }

    Branch(xercesc::DOMNode* const dom) : fDomNode(dom) { }

    //! Get first child of this Branch 
    /*! \return NullBranch if no child found */
    Branch GetFirstChild() const;
  
    //! Get child of this Branch by child name
    /*! \return NullBranch if no child found */
    Branch GetChild(const std::string& childName) const;
  
    //! Get child of this Branch by child name and an attribute called "id"
    /*! \return NullBranch if no child found */
    Branch GetChild(const std::string& childName, const std::string& multiID) const;

    //! Get child of this Branch by child name and any number of attribute-value pairs.
    /*! This method allows the user to specify values for any number
      of attributes. A (non-null) Branch is returned only if the
      requested name is found in the XML file and <b>all</b> the
      attributes are found <b>with</b> their specified value, <b>with
      the single exception of the unit attribute</b>.  Unit is handled
      as a special case by the reader, and any unit attribute which
      may be present, either in the XML file or the second argument of
      this method, will be ignored by this method. The requested
      attributes are specified in a map\<string, string\> where the
      first string is the attribute name and the second string is the
      attribute value.
    */
    Branch GetChild(const std::string& childName,
		    std::map<std::string, std::string> attributeMap) const;

    //! Get next sibling of this branch
    /*! \return Null-Branch if no child found */
    Branch GetNextSibling() const;
  
    //! Get sibling by name
    /*! \return Null-Branch if no sibling found */
    Branch GetSibling(const std::string& childName) const;
  
    //! Get sibling by name and ID
    /*! \return Null-Branch if no sibling found */
    Branch GetSibling(const std::string& childName, const std::string& multiID) const;

    //! Get sibling of this Branch by child name and any number of attribute-value pairs.
    /*! \return Null-Branch if no sibling found */
    Branch GetSibling(const std::string& childName,
		      std::map<std::string, std::string>& attributeMap) const;

    //! Get a map\<list, list\> containing all the attributes of this Branch
    /*! Unlike the GetData() method, GetAttributes() makes no attempt
      to cast the attributes */
    std::map<std::string, std::string> GetAttributes() const;

    // GetData methods.  These methods find data in the current branch
    // and attempt to cast it according to the type of argument in the
    // GetData argument list.  Note that, if you look at the actual
    // implementation of GetData methods, most of them simply invoke
    // the (private) castData template function.  The reason for the
    // intermediate GetData method is that is some cases, one might
    // want to deal with special cases that are not dealt with by a
    // simple template.  For example, suppose one has an XML element:
    // <someData> 13 14 15 16 </someData> If the user requests to
    // retrieve <someData> as an int, probably one should return just
    // 13.  If one asks for a string, probably one should return the
    // whole list of numbers.  These two different interpretations can
    // be dealt with using the intermediate GetData methods to
    // compliment the templated castData method.  Note that,
    // currently, scaling by the appropriate unit factor is handled in
    // the GetData methods, not the castData method.  (in principle
    // this could change, but keep in mind that for the case of
    // strings you don't want to try to multiply by a scale factor.)

    //! Overloads of the GetData member template function
    // for cases that do not need unit conversion

    // bool
    void GetData(bool& b) const;
    void GetData(std::vector<bool>& b) const;
    void GetData(std::list<bool>& b) const;
    // string
    void GetData(std::string& s) const;
    void GetData(std::vector<std::string>& s) const;
    void GetData(std::list<std::string>& s) const;
    // char*
    void GetData(char*& c) const;
    // TimeStamp
    //void GetData(utl::TimeStamp& t) const;
    //void GetData(std::vector<utl::TimeStamp>& vt) const;
    
    /// Get data in the current branch into an atomic type.
    /*! Data in the Branch gets cast to type T */
    template<typename T> void GetData(T& a) const
    { CastData(a); a *= static_cast<T>(GetUnit()); }
  
    //! Get data in the current Branch into an STL list or vector.
    /*! Data are loaded into an STL container of type W. Atomic types
     *  in the XML Branch can be space, tab, or CR delimited. For
     *  example, data can be provided in an XML file like so:
     *
     *  \code
     *  <someData> 12.0 16.1 18.4 </someData>
     *  \endcode
     *
     * Since there are 3 space delimited floating point numbers
     * between the \<someData\> tags, one could read these data with
     * something like:
     *
     * \code
     * vector<double> someData;
     * someDataBranch.GetData(someData);
     * // someDataBranch assumed to point to <someData> element
     * assert(someData.size() == 3);
     * // the vector should have three numbers in it
     * \endcode
     */
    template<typename T, class A, template<typename, typename> class W>
    void
    GetData(W<T, A>& a)
      const
    {
      CastData(a);
      T u = static_cast<T>(GetUnit());    
      for (typename W<T, A>::iterator it = a.begin(); it != a.end(); ++it)
	*it *= u;
    }

    /// Get data in the current branch into a pair\<\>
    template<typename T1, typename T2>
    void
    GetData(std::pair<T1, T2>& p)
      const
    {
      std::istringstream is(GetDataString());
      is >> p.first >> p.second;
      p.first  *= static_cast<T1>(GetUnit());
      p.second *= static_cast<T2>(GetUnit());
    }

    /// Retrieve branch name as a string  
    std::string GetBranchNameString() const;

    xercesc::DOMNode* GetDOMNode() { return fDomNode; }

    Branch& operator=(const Branch& b);

    // the following two operators are for backward compatibility for the
    // constructs like "Branch == NULL" and "Branch != NULL"
    // note that compiler warning is issued with -Wall since this
    // is only half legal (comparing an object to a pointer)
    // all comparison with NULL will have to be removed in favour of
    // "if (Branch)" or "if (!Branch)"
    template<typename T> bool operator==(T* p) const
    { return fDomNode == p; }

    template<typename T> bool operator!=(T* p) const
    { return fDomNode != p; }

    /// Branches are equal if they point to the same memory (underlying DOMNode)
    /*! Note that branches can have the same contents and be different branches */
    bool operator==(Branch& b) const
    { return fDomNode == b.GetDOMNode(); }

    /// Branches are not equal if they point to different memory
    bool operator!=(Branch& b) const
    { return fDomNode != b.GetDOMNode(); }

    // the "legal" operators
    operator bool() const
    { return fDomNode; }

    bool operator!() const
    { return !fDomNode; }
  
  private:
    mutable xercesc::DOMNode* fDomNode;

    static std::string fgWarning;

    // Return data for the Branch, attempt to cast it as
    // the requested type, and multiply by the unit multiplier.
    // These methods use the getDataAndUnit method to find the
    // data string and unit multiplier
    template<typename T>
    void CastData(T& dataT) const {
      std::istringstream is(GetDataString());
      is >> dataT;
    }
  
    template<typename T, class A, template<typename, typename> class W>
    void
    CastData(W<T, A>& v)
      const
    {
      const std::string dataString = GetDataString();
      std::istringstream is(dataString);
      T value;
      while (!is.eof())
	if (is >> value)
	  v.push_back(value);
    }

    // helper function to get the data inside an element as one big string
    std::string GetDataString() const;
  
    // helper function to the (optional) unit attribute and return the
    // appropriate scale factor.
    double GetUnit() const;

    friend class Reader;
    friend std::ostream& operator<<(std::ostream& os, const Branch b);

  };


  /*!
    \class ReaderStringInput

    \brief This just defines a type which holds some character data
           to be parsed by the Reader.  

    This class is for use in cases where the Reader should parse a
    string of XML information in memory instead of information in a
    file.  It takes care of some of the pain of using the
    MemBufInputSource class of Xerces.
  */
 
  class ReaderStringInput {

  public:
    ReaderStringInput() { }

    ReaderStringInput(const std::string& inputString)
    { fInputString = inputString; }

    ReaderStringInput(const std::ostringstream& inputStream)
    { fInputString = inputStream.str(); }

    ~ReaderStringInput() { }
 
    const std::string& GetInputString() { return fInputString; }
   
  private:
    std::string fInputString;

  };


  class Reader {

  public:
    enum Validation {
      eDTD,
      eSCHEMA,
      eNONE
    };

    //! Constructor with arguments for XML file name and an optional flag
    // to switch on or off validation
    Reader(const std::string&, bool validation);

    Reader(const std::string&, Validation validationType = Reader::eNONE);

    Reader(const ReaderStringInput& inputString, Validation vaidationType = Reader::eNONE);

    //! Get the top Branch (represents same entity as document node)
    Branch GetTopBranch() { return fTopBranch; }

    //! Get the document node
    xercesc::DOMDocument* GetDocument() { return fDocumentNode; }

    //! Utility to make ASCII dump of the XML tree
    void DumpTree();

    //! Utility to convert strings from XML Schema dateTime format
    /*! into separate integers giving date and time constituents */
    /*static void BreakdownXMLTime(const std::string& xmlTime,
                                 unsigned int& year,
                                 unsigned int& month,
                                 unsigned int& day,
                                 unsigned int& hour,
                                 unsigned int& minute,
                                 double&       second,
                                 int&          utmOffsetHour,
                                 unsigned int& utmOffsetMinute);

    static void GetXMLTime(const std::string& str, utl::TimeStamp& ts);
    */
  private:
    Reader(const Reader&);
    Reader& operator=(const Reader&);

    void SetUp(const Validation v);

    std::string fXMLFileName;         // For case of reading from a file
    ReaderStringInput fMemInput;      // For case of reading from an input string

    ShadowPtr<xercesc::XercesDOMParser> fParser;
    ShadowPtr<xercesc::ErrorHandler> fErrReporter;

    // Name of XML "data-cards" (for the FFREAD-style API)
    std::string fCardname;
    // Branch where the XML "data-cards" live (for the FFREAD-style API)
    Branch fDataCardBranch;

    // root node of the document, owned by the parser, do not delete
    xercesc::DOMDocument* fDocumentNode;

    // set type of validation requested (DTD, Schema, or None)
    Validation fValidationMode;

    // Branch at the top of the document tree.
    Branch fTopBranch;

    void Initialize();
    void Parse();

    // helper function for dumping out the XML tree (recursive)
    void DumpMe(const Branch& b);
    std::string fSpaces;

  };


} // namespace utl


#endif

// Configure (x)emacs for this file ...
// Local Variables:
// mode: c++
// compile-command: "make -C .. -k"
// End:
