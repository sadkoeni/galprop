#ifndef _utl_Reader_imp_h_
#define _utl_Reader_imp_h_

/*!
  \file Reader_imp.h

  \brief Implementation file for specialization of the template member
  functions in the Reader class
  
  \author T. Paul
  \author P. Cattaneo
  \version $Id: Reader_imp.h,v 1.4 2004/10/30 01:53:03 veberic Exp $
*/

namespace utl {

  //! Getting a boolean
  template<>
  void
  Branch::GetData(bool& b)
    const
  {
    // note we don't bother with the unit here.
    CastData(b);
  }


  //! Getting a std::string
  template<>
  void
  Branch::GetData(std::string& s)
    const
  {
    // do not use castData, as in this case
    // we want to return the ENTIRE string, including
    // possible white space.
    s = GetDataString();
  }


  //! Getting a char*
  template<>
  void
  Branch::GetData(char*& c)
    const
  {
    // do not use castData, as in this case
    // we want to return the ENTIRE string, including
    // possible white space.
    const std::string s = GetDataString();
    c = strdup(s.c_str());
  }


  //! Getting a vector of string
  template<>
  void
  Branch::GetData(std::vector<std::string>& s)
    const
  {
    // ignore units (if they exist) in this case
    CastData(s);
  }


  //! Getting a vector of bool
  template<>
  void
  Branch::GetData(std::vector<bool>& b)
    const
  {
    CastData(b);
  }


  //! Getting a list of string
  template<>
  void
  Branch::GetData(std::list<std::string>& s)
    const
  {
    // ignore units (if they exist) in this case
    CastData(s);
  }


  //! Getting a list of bool
  template<>
  void
  Branch::GetData(std::list<bool>& b)
    const
  {
    CastData(b);
  }
 
} // namespace utl

#endif

// Configure (x)emacs for this file ...
// Local Variables:
// mode: c++
// End:
