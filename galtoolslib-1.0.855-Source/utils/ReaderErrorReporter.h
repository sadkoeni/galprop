/*
 * The Apache Software License, Version 1.1
 *
 * Copyright (c) 1999-2002 The Apache Software Foundation.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *    if any, must include the following acknowledgment:
 *       "This product includes software developed by the
 *        Apache Software Foundation (http://www.apache.org/)."
 *    Alternately, this acknowledgment may appear in the software itself,
 *    if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Xerces" and "Apache Software Foundation" must
 *    not be used to endorse or promote products derived from this
 *    software without prior written permission. For written
 *    permission, please contact apache\@apache.org.
 *
 * 5. Products derived from this software may not be called "Apache",
 *    nor may "Apache" appear in their name, without prior written
 *    permission of the Apache Software Foundation.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.  IN NO EVENT SHALL THE APACHE SOFTWARE FOUNDATION OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of the Apache Software Foundation, and was
 * originally based on software copyright (c) 1999, International
 * Business Machines, Inc., http://www.ibm.com .  For more information
 * on the Apache Software Foundation, please see
 * <http://www.apache.org/>.
 */

#ifndef _ReaderErrorReporter_h_
#define _ReaderErrorReporter_h_

#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <iostream>

namespace utl {
  
  class ReaderErrorReporter : public xercesc::ErrorHandler {
  public:
    
  ReaderErrorReporter() : fSawErrors(false) { }
    
    // -----------------------------------------------------------------------
    //  Implementation of the error handler interface
    // -----------------------------------------------------------------------
    void warning(const xercesc::SAXParseException& toCatch);
    void error(const xercesc::SAXParseException& toCatch);
    void fatalError(const xercesc::SAXParseException& toCatch);
    void resetErrors();
    
    // -----------------------------------------------------------------------
    //  Getter methods
    // -----------------------------------------------------------------------
    bool getSawErrors() const { return fSawErrors; }
    
    // -----------------------------------------------------------------------
    //  Private data members
    //
    //  fSawErrors
    //      This is set if we get any errors, and is queryable via a getter
    //      method. Its used by the main code to suppress output if there are
    //      errors.
    // -----------------------------------------------------------------------
    bool fSawErrors;
  };

  // ---------------------------------------------------------------------------
  //  This is a simple class that lets us do easy (though not terribly efficient)
  //  trancoding of strings to XMLCh data.
  // ---------------------------------------------------------------------------
  class XStr {
  public:
    
  XStr(const std::string &toTranscode)
    : fUnicodeForm(0) {
      
      // Call the private transcoding method
      fUnicodeForm = xercesc::XMLString::transcode(toTranscode.c_str());
    }
    
    ~XStr() {
      
      if (fUnicodeForm)
	xercesc::XMLString::release(&fUnicodeForm);
    }
    
    // -----------------------------------------------------------------------
    //  Getter methods
    // -----------------------------------------------------------------------
    const XMLCh* unicodeForm() const { return fUnicodeForm; }
    
  private:
    // since we have pointer data-members we should be careful about copying:
    XStr(const XStr&);
    XStr& operator=(const XStr&);
    // -----------------------------------------------------------------------
    //  Private data members
    //
    //  fUnicodeForm
    //      This is the Unicode XMLCh format of the string.
    // -----------------------------------------------------------------------
    XMLCh* fUnicodeForm;
    
  };
  
  // ---------------------------------------------------------------------------
  //  This is a simple class that lets us do easy (though not terribly efficient)
  //  trancoding of XMLCh data to local code page for display.
  // ---------------------------------------------------------------------------
  class StrX {
  public:
    
  StrX(const XMLCh* const toTranscode)
    : fLocalForm(0) {
      
      // Call the private transcoding method
      fLocalForm = xercesc::XMLString::transcode(toTranscode);
    }
    
    ~StrX() {
      
      if (fLocalForm)
	xercesc::XMLString::release(&fLocalForm);
    }
    
    // -----------------------------------------------------------------------
    //  Getter methods
    // -----------------------------------------------------------------------
    const char* localForm() const { return fLocalForm; }
    
  private:
    // since we have pointer data-members we should be careful about copying:
    StrX(const StrX&);
    StrX& operator=(const StrX&);
    // -----------------------------------------------------------------------
    //  Private data members
    //
    //  fLocalForm
    //      This is the local code page form of the string.
    // -----------------------------------------------------------------------
    char* fLocalForm;
    
  };
  
  inline std::ostream&
    operator<<(std::ostream& target, const StrX& toDump)
  { return target << toDump.localForm(); }
  
}

#endif
