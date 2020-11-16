
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Spectrum.h *                                  galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#ifndef Spectrum_h
#define Spectrum_h

#include <valarray>

struct Spectrum {

  //public: 
  size_t nElems;
  //double* s; //AWS20050624
  std::valarray<double> s;

Spectrum() : nElems(0) {}

Spectrum(const size_t n, const double val = 0.) : nElems(n) {

  //s = new double[nElems];
  s.resize(nElems);
  for (size_t i = 0; i < nElems; ++i)
    s[i] = val;

}

  Spectrum(const Spectrum& spec) { 

    //delete[] s; 
    nElems = spec.nElems; 
    //s = new double[nElems]; 
    s.resize(nElems);
    s = spec.s;
    //for (size_t i = 0; i < nElems; ++i) 
    //s[i] = spec.s[i];

  }

  Spectrum(const std::valarray<double>& spec) { 

    //delete[] s; 
    nElems = spec.size(); 
    //s = new double[nElems]; 
    s.resize(nElems);
    s = spec;
    //for (size_t i = 0; i < nElems; ++i) 
    //s[i] = spec[i]; 

  }

  Spectrum& operator=(const Spectrum& spec) { 
  
    if (this != &spec) { 
    
      //delete[] s; 
      nElems = spec.nElems; 
      //s = new double[nElems]; 
      s.resize(nElems);
      s = spec.s;
      //for (size_t i = 0; i < nElems; ++i) 
      //s[i] = spec.s[i]; 
    
    } 
    
    return *this; 
  
  }
   
    Spectrum& operator=(const double val) {

      for (size_t i = 0; i < nElems; ++i)
	s[i] = val;

      return *this;

    }

    const std::valarray<double>& GetSpec() { 
      
      return s;

      //std::valarray<double> result(0., nElems); 
      //for (size_t i = 0; i < nElems; ++i) 
      //result[i] = s[i]; 
      //return result; 
    
    }

  ~Spectrum() { 
    
    //delete[] s; 
    nElems = 0; 

  }

};

#endif
