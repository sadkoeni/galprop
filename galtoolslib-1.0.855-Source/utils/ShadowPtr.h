#ifndef _ShadowPtr_h_
#define _ShadowPtr_h_

#include <sys/types.h>
#include <sstream>
#include <iostream>
#ifdef __GNUC__
#include <cxxabi.h>
#endif

#include <ShadowPtr_fwd.h>
  
namespace utl {
  
  template<typename Type>
    class Meta {
    
  public:
    typedef class { } False;
    typedef class { False f[256]; } True;
    
    /* DV: note that the signature of a pointer to virtual const method Clone()
       is T* (T::*)() const
       vvvvvvvvvvvvvvvvv                                   */
    /*template<typename T, T* (T::*)() const>
      struct MatchCloneSignature {
      typedef True type;
      };
      
      static False MatchClone(...);
      
      template<typename T>
      static typename MatchCloneSignature<T, &T::Clone>::type
      MatchClone(const T&);
      
      template<typename T>
      struct HasClone {
      static const bool value =
      sizeof(MatchClone(*(T*)0)) == sizeof(True);
      };*/
    
    template<typename T>
      static True MatchCloneTag(typename T::IsClonableTag* const);
    
    template<typename T>
      static False MatchCloneTag(...);
    
    template<typename T>
      struct HasClone {
	static const bool value =
	  sizeof(Meta<Type>::template MatchCloneTag<T>(0)) == sizeof(True);
      };
    
    template<bool v, int = 0>
      struct BoolSwitch {
	static Type* GetCopy(const Type& t) { return new Type(t); }
      };
    
    template<int i>
      struct BoolSwitch<true, i> {
      static Type* GetCopy(const Type& t) { return t.Clone(); }
    };
    
  public:
    static Type* GetCopy(const Type& t)
    { return BoolSwitch<HasClone<Type>::value>::GetCopy(t); }
  
  };
  
  template<typename Type>
    class ShadowPtr {
 
  public:
  ShadowPtr() : fSpook(0) { }
    
    // shallow ctor
    explicit ShadowPtr(Type* const phantom) : fSpook(phantom) { }
    
    // deep copy ctor
    explicit ShadowPtr(const ShadowPtr& shadow)
    { DeepCopy(shadow.fSpook); }
    
    ~ShadowPtr() { Delete(); }
    
    Type* Get() { return fSpook; }
    const Type* Get() const { return fSpook; }
    
    // deep copy assignment
    ShadowPtr& operator=(const ShadowPtr& shadow)
      { Delete(); DeepCopy(shadow.fSpook); return *this; }
    
    // shallow assignment
    ShadowPtr& operator=(Type* const spook)
      { Delete(); fSpook = spook; return *this; }
    
    Type& operator*() { CheckPtr(); return *fSpook; }
    const Type& operator*() const { CheckPtr(); return *fSpook; }
    
    Type* operator->() { CheckPtr(); return fSpook; }
    const Type* operator->() const { CheckPtr(); return fSpook; }
    
    bool operator==(const ShadowPtr& shadow) const { return fSpook == shadow.fSpook; }
    bool operator!=(const ShadowPtr& shadow) const { return fSpook != shadow.fSpook; }
    bool operator==(const Type* const spook) const { return fSpook == spook; }
    bool operator!=(const Type* const spook) const { return fSpook != spook; }
    operator bool() const { return fSpook; }
    bool operator!() const { return !fSpook; }
    
  private:
    void CheckPtr() const {
      
      if (!fSpook) {
	std::stringstream err;
	err << "Dereferencing zero pointer to class ";
#ifdef __GNUC__
	int status;
	err << abi::__cxa_demangle(typeid(Type).name(), 0, 0, &status);
#else
	err << typeid(Type).name();
#endif
	std::cerr << err.str() << std::endl;
	//throw utl::NonExistentComponentException(err.str());
      
      }
     
    }
 
    void Delete() { delete fSpook; }
 
  protected:
    void DeepCopy(const Type* const spook)
    { fSpook = (spook ? Meta<Type>::GetCopy(*spook) : 0); }
    
    Type* fSpook;
    
  };
  
  template<typename Type>
    class InitializedShadowPtr : public ShadowPtr<Type> {
    
  public:
  InitializedShadowPtr() : ShadowPtr<Type>(new Type) { }
    
  };
  
  
  template<typename T>
    bool DeepEqual(const ShadowPtr<T>& s1, const ShadowPtr<T>& s2)
    { return (!s1 && !s2) || (s1 && s2 && *s1 == *s2); }
  
}

#endif
