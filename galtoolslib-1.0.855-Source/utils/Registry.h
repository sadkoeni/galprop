#ifndef _utl_Registry_h_
#define _utl_Registry_h_

#include <string>
#include <map>
#include <memory>
#include <ErrorLogger.h>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <memory>
#include <utility>

namespace utl {

   /** \brief Template class to implement a class registry for subclasses
    *
    * No instance is required, all the methods and registrations are static.
    * To register a subclass of a baseclass, you create an instance of
    * utl::Registry0<baseclass>::Registrar<subclass>( name )
    * Then an instance can be created with
    * utl::Registry0<baseclass>::create(name)
    */
  template <typename T>
    class Registry0 {
    
  public:
    static std::unique_ptr<T> create ( const std::string &name ) {
      typename std::map<std::string, Registry0<T>::create_fp>::iterator it = creators().find(name);
      
      if ( it == creators().end() ) {
	std::ostringstream os;
	os << "Error in registry, no subclass named \""<<name<<"\" registered for class \""<<typeid(T).name()<<"\"."<<std::endl;
	os << "Available names are:"<<std::endl;
	for (auto x : creators())// it = creators().begin(); it != creators().end(); ++it )
	  os << "   "<< x.first << std::endl;
	FATAL(os.str());
	throw(std::runtime_error("Name not found in registry"));
      }
      
      return (it->second)();
    
    }
    
    typedef std::unique_ptr<T> (*create_fp) ( );
    
    template <typename C>
      class Registrar {
    public:
      explicit Registrar (const std::string &name) {
	Registry0<T>::Register(name, &Registry0<T>::Registrar<C>::createNew);
      }
      static std::unique_ptr<T> createNew() {
	return std::unique_ptr<T> (new C);
      }
    };
    
    static void Register(const std::string &name, create_fp fp) {
      typename std::map<std::string, Registry0<T>::create_fp>::iterator it = creators().find(name);
      
      if ( it != creators().end() ) {
	std::ostringstream os;
	os << "Error in registry, subclass named \""<<name<<"\" already registered for class \""<<typeid(T).name()<<"\".";
	FATAL(os.str());
	throw(std::runtime_error("Name already in registry"));
      }
      
      creators()[name] = fp;
    }
    
    static std::vector<std::string> GetRegisteredNames() {
      std::vector<std::string> names;
      for (auto x : creators())
    	names.push_back(x.first);
      return names;
    }

    static void QueryRegistered() {
      std::ostringstream os;
      os << "Available names are:"<<std::endl;
      for (auto x : creators())
	os << "   "<< x.first << std::endl;
      INFO(os.str());
    }

  private:
    static std::map<std::string, Registry0<T>::create_fp> & creators() {
      static std::map<std::string, Registry0<T>::create_fp> impl;
      return impl;
    }
    
  };
  
  /** \brief Template class to implement a class registry for subclasses that take a single argument in the creator
   *
   * No instance is required, all the methods and registrations are static.
   * To register a subclass of a baseclass that requires an argument in the creator, 
   * you create an instance of
   * utl::Registry1<baseclass,argument>::Registrar<subclass>( name )
   * Then an instance can be created with
   * utl::Registry1<baseclass,argument>::create(name)
   */
  /*
  template <typename T, typename I>
    class Registry1 {
    
  public:
    typedef std::unique_ptr<T> (*create_fp) ( I );
    
    static std::unique_ptr<T> create ( const std::string &name, I i ) {
      
      typename std::map<std::string, create_fp>::iterator it = creators().find(name);
      
      if ( it == creators().end() ) {
	std::ostringstream os;
	os << "Error in registry, no subclass named \""<<name<<"\" registered for class \""<<typeid(T).name()<<"\"."<<std::endl;
	os << "Available names are:"<<std::endl;
	for (auto x : creators())// it = creators().begin(); it != creators().end(); ++it )
	  os << "   "<< x.first << std::endl;
	FATAL(os.str());
	throw(std::runtime_error("Name not found in registry"));
      }
      
      return (it->second)(std::forward<I>(i));
    }
    
    template <typename C>
      class Registrar {
    public:
      explicit Registrar (const std::string &name) {
	Registry1<T,I>::Register(name, &Registrar<C>::createNew);
      }
      static std::unique_ptr<T> createNew(I i) {
	return std::unique_ptr<T>(new C(std::forward<I>(i)));
      }
    };
    
    static void Register(const std::string &name, create_fp fp) {
      typename std::map<std::string, create_fp>::iterator it = creators().find(name);
      
      if ( it != creators().end() ) {
	std::ostringstream os;
	os << "Error in registry, subclass named \""<<name<<"\" already registered for class \""<<typeid(T).name()<<"\".";
	FATAL(os.str());
	throw(std::runtime_error("Name already in registry"));
      }
      
      creators()[name] = fp;
    }
    
    static std::vector<std::string> GetRegisteredNames() {
      std::vector<std::string> names;
      for (auto x : creators())
    	names.push_back(x.first);
      return names;
    }

    static void QueryRegistered() {
      std::ostringstream os;
      os << "Available names are:"<<std::endl;
      for (auto x : creators())
	os << "   "<< x.first << std::endl;
      INFO(os.str());
    }

  private:
    static std::map<std::string, Registry1<T,I>::create_fp> & creators() {
      static std::map<std::string, create_fp> impl;
      return impl;
    }
    
  };
  */

  /** \brief Template class to implement a class registry for subclasses that take a variable number of arguments in the creator
   *
   * No instance is required, all the methods and registrations are static.
   * To register a subclass of a baseclass that requires an argument in the creator, 
   * you create an instance of
   * utl::RegistryN<baseclass,argument1,argument2,...>::Registrar<subclass>( name )
   * Then an instance can be created with
   * utl::RegistryN<baseclass,argument1,argument2,...>::create(name)
   */
  template <typename T, typename... Args>
    class RegistryN {
    
  public:
    using create_fp = std::unique_ptr<T> (*) ( Args... args );
    
    static std::unique_ptr<T> create ( const std::string &name, Args... args ) {
      
      const auto it = creators().find(name);
      
      if ( it == creators().end() ) {
	std::ostringstream os;
	os << "Error in registry, no subclass named \""<<name<<"\" registered for class \""<<typeid(T).name()<<"\"."<<std::endl;
	os << "Available names are:"<<std::endl;
	for (auto x : creators())// it = creators().begin(); it != creators().end(); ++it )
	  os << "   "<< x.first << std::endl;
	FATAL(os.str());
	throw(std::runtime_error("Name not found in registry"));
      }
      
      //return (it->second)(args...);
      return (it->second)(std::forward<Args>(args)...);
    }
    
    template <typename C>
      class Registrar {
    public:
      explicit Registrar (const std::string &name) {
	RegistryN<T,Args...>::Register(name, &Registrar<C>::createNew);
      }
      static std::unique_ptr<T> createNew(Args... args) {
	//return std::unique_ptr<T>(new C(args...));
	return std::unique_ptr<T>(new C(std::forward<Args>(args)...));
      }
    };
    
    static void Register(const std::string &name, create_fp fp) {
      auto it = creators().find(name);
      
      if ( it != creators().end() ) {
	std::ostringstream os;
	os << "Error in registry, subclass named \""<<name<<"\" already registered for class \""<<typeid(T).name()<<"\".";
	FATAL(os.str());
	throw(std::runtime_error("Name already in registry"));
      }
      
      creators()[name] = fp;
    }
    
    static std::vector<std::string> GetRegisteredNames() {
      std::vector<std::string> names;
      for (auto x : creators())
    	names.push_back(x.first);
      return names;
    }

    static void QueryRegistered() {
      std::ostringstream os;
      os << "Available names are:"<<std::endl;
      for (auto x : creators())
	os << "   "<< x.first << std::endl;
      INFO(os.str());
    }

  private:
    static std::map<std::string, RegistryN<T,Args...>::create_fp > & creators() {
      static std::map<std::string, create_fp > impl;
      return impl;
    }
    
  };

  template <typename T, typename I>
     using Registry1 = RegistryN<T, I>;

}

#endif
