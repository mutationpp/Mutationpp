/**
 * @file AutoRegistration.h
 *
 * Provides classes which implement automatic registration of polymorphic types
 * from a name which identifies a subtype.  The code in this file has been 
 * entirely taken from the SConfig library written by Andrea Lani and the
 * documentation and formatting has only been altered to match the standards
 * used in the Mutation++ library.  In addition, the use of exceptions has been
 * removed from the original code and replaced with simple cout statements.
 *
 * @author Andrea Lani
 * @author J.B. Scoggins
 */
#ifndef AUTO_REGISTRATION_H
#define AUTO_REGISTRATION_H

#include <iostream>
#include <map>
#include <string>

namespace Mutation {
    namespace Utilities {
        namespace Config {

template <class PTYPE> class Provider;

/**
 * This class stores an archive of available Provider's for a generic 
 * polymorphic type PTYPE. Once registered, they can be accessed by name.
 *
 * @author Andrea Lani
 */
template <class PTYPE>
class Factory 
{
public:
  
    /**
     * Resister the Provider.
     */
    void add(Provider<PTYPE>* ptr)
    {
        if (m_providers.count(ptr->getName()) == 0) {
            m_providers[ptr->getName()] = ptr;
        } else {
            std::cout << "Provider < " << ptr->getName()
                      << " > is already registered" << std::endl;            
        }
    }
  
    /**
     * Get the provider corresponding to the given key.
     */
    Provider<PTYPE>* getProvider(const std::string& key)
    {    
        if (m_providers.count(key) == 0) {
            std::cout << "Provider < " <<  key << " > is not registered" 
                      << std::endl;
        }
        
        return dynamic_cast<Provider<PTYPE>*>(
            m_providers.find(key)->second);
    }
  
    /**
     * Instance of the singleton Factory.
     */
    static Factory<PTYPE>& getInstance() 
    {
        static Factory<PTYPE> f;
        return f;
    }
    
    /**
     * Returns a pointer to an object of type registered by the given name.
     */
    static PTYPE* create(const std::string& name, typename PTYPE::ARGS args)
    {
        return getInstance().getProvider(name)->create(args);
    }
  
private:
  
    /**
     * Constructor.
     */ 
    Factory() {}
  
    /**
     * Destructor.
     */
    ~Factory() {}
  
private: // data
  
  // archive of providers with type <name, Provider*>
  std::map<std::string, Provider<PTYPE>*> m_providers;
  
}; // class Factory


/**
 * This class provides a name to a deriving object.
 *
 * @author Andrea Lani
 */
class NamedObject 
{
public:
  
    /**
     * Constructor.
     */
    NamedObject(const std::string& name) 
        : m_name(name) 
    { }
  
    /**
     * Destructor.
     */
    ~NamedObject() {}
  
    /**
     * Get the name of the object.
     */
    std::string getName() const {
        return m_name;
    }
  
private:
  
    // name of the object
    std::string m_name;
    
}; // class NamedObject


/**
 * This class implements an object provider.
 *
 * @author Andrea Lani
 */
template <class PTYPE>
class Provider : public NamedObject 
{
public:
  
    /**
     * Constructor.
     */
    Provider(const std::string& name) 
        : NamedObject(name)
    {
        Factory<PTYPE>::getInstance().add(this);
    }
  
    /**
     * Destructor.
     */ 
  virtual ~Provider() {}
  
  /**
   * Create the chosen concrete object.
   */
  virtual PTYPE* create(typename PTYPE::ARGS args) = 0;
  
}; // class Provider


/**
 * This class implements a generic provider allowing client code to create 
 * arbitrary types of polymorphic objects. 
 *
 * @author Andrea Lani
 */
template <class CTYPE, class PTYPE>
class ObjectProvider : public Provider<PTYPE>
{
public:
  
    /**
     * Constructor.
     */
    ObjectProvider(const std::string& name) 
        : Provider<PTYPE>(name) 
    { }
  
    /**
     * Destructor.
     */
    ~ObjectProvider() {}
  
    /**
     * Create an object of polymorphic type PTYPE (parent) and static type 
     * CTYPE (concrete).
     */
    PTYPE* create(typename PTYPE::ARGS args) {
        return new CTYPE(args);
    }
  
}; // class ObjectProvider

        } // namespace Config
    } // namespace Utilities
} // namespace Mutation

#endif // AUTO_REGISTRATION_H


