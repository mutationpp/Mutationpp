/**
 * @file AutoRegistration.h
 *
 * Provides classes which implement automatic registration of polymorphic types
 * from a name which identifies a subtype.  The code in this file has been 
 * entirely taken from the SConfig library written by Andrea Lani and the
 * documentation and formatting has only been altered to match the standards
 * used in the Mutation++ library.  In addition, the use of exceptions has been
 * changed to match the Mutation++ exception policy.
 *
 * @author Andrea Lani
 * @author J.B. Scoggins
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef AUTO_REGISTRATION_H
#define AUTO_REGISTRATION_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Errors.h"

namespace Mutation {
    namespace Utilities {
        namespace Config {

template <class PTYPE> class Provider;

/**
 * This class stores an archive of available Provider's for a generic 
 * polymorphic type PTYPE. Once registered, they can be accessed by name.
 * This class is originally taken and modified from the VKI CoolFLUID package.
 */
template <class PTYPE>
class Factory 
{
public:
  
    /**
     * Register a Provider.
     */
    void add(Provider<PTYPE>* p_provider)
    {
        typename std::map<std::string, Provider<PTYPE>*>::const_iterator it;
        it = m_providers.find(p_provider->getName());

        if (it == m_providers.end())
            m_providers[p_provider->getName()] = p_provider;
        else {
            throw LogicError()
                << "Provider <" << p_provider->getName() << "> has already "
                << "been registered for type " << PTYPE::typeName();
        }
    }
  
    /**
     * Get the provider corresponding to the given key.
     */
    Provider<PTYPE>* getProvider(const std::string& key)
    {    
        typename std::map<std::string, Provider<PTYPE>*>::const_iterator it;
        it = m_providers.find(key);

        if (it != m_providers.end())
            return dynamic_cast<Provider<PTYPE>*>(it->second);
        
        // Error because key is not a valid provider
        InvalidInputError error("key", key);
        error << "Provider <" <<  key << "> for type " << PTYPE::typeName()
              << " is not registered.  Possible providers are:\n";
        for (it = m_providers.begin() ; it != m_providers.end(); it++)
            error << "  " << it->first << "\n";
        throw error;
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
        Provider<PTYPE>* p_provider = getInstance().getProvider(name);
        if (p_provider == NULL)
        	return NULL;
        return p_provider->create(args);
    }

    /**
     * Returns a list of all object names registered to this factory.
     */
    static std::vector<std::string> names() {
        std::vector<std::string> names;
        typename std::map<std::string, Provider<PTYPE>*>& providers =
            getInstance().m_providers;
        typename std::map<std::string, Provider<PTYPE>*>::const_iterator iter;
        for (iter = providers.begin(); iter != providers.end(); ++iter)
            names.push_back(iter->first);
        return names;
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


