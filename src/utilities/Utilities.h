/**
 * @file Utilities.h
 *
 * @todo Packages all of the Mutation++ utilities together to reduce complexity
 * in other namespaces.
 */

#ifndef UTILITIES_UTILITIES_H
#define UTILITIES_UTILITIES_H

#include <string>
#include <cstdlib>

#include "AutoRegistration.h"
#include "IteratorWrapper.h"
#include "LookupTable.h"
#include "ReferenceServer.h"
#include "StringUtils.h"
#include "Units.h"
#include "XMLite.h"

// Define a useful debug output function which does nothing unless VERBOSE is
// defined
#ifndef VERBOSE
#define DEBUG(__out__) ;
#else
#define DEBUG(__out__) std::cout << __out__;
#endif

namespace Mutation {
    namespace Utilities {

/**
 *  Wraps the std::getenv() function to return a string.
 */
static std::string getEnvironmentVariable(const std::string& key)
{
    char* value = std::getenv(key.c_str());
    return (value == NULL ? std::string() : std::string(value));
}

    } // namespace Utilities
} // namespace Mutation

#endif // UTILITIES_UTILITIES_H
