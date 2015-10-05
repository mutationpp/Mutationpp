/**
 * @file Utilities.h
 *
 * @brief Packages all of the Mutation++ utilities together to reduce complexity
 * in other namespaces.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
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


#ifndef UTILITIES_UTILITIES_H
#define UTILITIES_UTILITIES_H

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>

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
    if (value == NULL) {
        std::cout << "Warning: environment variable " << key << " not found" << std::endl;
        return std::string();
    }
    return std::string(value);
}

/**
 * Finds the full name to use for a database file name.
 */
static std::string databaseFileName(
    std::string name, const std::string& dir, const std::string& ext = ".xml")
{
    // Add extension if necessary
    if (name.substr(name.length()-ext.length()) != ext)
        name = name + ext;

    // Add directory if not in local
    {
        std::ifstream file(name.c_str(), std::ios::in);

        // If that doesn't work, look in MPP_DATA_DIRECTORY/transport
        if (!file.is_open())
            name = Utilities::getEnvironmentVariable("MPP_DATA_DIRECTORY") +
                "/" + dir + "/" + name;
    }

    return name;
}

    } // namespace Utilities
} // namespace Mutation

#endif // UTILITIES_UTILITIES_H
