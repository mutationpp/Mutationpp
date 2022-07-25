/**
 * @file Utilities.h
 *
 * @brief Packages all of the Mutation++ utilities together to reduce complexity
 * in other namespaces.
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

#ifndef UTILITIES_UTILITIES_H
#define UTILITIES_UTILITIES_H

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "AutoRegistration.h"
#include "GlobalOptions.h"
#include "IteratorWrapper.h"
#include "LookupTable.h"
#include "ReferenceServer.h"
#include "StringUtils.h"
#include "TemporaryFile.h"
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
 * Prepends a directory path to a file name.
 */
static std::string prependPath(
    const std::string& dir, const std::string& file)
{
    // For an empty directory string, return the file name only
    if (dir.length() == 0)
        return file;

    // Check if the separator is already appended on the directory name
    const char sep = GlobalOptions::separator();
    if (dir[dir.length()-1] == sep)
        return dir + file;
    else
        return dir + sep + file;
}

/**
 * Appends a file extension to the given file name.  If the extension is empty,
 * no change is made.  If the period character is not already included in the
 * given extension, then it is added as well.
 */
static std::string appendExtension(
     const std::string& name, const std::string& ext)
{
    // For an empty extension, return the file name as is
    if (ext.length() == 0)
        return name;

    // Check if period character is already included
    if (ext[0] == '.')
        return name + ext;
    else
        return name + '.' + ext;
}

/**
 * Finds the full name to use for a database file name.  If the given extension
 * is not included, then it is appended to the file name.  The file path is
 * located in the following order:
 *
 * 1) working_directory/name.ext
 * 2) working_directory/dir/name.ext
 * 3) data_directory/name.ext
 * 4) data_directory/dir/name.ext
 *
 * If the file is not found, then the last location searched is returned.
 */
static std::string databaseFileName(
    std::string name, const std::string& dir, const std::string& ext = ".xml")
{
    // Add extension if necessary
    if (name.length() < ext.length() + 1 ||
        name.substr(name.length()-ext.length()) != ext)
        name = appendExtension(name, ext);

    // Check the working directory first
    std::string path = prependPath(GlobalOptions::workingDirectory(), name);
    if (std::ifstream(path.c_str(), std::ios::in).is_open())
        return path;

    // Check the working directory with assumed folder hierarchy
    path = prependPath(GlobalOptions::workingDirectory(), dir);
    path = prependPath(path, name);
    if (std::ifstream(path.c_str(), std::ios::in).is_open())
        return path;

    // Check data directory
    path = prependPath(GlobalOptions::dataDirectory(), name);
    if (std::ifstream(path.c_str(), std::ios::in).is_open())
        return path;

    // If we didn't find it there then just return path with data directory and
    // assumed folder hierarchy
    path = prependPath(GlobalOptions::dataDirectory(), dir);
    return prependPath(path, name);
}

/**
 * Returns the ordinal suffix of a given integer.
 */
static const char* ordinalSuffix(int n)
{
    switch (n%100) {
    case 11: case 12: case 13: return "th";
    default:
        switch (n%10) {
        case 1: return "st";
        case 2: return "nd";
        case 3: return "rd";
        default:
            return "th";
        }
    }
}


    } // namespace Utilities
} // namespace Mutation

#endif // UTILITIES_UTILITIES_H
