/**
 * @file StringUtils.h
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

#ifndef UTILITIES_STRING_UTILS_H
#define UTILITIES_STRING_UTILS_H

#include <string>
#include <vector>

namespace Mutation {
    namespace Utilities {
        namespace String {

/**
 * Tokenizes a string into a vector of sub strings that are separated by the 
 * characters in teh delim string.
 */
std::vector<std::string>& tokenize(
    const std::string &str, std::vector<std::string> &tokens, 
    const std::string &delim = " ", const bool multi_delim = true);


/**
 * Removes all characters from the string belonging to to_erase.
 */
std::string& eraseAll(
    std::string &str, const std::string &to_erase = " \n\t");

/**
 * Returns the string with the whitespace at the beginning and end of the string
 * removed.
 */
std::string trim(
    const std::string& str, const std::string &to_erase = " \t\f\v\n\r");
    
/**
 * Returns the string with all characters belonging to the set {' ', '\t', '\f',
 * '\v', '\r', '\n'} removed.
 */
std::string removeWhiteSpace(const std::string& str);

/**
 * Returns the string in upper case.
 */
std::string toUpperCase(const std::string& str);

/**
 * Returns the string in lower case.
 */
std::string toLowerCase(const std::string& str);

/**
 * Returns true if the string represents a number in ascii where number_base is
 * the base of the number space the string represents.
 */
bool isNumeric(const std::string& input, int number_base = 10);

/**
 * Returns true if all the strings in the vector represent numbers in ascii.
 * @see isNumeric()
 */
bool isNumeric(
    const std::vector<std::string>& array, int number_base = 10);


        } // namespace String
    } // namespace Utilities
} // namespace Mutation

#endif // UTILITIES_STRING_UTILS_H
