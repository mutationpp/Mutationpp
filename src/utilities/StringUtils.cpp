/**
 * @file StringUtils.cpp
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

#include <algorithm>
#include "StringUtils.h"

using std::string;
using std::vector;

namespace Mutation {
    namespace Utilities {
        namespace String {

//==============================================================================

std::vector<std::string>& tokenize(
    const string &str, vector<string> &tokens, const string &delim,
    const bool multi_delim)
{
    if (str.empty()) return tokens;
    
    string::size_type last_pos  = 0;
    string::size_type delim_pos;
    if (multi_delim)
        delim_pos = str.find_first_of(delim);
    else
        delim_pos = str.find(delim);
    
    tokens.clear();
    
    while (last_pos != string::npos) {
        if (last_pos == delim_pos) {
            // Skip over delimeter(s)
            if (multi_delim) {
                last_pos++;
                delim_pos = str.find_first_of(delim, last_pos);
            } else {
                last_pos += delim.length();
                delim_pos = str.find(delim, last_pos);
            }
            
            // Special end case to prevent adding empty string to tokens
            if (last_pos >= str.length())
                last_pos = string::npos;
        } else {
            // Copy token
            if (delim_pos == string::npos)
                tokens.push_back(str.substr(last_pos));
            else
                tokens.push_back(str.substr(last_pos, delim_pos-last_pos));
            
            last_pos = delim_pos;
        }
    }

    return tokens;
}

//==============================================================================

string& eraseAll(string& str, const string& to_erase)
{
    string::size_type pos = str.find(to_erase);
    
    while (pos != string::npos) {
        str.erase(pos, to_erase.length());
        pos = str.find(to_erase, pos);
    }
    
    return str;
}

//==============================================================================

string trim(const string& str, const string& to_erase)
{
    size_t i1 = str.find_first_not_of(to_erase);
    
    if (i1 != string::npos)
        return str.substr(i1, str.find_last_not_of(to_erase)-i1+1);
    else
        return string();
}

//==============================================================================

std::string removeWhiteSpace(const string& str)
{
    string::size_type pos = str.find_first_of(" \f\v\t\r\n");
    string to_return = str;
    
    while (pos != string::npos) {
        to_return.erase(pos, 1);
        pos = to_return.find_first_of(" \f\v\t\r\n", pos);
    }
    
    return to_return;
}

//==============================================================================

string toUpperCase(const string &str)
{
    string new_str;        
    std::transform(str.begin(), str.end(), std::back_inserter(new_str),
        static_cast<int(*)(int)>(toupper));
    return new_str;
}

//==============================================================================

string toLowerCase(const string &str)
{
    string new_str;        
    std::transform(str.begin(), str.end(), std::back_inserter(new_str),
        static_cast<int(*)(int)>(tolower));
    return new_str;
}

//==============================================================================

bool isNumeric(const string& input, int number_base)
{
	const string base = "0123456789ABCDEF";
 
	return (
	    input.find_first_not_of(
	        base.substr(0, number_base) + ".+-Ee") == 
	    string::npos);
}

//==============================================================================

bool isNumeric(
    const vector<string>& array, int number_base)
{
    for (int i = 0; i < array.size(); ++i)
        if (!isNumeric(array[0], number_base)) return false;
    
    return true;
}

        } // namespace String
    } // namespace Utilities
} // namespace Mutation

