
#include <iostream>
#include <algorithm>

#include "utilities.h"

using namespace std;

namespace utils {

void StringUtils::tokenize(
    const string &str, vector<string> &tokens, const string &delim,
    const bool multi_delim)
{
    if (str.empty()) return;
    
    string::size_type last_pos  = 0;
    string::size_type delim_pos;
    if (multi_delim)
        delim_pos = str.find_first_of(delim);
    else
        delim_pos = str.find(delim);
    
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
}

string& StringUtils::eraseAll(string &str, const string &to_erase)
{
    string::size_type pos = str.find(to_erase);
    
    while (pos != string::npos) {
        str.erase(pos, to_erase.length());
        pos = str.find(to_erase, pos);
    }
    
    return str;
}

string StringUtils::trim(const string &str, const string &to_erase)
{
    size_t i1 = str.find_first_not_of(to_erase);
    
    if (i1 != string::npos)
        return str.substr(i1, str.find_last_not_of(to_erase)-i1+1);
    else
        return string();
}

std::string StringUtils::removeWhiteSpace(const std::string& str)
{
    std::string::size_type pos = str.find_first_of(" \f\v\t\r\n");
    std::string to_return = str;
    
    while (pos != string::npos) {
        to_return.erase(pos, 1);
        pos = to_return.find_first_of(" \f\v\t\r\n", pos);
    }
    
    return to_return;
}

string StringUtils::toUpperCase(const string &str)
{
    std::string new_str;        
    std::transform(str.begin(), str.end(), std::back_inserter(new_str),
        static_cast<int(*)(int)>(toupper));
    return new_str;
}

string StringUtils::toLowerCase(const string &str)
{
    std::string new_str;        
    std::transform(str.begin(), str.end(), std::back_inserter(new_str),
        static_cast<int(*)(int)>(tolower));
    return new_str;
}

} // namespace utils


