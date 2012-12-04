/**
 * @file utilities.h
 *
 * @todo Document this file!
 */

#ifndef UTILITIES_UTILITIES_H
#define UTILITIES_UTILITIES_H

#include <string>
#include <vector>
#include <cstdlib>

namespace utils {

static std::string getEnvironmentVariable(const std::string& key)
{
    char* value = std::getenv(key.c_str());
    return (value == NULL ? std::string() : std::string(value));
}

static bool isNumeric(const std::string& input, int number_base = 10)
{
	const std::string base = "0123456789ABCDEF";
 
	return (
	    input.find_first_not_of(
	        base.substr(0, number_base) + ".+-Ee") == 
	    std::string::npos);
}

static bool isNumeric(
    const std::vector<std::string>& array, int number_base = 10)
{
    for (int i = 0; i < array.size(); ++i)
        if (!isNumeric(array[0], number_base)) return false;
    
    return true;
}

class StringUtils
{
public:

    static void tokenize(
        const std::string &str, std::vector<std::string> &tokens, 
        const std::string &delim = " ", const bool multi_delim = true);
    
    static std::string &eraseAll(
        std::string &str, const std::string &to_erase = " \n\t");
    
    static std::string trim(
        const std::string& str, const std::string &to_erase = " \t\f\v\n\r");
    
    /**
     * Removes all characters belonging to the set {' ', '\t', '\f', '\v', '\r',
     * '\n'}.
     */
    static std::string removeWhiteSpace(const std::string& str);
    
    static std::string toUpperCase(const std::string& str);
    
    static std::string toLowerCase(const std::string& str);
    
}; // class StringUtils

} // namespace utilities

#endif // UTILITIES_UTILITIES_H
