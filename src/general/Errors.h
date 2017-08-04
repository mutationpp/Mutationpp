/*
 * Copyright 2017 von Karman Institute for Fluid Dynamics (VKI)
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

/**
 * @file Errors.h
 * @brief Implements error handling for Mutation++.
 */

#ifndef ERRORS_H
#define ERRORS_H

#include <algorithm>
#include <exception>
#include <string>
#include <vector>
#include <sstream>

namespace Mutation {

/**
 * Base class for all Mutation exceptions.  Defines a consistent look-and-feel
 * for all error messages in the library.
 */
class Error : public std::exception
{
public:

    /**
     * Creates a new Error object with a given type string and error message.
     *
     * @param type A string representing the type of error.
     * @param message The error message.
     */
    Error(const std::string& type, const std::string& message) :
        m_type(type), m_message(message)
    { }

    /// Destructor.
    virtual ~Error() throw() { }

    /**
     * Creates a formatted error message.  The message is formatted as
     *
     * @code
     * M++ error: <type>.
     *   <name>: <value>
     *   <message>
     * @endcode
     *
     * where <type> and <message> correspond to the arguments in the constructor
     * and the optional (name, value) lines correspond to extra information
     * provided by addExtraInfo() or operator().
     */
    const char* what() throw()
    {
        try {
            std::stringstream ss;
            std::vector<std::pair<std::string, std::string> >::iterator it;

            // Format the error message
            ss << "M++ error: " << m_type << ".\n  ";
            for (it = m_extra_info.begin(); it != m_extra_info.end(); ++it)
                ss << it->first << ": " << it->second << "\n  ";
            ss << m_message << "\n";

            // Save the formatted message (so the pointer actually points to
            // something when the method returns!)
            m_formatted_message = ss.str();
            return m_formatted_message.c_str();

        } catch (std::exception&) {
            // In the unlikely event that the formatting actually throws an
            // error, we need to have a backup plan.  Return unformatted
            // message.  This is guaranteed not to throw an exception according
            // to the C++ standard library docs.
            return m_message.c_str();
        }
    }

    /**
     * Adds a (name, value) pair as an additional piece of information for the
     * error.  The name and value pair will be printed on an individual line
     * after the error type and before the error message.  The extra info lines
     * are printed in the order they were added to the Error object.
     *
     * @param name The name of the extra piece of info.
     * @param value The value of the extra piece of info (must be convertible to
     * a std::string using std::stringstream).
     */
    template <typename T>
    void addExtraInfo(const std::string& name, const T& value)
    {
        // This should be safe, any type T that isn't supported by the <<
        // operator should cause a compilation error.
        try {
            std::stringstream ss;
            ss << value;
            m_extra_info.push_back(make_pair(name, ss.str()));
        } catch (std::exception&) {
            // In the unlikely event this throws an exception, just ignore the
            // extra info.  It is better to protect this exception and get some
            // error message, then to get a message about what happened here.
        }
    }

    /// Provides a shorthand notation for addExtraInfo().
    template <typename T>
    Error& operator()(const std::string& name, const T& value) {
        addExtraInfo(name, value);
        return *this;
    }

    /**
     * Returns the value of the extra info with given name.
     *
     * @param name The name of the desired info.
     * @return The value of the info if it exists, otherwise an empty string
     */
    const std::string& getExtraInfo(const std::string& name) const
    {
        // This should be safe from exceptions.
        std::vector<std::pair<std::string, std::string> >::const_iterator it;
        for (it = m_extra_info.begin(); it != m_extra_info.end(); ++it)
            if (it->first == name) return it->second;

        // Could throw, but there is nothing we could really do because we have
        // to return something.  It is very very unlikely this would throw.
        static std::string empty;
        return empty;
    }

private:

    /// Type of error message.
    std::string m_type;

    /// The error message.
    std::string m_message;

    /// Vector of additional (name, value) pairs
    std::vector<std::pair<std::string, std::string> > m_extra_info;

    /// The fully formatted message which gets displayed by what()
    std::string m_formatted_message;

}; // class Error


/// Reports a "file not found" error.  Includes "file" info.
class FileNotFoundError : public Error
{
public:
    /// Constructor taking the file name and error message.
    FileNotFoundError(const std::string& file, const std::string& m) :
        Error("file not found", m)
    {
        addExtraInfo("file", file);
    }

    /// Returns the file name.
    const std::string& file() const {
        return getExtraInfo("file");
    }
}; // class FileNotFoundError


/// Reports an "error parsing file" error. Includes "file" and "line" info.
class FileParseError : public Error
{
public:
    /// Constructor taking the file, line number, and error message.
    FileParseError(const std::string& file, int line, const std::string& m) :
        Error("error parsing file", m)
    {
        addExtraInfo("file", file);
        addExtraInfo("line", line);
    }

    /// Returns the file name.
    const std::string& file() const {
        return getExtraInfo("file");
    }

    /// Returns the line number (as a string).
    const std::string& line() const {
        return getExtraInfo("line");
    }
}; // class FileParseError

} // namespace Mutation

#endif // ERRORS_H
