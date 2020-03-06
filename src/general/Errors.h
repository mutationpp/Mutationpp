/*
 * Copyright 2017-2020 von Karman Institute for Fluid Dynamics (VKI)
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

//#include <algorithm>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace Mutation {

void terminateOnError();

/**
 * Base class for all Mutation++ exceptions.  Defines a consistent look-and-feel
 * for all error messages in the library.  Also ensures that catching the Error
 * type will catch all exceptions thrown by Mutation++ code.  Finally, the last
 * Error (or derived) object which is constructed before std::terminate is
 * called, is guaranteed to be printed to standard output.  This ensures that a
 * user of M++ will be notified of an uncaught Error object.
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
    Error(const std::string& type) :
        m_type(type),
        m_extra_info(),
        m_formatted_message()
    {
        // Set this Error as the most recent
        lastError() = this;
        errorCount()++;

        // The first Error object created needs to save the current terminate
        // handler
        if (errorCount() == 1)
            terminateHandler() = std::set_terminate(terminateOnError);
        
        formatMessage();
    }

    /**
     * Copy constructor.
     *
     * @param copy The object to copy.
     */
    Error(const Error& copy) :
        m_type(copy.m_type),
        m_extra_info(copy.m_extra_info),
        m_formatted_message(copy.m_formatted_message)
    {
        m_message_stream << copy.m_message_stream.rdbuf();

        // Set this Error as the most recent
        lastError() = this;
        errorCount()++;
    }

    /// Destructor.
    virtual ~Error() throw()
    {
        // Reset the terminate function
        errorCount()--;
        if (errorCount() == 0)
            std::set_terminate(terminateHandler());
    }

    /// Input stream operator allows for message to be added to error.
    template <typename T>
    Error& operator<< (const T& in)
    {
        m_message_stream << in;
        formatMessage();
        return *this;
    }

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
    virtual const char* what() const throw() {
        return m_formatted_message.c_str();
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
            formatMessage();
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

    // Give private access to terminate function
    friend void terminateOnError();

private:

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
    void formatMessage()
    {
        m_formatted_message = "\nM++ error: " + m_type + ".\n";
        
        std::vector<std::pair<std::string, std::string> >::iterator it;
        for (it = m_extra_info.begin(); it != m_extra_info.end(); ++it)
            m_formatted_message += it->first + ": " + it->second + "\n";
        
        m_formatted_message += m_message_stream.str() + "\n";
    }

    /// Keeps track of the most recent Error.
    static Error*& lastError() {
        static Error* p_last_error = NULL;
        return p_last_error;
    }

    /// Keeps track of the number of instantiated Error objects.
    static int& errorCount() {
        static int count = 0;
        return count;
    }

    /// Keeps track of the default terminate handler.
    static std::terminate_handler& terminateHandler() {
        static std::terminate_handler handler;
        return handler;
    }

private:

    /// Type of error message.
    std::string m_type;

    /// The error message.
    std::stringstream m_message_stream;

    /// Vector of additional (name, value) pairs
    std::vector<std::pair<std::string, std::string> > m_extra_info;

    /// The fully formatted message which gets displayed by what()
    std::string m_formatted_message;

}; // class Error


/**
 * Reports the last uncaught Error.
 */
inline void terminateOnError()
{
    // Notify the user of the last error.
    std::cout << Error::lastError()->what() << std::endl;

    // Try to clean up before terminating
    std::exit(1);
}


/**
 * This class bridges between the basic Error class and subclasses of Error and
 * enables the correct polymorphism behavior while also providing function
 * overloading for methods which need to return a reference to the Error object.
 * Without these overloaded functions, the code
 * @code{c++}
 * throw DerivedError() << "some message";
 * @endcode
 * would actually throw an Error object, rather than the intended DerivedError
 * type.  This is not really a problem if you just want to catch all Error
 * types, however if you need to catch a specific derived type, this would not
 * work. This class simply provides a convenient alternative to having to write
 * the function overloads specifically in each subclass of Error.  All derived
 * exception types should extend from this type.
 */
template <typename ConcreteType>
class ErrorExtension : public Error
{
public:

    /// Destructor.
    virtual ~ErrorExtension() throw() { }

    /// Overloaded << operator.
    template <typename T>
    ConcreteType& operator<< (const T& in)
    {
        Error::operator<<(in);
        return static_cast<ConcreteType&>(*this);
    }

    /// Overloaded () operator.
    template <typename T>
    ConcreteType& operator()(const std::string& name, const T& value) {
        Error::addExtraInfo(name, value);
        return static_cast<ConcreteType&>(*this);
    }

protected:

    /// This class is only meant to be extended by concrete Error types.
    ErrorExtension(const std::string& type) :
        Error(type)
    { }
};


/// Reports a "file not found" error.  Includes "file" info.
class FileNotFoundError : public ErrorExtension<FileNotFoundError>
{
public:
    /// Constructor taking the file name.
    FileNotFoundError(const std::string& file) :
        ErrorExtension<FileNotFoundError>("file not found")
    {
        addExtraInfo("file", file);
    }

    /// Returns the file name.
    const std::string& file() const {
        return getExtraInfo("file");
    }
}; // class FileNotFoundError


/// Reports an "error parsing file" error. Includes "file" and "line" info.
class FileParseError : public ErrorExtension<FileParseError>
{
public:
    /// Constructor taking the file and line number.
    FileParseError(const std::string& file, int line) :
        ErrorExtension<FileParseError>("error parsing file")
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


/// Reports an "invalid input" error.  Includes input name and value.
class InvalidInputError : public ErrorExtension<InvalidInputError>
{
public:
    /// Constructor taking the input name and value.
    template <typename T>
    InvalidInputError(const std::string& name, const T& value) :
        ErrorExtension<InvalidInputError>("invalid input"),
        m_name(name)
    {
        std::stringstream ss; ss << value;
        m_value = ss.str();
        addExtraInfo("input", m_name + " = " + m_value);
    }

    /// Destructor.
    ~InvalidInputError() throw() { }

    /// Returns the name of the invalid input.
    const std::string& inputName() const {
        return m_name;
    }

    /// Returnst the value of the invalid input.
    const std::string& inputValue() const {
        return m_value;
    }
private:
    std::string m_name;
    std::string m_value;
}; // class InvalidInput


/**
 * Reports a "logic error".  Includes file and line number.  This exception type
 * should only be used in cases which fail due to programming logic errors which
 * should be preventable and cannot be handled with simple assert() statements.
 */
class LogicError : public ErrorExtension<LogicError>
{
public:
    LogicError(const char* file, int line) :
        ErrorExtension<LogicError>("logic error")
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
}; // class LogicError

#define LogicError() LogicError(__FILE__, __LINE__)


/// Reports a "function not implemented" error.  Includes file, line number, and
/// function name.
class NotImplementedError : public ErrorExtension<NotImplementedError>
{
public:
    NotImplementedError(const char* func, const char* file, int line) :
        ErrorExtension<NotImplementedError>("function not implemented")
    {
        addExtraInfo("function", func);
        addExtraInfo("file", file);
        addExtraInfo("line", line);
    }

    /// Returns the file name.
    const std::string& file() const {
        return getExtraInfo("file");
    }

    /// Returns the function.
    const std::string& function() const {
        return getExtraInfo("function");
    }

    /// Returns the line number (as a string).
    const std::string& line() const {
        return getExtraInfo("line");
    }
};

/// Reports a "missing data" error.  
class MissingDataError : public ErrorExtension<MissingDataError>
{
public:
    MissingDataError() :
        ErrorExtension<MissingDataError>("missing data")
    { }
};

#define NotImplementedError(_func_) NotImplementedError(_func_, __FILE__, __LINE__)

} // namespace Mutation

#endif // ERRORS_H
