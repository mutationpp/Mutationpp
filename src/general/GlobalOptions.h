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

#ifndef GLOBAL_OPTIONS_H
#define GLOBAL_OPTIONS_H

namespace Mutation {

/// Singleton class providing access to global options.
class GlobalOptions
{
public:

    /// Access to singleton GlobalOptions object.
    static GlobalOptions& getInstance() {
        static GlobalOptions opts;
        return opts;
    }

    /// Resets all options to the default values.
    static void reset() {
        getInstance().resetOptions();
    }

    /// Gets the data directory.
    static const std::string& dataDirectory() {
        return getInstance().m_data_directory;
    }

    /// Sets the data directory.
    static void dataDirectory(const std::string& dir) {
        getInstance().m_data_directory = dir;
    }

    /// Gets the data directory.
    static const std::string& workingDirectory() {
        return getInstance().m_working_directory;
    }

    /// Sets the data directory.
    static void workingDirectory(const std::string& dir) {
        getInstance().m_working_directory = dir;
    }

    /// Gets the file separator character.
    static const char separator() {
        return getInstance().m_separator;
    }

    /// Sets the file separator character.
    static void separator(char sep) {
        getInstance().m_separator = sep;
    }

private:

    /// Sets the default options and makes class a singleton.
    GlobalOptions() { resetOptions(); }

    /// Resets all options to the default values.
    void resetOptions() {
        m_data_directory = getEnvironmentVariable("MPP_DATA_DIRECTORY");
        m_working_directory = "";
#ifdef _WIN32
        m_separator = '\\';
#else
        m_separator = '/';
#endif
    }

    /**
     * Returns the value for the environment variable key.  If the key is not
     * found, then an empty string is returned.
     */
    std::string getEnvironmentVariable(const std::string& key)
    {
        char* value = std::getenv(key.c_str());
        return std::string(value == NULL ? "" : value);
    }

private:

    /// Data directory
    std::string m_data_directory;

    /// Working directory
    std::string m_working_directory;

    /// File separator character
    char m_separator;

}; // class GlobalOptions

} // namespace Mutation

#endif // GLOBAL_OPTIONS_H
