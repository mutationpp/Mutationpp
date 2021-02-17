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

#ifndef TEMPORARY_FILE_H
#define TEMPORARY_FILE_H

#include <fstream>
#include <string>

namespace Mutation {
namespace Utilities {
namespace IO {

/**
 * Returns a unique filename.
 */
static std::string uniqueFileName();

/**
 * Creates a temporary file with a random filename which is automatically
 * deleted when the object leaves scope.
 */
class TemporaryFile
{
public:

    /// Opens a temporary filename which can be read and written to.
    TemporaryFile(const char* const extension = "");

    /// Close and delete the temporary file.
    virtual ~TemporaryFile();

    /// Returns the name of the temporary file.
    const std::string& filename() const { return m_filename; }

    /// Closes the file.
    void close() { m_fstream.close(); }

    /**
     * @brief Sets whether or not to delete the temporary file when the object
     * is destroyed.
     * 
     * This can be useful for testing.
     * 
     * @param del Defaults to true.
     */
    void deleteOnDestruct(bool del) { m_delete = del; }

    /// Reopens the file
    void open();

    /// Provides access to output stream operator of underlying fstream object.
    template <typename T>
    std::ostream& operator<< (const T& out) {
        return (m_fstream << out);
    }

    /// Provides access to input stream operator of underlying fstream object.
    template <typename T>
    std::istream& operator>> (const T& in) {
        return (m_fstream >> in);
    }

private:

    std::string m_filename;
    std::fstream m_fstream;
    bool m_delete;

}; // class TemporaryFile

} // namespace IO
} // namespace Utilities
} // namespace Mutation

#endif // TEMPORARY_FILE_H
