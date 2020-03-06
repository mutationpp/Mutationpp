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

#include "Errors.h"
#include "TemporaryFile.h"
#include <ctime>

namespace Mutation {
namespace Utilities {
namespace IO {

//==============================================================================

static std::string uniqueFileName(const char* const extension)
{
    // List of characters to use in the file name
    const std::string char_set =
        "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

    static unsigned short lfsr = std::clock();
    unsigned short bit;
    std::string name;

    do {
        name.clear();

        // Use simple linear feedback shift register to generate name
        for (int i = 0; i < 10; ++i) {
            bit = ((lfsr >> 0) ^ (lfsr >> 2) ^ (lfsr >> 3) ^ (lfsr >> 5)) & 1;
            lfsr =  (lfsr >> 1) | (bit << 15);
            name += char_set[lfsr % char_set.length()];
        }
        name += extension;
    } while (std::ifstream(name.c_str()).is_open());

    return name;
}

//==============================================================================

TemporaryFile::TemporaryFile(const char* const extension) :
    m_filename(uniqueFileName(extension)), m_delete(true)
{
    // Create the file
    if (!std::ofstream(m_filename.c_str()).is_open())
        throw Error("cannot create file")
            << "Trying to create temporary file \"" << m_filename << "\".";

    // Open with read/write access
    open();
}

//==============================================================================

TemporaryFile::~TemporaryFile()
{
    // Close and delete the file
    if (m_fstream.is_open()) m_fstream.close();
    if (m_delete) std::remove(m_filename.c_str());
}

//==============================================================================

void TemporaryFile::open()
{
    if (!m_fstream.is_open())
        m_fstream.open(m_filename.c_str());

    if (!m_fstream.is_open())
        throw Error("cannot open file")
            << "Trying to open temporary file \"" << m_filename << "\".";
}

} // namespace IO
} // namespace Utilities
} // namespace Mutation

