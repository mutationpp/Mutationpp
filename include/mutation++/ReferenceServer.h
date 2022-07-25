/**
 * @file ReferenceServer.h
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

#ifndef REFERENCE_SERVER_H
#define REFERENCE_SERVER_H

namespace Mutation {
    namespace Utilities {

/**
 * Non-instantiable class which serves up revolving references to objects of 
 * class T.  At start up, an array of size elements is created using the 
 * default constructor of type T and an index initially points to the first
 * object in the array.  The serve() method returns a reference to the object
 * which is next in the list and updates the index counter.  When the counter
 * hits the end of the list, it is reset to zero.  This class is useful when
 * you want to be able to access temporary objects of some type without having
 * to instantiate them over and over.  The size template parameter should be set
 * to the maximum number of objects that could be required simultaneously.
 */
template <typename T, int size>
class ReferenceServer
{
public:
    static T& serve() {
        return m_refs[(++m_current < size ? m_current : m_current = 0)];
    }

private:
    static T m_refs [size];
    static int m_current;
}; // class ReferenceServer

template <typename T, int size>
int ReferenceServer<T, size>::m_current = 0;

template <typename T, int size>
T ReferenceServer<T, size>::m_refs [size];

    } // namespace Utilities
} // namespace Mutation

#endif // REFERENCE_SERVER_H

