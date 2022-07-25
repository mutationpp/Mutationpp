/**
 * @file IteratorWrapper.h
 *
 * @brief Provides IteratorWrapper class.
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

#ifndef ITERATOR_BASE_H
#define ITERATOR_BASE_H

namespace Mutation {
    namespace Utilities {

/**
 * Simple class that acts to wrap a predefined iterator class.  
 * 
 * This is useful if you would like to provide an iterator to some container 
 * inside your class but you don't want to have to provide the container's 
 * iterator directly in order to maintain data encapsulation inside your class.
 * Using IteratorWrapper class allows you to provide an iterator for your class
 * without having to reimplement the container's iterator but also not making
 * a user of your class dependent on which type of container is being used.  For
 * example if you had a class defined as
 *
 * @code
 * class MyClass
 * {
 * public:
 *
 *     ...
 *
 *     class Iterator 
 *         : public IteratorWrapper<double, std::vector<double>::iterator>
 *     {
 *         friend class MyClass;
 *     }
 *     
 *     Iterator &begin()
 *     {
 *         m_iterator.m_iter = m_container.begin();
 *         return m_iterator;
 *     }
 *
 *     Iterator &end()
 *     {
 *         m_iterator.m_iter = m_container.end();
 *         return m_iterator;
 *     }
 *
 * private:
 *     
 *     ...
 *
 *     std::vector<double> m_container;
 *     Iterator            m_iterator;
 * } // class MyClass
 * @endcode
 *
 * You could then use the iterator as follows
 *
 * @code
 * MyClass myclass(...);
 *
 * MyClass::Iterator iter = myclass.begin();
 * Myclass::Iterator end  = myclass.end();
 *
 * for ( ; iter != end; ++iter)
 *     // do something with iter
 * @endcode
 */
template <typename T, typename Iterator>
class IteratorWrapper
{
public:
    // assignment operator
    IteratorWrapper &operator=(const IteratorWrapper &to_copy)
    {
        if (*this != to_copy)
            m_iter = to_copy.m_iter;
        return *this;
    }
    
    // == and != operators
    bool operator==(const IteratorWrapper &to_check) {
        return (m_iter == to_check.m_iter);
    }        
    bool operator!=(const IteratorWrapper &to_check) {
        return !(*this == to_check);
    }
    
    // ++ and -- operators
    IteratorWrapper& operator++()
    {
        m_iter++;
        return *this;
    }
    IteratorWrapper& operator--()
    {
        m_iter--;
        return *this;
    }
    
    // * and -> operators
    T &operator*() {
        return *m_iter;
    }
    T *operator->() {
        return m_iter.operator->();
    }
    
protected:
    Iterator m_iter;
};

    } // namespace Utilities
} // namespace Mutation

#endif // ITERATOR_BASE_H
