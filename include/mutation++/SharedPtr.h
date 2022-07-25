/**
 * @file SharedPtr.h
 *
 * Provides the SharedPtr type.
 *
 * @author J.B. Scoggins
 */

/*
 * Copyright 2015-2020 von Karman Institute for Fluid Dynamics (VKI)
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

#ifndef UTILITIES_SHARED_PTR_H
#define UTILITIES_SHARED_PTR_H

// Use std::shared_ptr implementation available in C++11 if supported
#if __cplusplus > 199711L
    //#pragma message("using std::shared_ptr")
    #include <memory>
    namespace Mutation {
        template <typename T> using SharedPtr = std::shared_ptr<T>;
    }

// Otherwise, use our own implementation (probably slower)
#else

#include <algorithm> // for std::swap

namespace Mutation {

/**
 * Very simple reference counting smart pointer.
 */
template <typename T>
class SharedPtr
{
public:
    SharedPtr() : mp_ref(NULL) { }

    /**
     * Constructs a new SharedPtr owning the object pointed to by the pointer.
     */
    explicit SharedPtr(T* ptr) :
        mp_ref(new RefCount())
    {
        mp_ref->ptr = ptr;
        mp_ref->count = 1;
    }

    /**
     * Copy constructor.
     */
    SharedPtr(const SharedPtr& other) :
        mp_ref(other.mp_ref)
    {
        mp_ref->count++;
    }

    /**
     * Decreases reference count to the object.  If this is the last reference,
     * then the object is deleted.
     */
    ~SharedPtr()
    {
        // Protect against non-initialized pointers
        if (mp_ref == NULL) return;

        // Reduce reference count and delete object
        if (--(mp_ref->count) == 0) {
            delete mp_ref->ptr;
            delete mp_ref;
        }
    }

    /// Assignment operator.
    SharedPtr<T>& operator = (SharedPtr<T> other) {
        std::swap(mp_ref, other.mp_ref);
        return *this;
    }

    /// Equality operator.
    bool operator == (SharedPtr<T> other) {
        return (mp_ref == other.mp_ref);
    }

    /// Returns reference to object owned by this shared pointer.
    T& operator*() const { return *(mp_ref->ptr); };

    /// Returns pointer to object owned by this shared pointer.
    T* operator->() const { return mp_ref->ptr; };

    template <class U, class V>
    friend bool operator == (const SharedPtr<U>&, const SharedPtr<V>&);

    template <class U, class V>
    friend bool operator != (const SharedPtr<U>&, const SharedPtr<V>&);

private:

    // Helper class for managing the reference count.
    struct RefCount
    {
        T*  ptr;
        int count;
    };

private:

    RefCount* mp_ref;

}; // class SharedPtr

// Compare two SharedPtr objects.
template <class T, class U>
bool operator == (const SharedPtr<T>& lhs, const SharedPtr<U>& rhs) {
    return (lhs.mp_ref == rhs.mp_ref);
}
template <class T, class U>
bool operator != (const SharedPtr<T>& lhs, const SharedPtr<U>& rhs) {
    return (lhs.mp_ref != rhs.mp_ref);
}



} // namespace Mutation

#endif // __cplusplus > 199711L
#endif // UTILITIES_SHARED_PTR_H

