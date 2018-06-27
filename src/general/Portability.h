/*
 * Copyright 2018 von Karman Institute for Fluid Dynamics (VKI)
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

#ifndef PORTABILITY_H
#define PORTABILITY_H

// Windows support
#if defined(WIN32) || defined(WNT)
#   ifdef MUTATION_EXPORTS
#       define MUTATION_API __declspec( dllexport )
#   else
#       define MUTATION_API __declspec( dllimport )
#   endif
#   define NAME_MANGLE(__name__)  MUTATION_API  mpp_##__name__##_
#   pragma warning(disable: 4251) // 'identifier' : class 'type' needs to have dll-interface to be used by clients of class 'type2'

// Other OS
#else
#   include <alloca.h>

#   define MUTATION_API
#   define NAME_MANGLE(__name__) mpp_##__name__##_
#endif

#endif // PORTABILITY_H
