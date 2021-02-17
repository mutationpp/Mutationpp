###############################################################################
#
# High level cmake configuration file for the Mutation++ project and associated
# helper tools and binaries.  All of the various configuration options are
# defined here.
#
# author: J.B. Scoggins (scoggins@vki.ac.be)
#
###############################################################################
#
# Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
#
# This file is part of MUlticomponent Thermodynamic And Transport
# properties for IONized gases in C++ (Mutation++) software package.
#
# Mutation++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Mutation++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with Mutation++.  If not, see
# <http://www.gnu.org/licenses/>.
#
cmake_minimum_required(VERSION 3.5)
cmake_policy(SET CMP0048 NEW)
set(CMAKE_CXX_STANDARD 11)

project(mutation++
    VERSION 1.0.1
    LANGUAGES CXX
)


# Add include path for cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Setup languages
enable_language(CXX)
if (BUILD_FORTRAN_WRAPPER)
    enable_language(Fortran)
endif()


###############################################################################
# Build types
###############################################################################

# Profile
set (CMAKE_CXX_FLAGS_PROFILE "-g3 -Wall -O3 -DNDEBUG" CACHE STRING
    "Flags used by the C++ compiler during Profile builds."
    FORCE )
set (CMAKE_C_FLAGS_PROFILE "-g3 -Wall -pedantic -O3 -DNDEBUG" CACHE STRING
    "Flags used by the C compiler during Profile builds."
    FORCE )
set (CMAKE_EXE_LINKER_FLAGS_PROFILE
    "" CACHE STRING
    "Flags used for linking binaries during Profile builds."
    FORCE )
set (CMAKE_SHARED_LINKER_FLAGS_PROFILE
    "" CACHE STRING
    "Flags used by the shared libraries linker during Profile builds."
    FORCE )
mark_as_advanced(
    CMAKE_CXX_FLAGS_PROFILE
    CMAKE_C_FLAGS_PROFILE
    CMAKE_EXE_LINKER_FLAGS_PROFILE
    CMAKE_SHARED_LINKER_FLAGS_PROFILE)


# Update the documentation string of CMAKE_BUILD_TYPE for GUIs and set the
# default build type as Release
if (NOT CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: None Debug Release
        RelWithDebInfo MinSizeRel Profile Coverage." FORCE)
endif ()


###############################################################################
# Install prefix settings
###############################################################################
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH
    "Install path prefix, prepended onto install directories.")

if(NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX})
    message(FATAL_ERROR 
        "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}"
        "You must provide an absolute install prefix...
        Try to use `readlink -f` or `realpath` commands")
endif()

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include/mutation++ CACHE PATH
  "Installation directory for header files")
set(INSTALL_CMAKE_DIR share/cmake/mutation++ CACHE PATH
  "Installation directory for CMake files")

foreach (dir LIB BIN INCLUDE CMAKE)
    set(inst_dir INSTALL_${dir}_DIR)
    # Mark install dirs preferences advanced
    mark_as_advanced(${inst_dir})

    # Make the paths absolute (the ${${var}} notation is evaluating the
    # evaluated variable. :: I.e. ${inst_dir} first evaluates, for example, 
    # to INSTALL_LIB_DIR. Then ${${inst_dir}} evaluates to ${INSTALL_LIB_DIR}
    # that respectively evaluates to the directory path containing libraries
    if(NOT IS_ABSOLUTE "${${inst_dir}}")
        set(${inst_dir} "${CMAKE_INSTALL_PREFIX}/${${inst_dir}}")
    endif()
endforeach()


if (CMAKE_COMPILER_IS_GNUCXX)
    #set (CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wold-style-cast")
    set (CMAKE_CXX_FLAGS "-g")
endif ()

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_MACOSX_RPATH 1)
endif()

################################################################################

# Fortran wrapper options
###############################################################################
option (BUILD_FORTRAN_WRAPPER
    "Generate the wrapper library for using mutation++ with Fortran" OFF)

if (BUILD_FORTRAN_WRAPPER)
    # FFLAGS depend on the compiler
    get_filename_component (
        Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

    if (Fortran_COMPILER_NAME STREQUAL "gfortran")
        # gfortran
        set (CMAKE_Fortran_FLAGS         "-fdefault-real-8")
        set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
        set (CMAKE_Fortran_FLAGS_DEBUG   "-g")
    elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
        # ifort (untested)
        set (CMAKE_Fortran_FLAGS         "-r8")
        set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
        set (CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback -fpe0 -check all")
    endif()
endif()



###############################################################################
# Doxygen documentation generation
###############################################################################
option (BUILD_DOCUMENTATION
    "Use Doxygen to create the HTML based API documentation" OFF)

if (BUILD_DOCUMENTATION)
    FIND_PACKAGE(Doxygen)

    if (NOT DOXYGEN_FOUND)
        message(FATAL_ERROR
            "Doxygen is needed to build the documentation. Please install it
             correctly")
    endif()
    # Configure the Template Doxyfile for our specific project
    configure_file(Doxyfile.in
                    ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
    # Add a custom target to run Doxygen when ever the project is built
    add_custom_target (docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
        SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
    add_custom_target(test_docs
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMAND python scripts/test_docs.py)
    add_dependencies(docs test_docs)
endif()

###############################################################################
# ThirdParty Libraries
###############################################################################

# Eigen
find_package(Eigen3)


###############################################################################
# Source code
################################################################################
option(ENABLE_COVERAGE "Generate coverage for codecov.io" OFF)
list(APPEND LCOV_REMOVE_PATTERNS 
    "'/usr/*'"
    "'*/tests/*'"
    "'*/install/*'"
    "'*/examples/*'"
    "'*/thirdparty/*'"
)
include(CodeCoverage)

# Descend into the src directory to build all targets and libraries
add_subdirectory(src)

###############################################################################
# CTest
################################################################################
option(ENABLE_TESTING "Enable regression testing with CTest" OFF)

if (ENABLE_TESTING OR ENABLE_COVERAGE)
    find_package(Catch2)
    enable_testing()
    add_subdirectory(tests)
    add_subdirectory(examples)
endif()


