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

cmake_minimum_required(VERSION 2.6)
project(mutation++)

set(MUTATIONPP_MAJOR_VERSION 0)
set(MUTATIONPP_MINOR_VERSION 2)
set(MUTATIONPP_VERSION 
    ${MUTATIONPP_MAJOR_VERSION}.${MUTATIONPP_MINOR_VERSION})

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
find_package(Eigen REQUIRED)
include_directories(${Eigen_INCLUDE_DIRS})


###############################################################################
# Source code
<<<<<<< HEAD
################################################################################
option(ENABLE_COVERAGE "Generate coverage for codecov.io" OFF)
find_package(codecov)

=======
###############################################################################
<<<<<<< HEAD
>>>>>>> WIP
=======

>>>>>>> Add exporting of Configuration File
# Descend into the src directory to build all targets and libraries
include_directories(
    ${CMAKE_SOURCE_DIR}/install/include
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/general
    ${CMAKE_CURRENT_SOURCE_DIR}/src/kinetics
    ${CMAKE_CURRENT_SOURCE_DIR}/src/numerics
    ${CMAKE_CURRENT_SOURCE_DIR}/src/thermo
    ${CMAKE_CURRENT_SOURCE_DIR}/src/transfer
    ${CMAKE_CURRENT_SOURCE_DIR}/src/transport
    ${CMAKE_CURRENT_SOURCE_DIR}/src/utilities
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gsi
)
add_subdirectory(src)

###############################################################################
# Exporting configuration files to be used by other projects that need 
# mutation++
###############################################################################

# Create the configure file adding the correct paths
# tutorials/How-to-create-a-ProjectConfig.cmake-file 
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")

# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(mutation++Config.cmake.in
  "${PROJECT_BINARY_DIR}/mutation++Config.cmake" @ONLY)

# ... for the install tree
set(CONF_INCLUDE_DIRS "\${MUTATIONPP_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(mutation++Config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mutation++Config.cmake" @ONLY)

# ... for both
configure_file(mutation++ConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/mutation++ConfigVersion.cmake" @ONLY)

# Install the mutation++Config.cmake and mutation++ConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mutation++Config.cmake"
  "${PROJECT_BINARY_DIR}/mutation++ConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)


# Export targets to be used by projects that will use mutation++
export(TARGETS mutation++ ${mutation_exes} FILE 
    "${PROJECT_BINARY_DIR}/mutation++Targets.cmake")

# Install the export set for use with the install-tree
install(EXPORT mutation++Targets DESTINATION
  "${INSTALL_CMAKE_DIR}" COMPONENT dev)

###############################################################################
# CTest
################################################################################
option(ENABLE_TESTING "Enable regression testing with CTest" OFF)

if (ENABLE_TESTING OR ENABLE_COVERAGE)
    # Look for Catch2
    find_package(Catch2 REQUIRED)
    enable_testing()
    add_subdirectory(tests)
    add_subdirectory(examples)

endif(ENABLE_TESTING OR ${CMAKE_BUILD_TYPE} MATCHES "Coverage")

list(APPEND LCOV_REMOVE_PATTERNS 
    "/usr/*"
    "'${PROJECT_SOURCE_DIR}/tests/*'"
    "'${PROJECT_SOURCE_DIR}/install/*'"
    "'${PROJECT_SOURCE_DIR}/thirdparty/*'"
)
coverage_evaluate()
