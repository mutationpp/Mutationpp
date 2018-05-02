# Distributed under the BSD 2-clause License

# FindCatch2
# --------
#
# Find Catch2 header-only library
#
# This module tries before to search for cmake-exported configuration files
# using the standard 
#
# This module sets the following variables:
#
# ::
#
#   Catch2_FOUND - set to true if a library implementing the BLAS interface
#     is found
#
#   Catch2_INCLUDE_DIRS - the directory containing the Catch2 header files
#
# 
# To give an hint for the directory to look for the user must set Catch23_DIR
#
#
# Author: Ruben Di Battista <rubendibattista@gmail.com>

message(STATUS "Searching for Catch2...")

set(Catch2_CHECK_PATHS 
    "/usr/local/include"
    "/usr/local/homebrew/include" # Mac OS X
    "/opt/local/var/macports/software" # Mac OS X.
    "/opt/local/include"
    "/usr/include"
    "$ENV{CPLUS_INCLUDE_PATH}"
    "$ENV{CPATH}"
    "${CMAKE_SOURCE_DIR}/thirdparty/catch"
    )

# First try to use the standard find_package that should be able to find the
# .cmake configuration files if installed in standard directories. 

find_package(Catch2 QUIET
    CONFIG)

if (NOT Catch2_FOUND)
    # Standard stuff did not work, so we try to locate Catch2 in include paths
    find_path(Catch2_INCLUDE_DIR NAMES catch.hpp
        HINTS ${Catch2_DIR}
        PATHS ${Catch2_CHECK_PATHS}
        PATH_SUFFIXES catch include)

    if (Catch2_INCLUDE_DIR)
        # We found the include dir, we add it to the standard variable
        # ${<library>_INCLUDE_DIRS}
        set (Catch2_INCLUDE_DIRS ${Catch2_INCLUDE_DIR})
    endif (Catch2_INCLUDE_DIR)
endif (NOT Catch2_FOUND)


# Handle REQUIRED/QUIET optional arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Catch2
  REQUIRED_VARS Catch2_INCLUDE_DIRS)

