# Distributed under the BSD 2-clause License

# FindEigen
# --------
#
# Find Eigen header-only library
#
# This module tries before to search for cmake-exported configuration files
# using the standard 
#
# This module sets the following variables:
#
# ::
#
#   Eigen_FOUND - set to true if a library implementing the BLAS interface
#     is found
#
#   Eigen_INCLUDE_DIRS - the directory containing the Eigen header files
#
# 
# To give an hint for the directory to look for the user must set Eigen3_DIR
#
#
# Author: Ruben Di Battista <rubendibattista@gmail.com>

message(STATUS "Searching for Eigen...")

set(Eigen_CHECK_PATHS 
    "/usr/local/include"
    "/usr/local/homebrew/include" # Mac OS X
    "/opt/local/var/macports/software" # Mac OS X.
    "/opt/local/include"
    "/usr/include"
    "$ENV{CPLUS_INCLUDE_PATH}"
    "$ENV{CPATH}"
    "${CMAKE_SOURCE_DIR}/thirdparty/eigen"
    )

# First try to use the standard find_package that should be able to find the
# .cmake configuration files if installed in standard directories. 

find_package(Eigen QUIET
    CONFIG)

if (NOT Eigen_FOUND)
    # Standard stuff did not work, so we try to locate Eigen in include paths
    find_path(Eigen_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
        HINTS ${Eigen_DIR}
        PATHS ${Eigen_CHECK_PATHS}
        PATH_SUFFIXES eigen3 Eigen3 eigen Eigen)

    if (Eigen_INCLUDE_DIR)
        # We found the include dir, we add it to the standard variable
        # ${<library>_INCLUDE_DIRS}
        set (Eigen_INCLUDE_DIRS ${Eigen_INCLUDE_DIR})
    endif (Eigen_INCLUDE_DIR)
endif (NOT Eigen_FOUND)


# Handle REQUIRED/QUIET optional arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen
  REQUIRED_VARS Eigen_INCLUDE_DIRS)

