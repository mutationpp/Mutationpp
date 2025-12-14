#                                               -*- cmake -*-
#
#  UseEigen3.cmake
#
#  ----------------------------------------------------------------------
#  This file provides legacy support for including Eigen in CMake.
#  The modern and preferred way is to use:
#
#      find_package (Eigen3 3.4 REQUIRED NO_MODULE)
#      add_executable (example example.cpp)
#      target_link_libraries (example Eigen3::Eigen)
#
#  Use the commands below only if you must maintain compatibility
#  with older CMake/Eigen projects.
#  For more information: 
#  https://libeigen.gitlab.io/eigen/docs-3.4/TopicCMakeGuide.html
#  ----------------------------------------------------------------------

add_definitions     ( ${EIGEN3_DEFINITIONS} )
include_directories ( ${EIGEN3_INCLUDE_DIRS} )
