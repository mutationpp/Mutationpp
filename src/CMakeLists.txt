#
# Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
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

#- Add sources for a target
#
#  ADD_SOURCES(<target> <source1> [<source2> ...])
#
function (add_sources target)
  # define the <target>_SRCS properties if necessary
  get_property(prop_defined GLOBAL PROPERTY ${target}_SRCS DEFINED)
  if(NOT prop_defined)
    define_property(GLOBAL PROPERTY ${target}_SRCS
      BRIEF_DOCS "Sources for the ${target} target"
      FULL_DOCS "List of source files for the ${target} target")
  endif()
  # create list of sources (absolute paths)
  set(SRCS)
  foreach(src IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${src}")
      get_filename_component(src "${src}" ABSOLUTE)
    endif()
    list(APPEND SRCS "${src}")
  endforeach()
  # append to global property
  set_property(GLOBAL APPEND PROPERTY "${target}_SRCS" "${SRCS}")
endfunction()

# - Add headers for a target (in order to be installed later)
function (add_headers target)
  # define the <target>_HDRS properties if necessary
  get_property(prop_defined GLOBAL PROPERTY ${target}_HDRS DEFINED)
  if(NOT prop_defined)
      define_property(GLOBAL PROPERTY ${target}_HDRS
          BRIEF_DOCS "Headers for the ${target} target"
      FULL_DOCS "List of header files for the ${target} target")
  endif()
  # create list of sources (absolute paths)
  set(HDRS)
  foreach(hfile IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${hfile}")
      get_filename_component(hfile "${hfile}" ABSOLUTE)
    endif()
    list(APPEND HDRS "${hfile}")
  endforeach()
  # append to global property
  set_property(GLOBAL APPEND PROPERTY "${target}_HDRS" "${HDRS}")
endfunction()

add_subdirectory(general)
add_subdirectory(kinetics)
add_subdirectory(numerics)
add_subdirectory(thermo)
add_subdirectory(transfer)
add_subdirectory(transport)
add_subdirectory(utilities)
add_subdirectory(gsi)

if (BUILD_FORTRAN_WRAPPER)
    add_subdirectory(fortran)
endif()



get_property(mutation++_SRCS GLOBAL PROPERTY mutation++_SRCS)
add_library(mutation++ SHARED ${mutation++_SRCS})
add_library(mutation++::mutation++ ALIAS mutation++)
target_include_directories(mutation++
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/general>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/kinetics>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/numerics>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/thermo>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/transfer>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/transport>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/utilities>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gsi>
        $<INSTALL_INTERFACE:include/mutation++>
)
target_link_libraries(mutation++ 
    PUBLIC
        Eigen3::Eigen
)

# Evaluate coverage
if(ENABLE_COVERAGE)
    add_coverage(mutation++)
endif()

# Install and export targets
include(GNUInstallDirs)

# Version configuration file
include(CMakePackageConfigHelpers)

# This actually does nothing since we do not have anything fancy or 
# replaceable in `mutation++Config.cmake.in`. But here to comply with 
# standard practices (and also to copy the config file in the 
# CMAKE_CURRENT_BINARY_DIR)
configure_package_config_file(
    "mutation++Config.cmake.in"
    "mutation++Config.cmake"
    INSTALL_DESTINATION 
        ${CMAKE_INSTALL_LIBDIR}/cmake/mutation++
)


write_basic_package_version_file("mutation++ConfigVersion.cmake"
    VERSION ${mutation++_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES 
    "${CMAKE_CURRENT_BINARY_DIR}/mutation++Config.cmake" 
    "${CMAKE_CURRENT_BINARY_DIR}/mutation++ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mutation++
)

# Install Headers
get_property(mutation++_HDRS GLOBAL PROPERTY mutation++_HDRS)
set_target_properties(mutation++ PROPERTIES PUBLIC_HEADER "${mutation++_HDRS}")
set(MUTATIONPP_INCLUDE_DIR "${CMAKE_INSTALL_INCLUDEDIR}/mutation++")

install(TARGETS mutation++ ${mutation_exes} EXPORT mutation++Targets
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${MUTATIONPP_INCLUDE_DIR}
    PUBLIC_HEADER DESTINATION ${MUTATIONPP_INCLUDE_DIR}
)

export(EXPORT mutation++Targets
    FILE mutation++Targets.cmake
    NAMESPACE mutation++::
)

install(EXPORT mutation++Targets
    FILE        mutation++Targets.cmake
    NAMESPACE   mutation++::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/mutation++"
    EXPORT_LINK_INTERFACE_LIBRARIES
)
