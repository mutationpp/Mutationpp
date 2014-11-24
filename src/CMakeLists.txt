cmake_minimum_required(VERSION 2.6)

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

add_subdirectory(general)
add_subdirectory(kinetics)
add_subdirectory(numerics)
add_subdirectory(thermo)
add_subdirectory(transfer)
add_subdirectory(transport)
add_subdirectory(utilities)
add_subdirectory(gsi)

#if (ENABLE_TESTING)
#    add_subdirectory(tests)
#endif()

if (BUILD_FORTRAN_WRAPPER)
    add_subdirectory(fortran)
endif()

get_property(mutation++_SRCS GLOBAL PROPERTY mutation++_SRCS)
add_library(mutation++ SHARED ${mutation++_SRCS})
install(TARGETS mutation++ DESTINATION lib)
