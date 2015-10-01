cmake_minimum_required(VERSION 2.6)

# Adds sources to a given target
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

# Source code here
add_sources(shockingpp
    shocking.cpp
    InputData.cpp
    InputFileParser.cpp
    MeshOptions.cpp
    PostShockConditions.cpp
    SolveSystem.cpp
)


install(FILES shocking.h DESTINATION bin/include)
install(FILES InputData.h DESTINATION bin/include)
install(FILES InputFileParser.h DESTINATION bin/include)
install(FILES MeshOptions.h DESTINATION bin/include)
install(FILES PostShockConditions.h DESTINATION bin/include)
install(FILES SolveSystem.h DESTINATION bin/include)
