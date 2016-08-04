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
    main.cpp
    Solver.cpp
    DataShocking1T.cpp
    DataShockingNT.cpp
    DataLarsen1T.cpp
    DataLarsenTTv.cpp
    Mesh.cpp
    Shocking1T.cpp
    ShockingTTv.cpp
    Larsen1T.cpp
    LarsenTTv.cpp
    RankineHugoniot1T.cpp
    RankineHugoniotNT.cpp
    SetupProblem.cpp
    SetupShocking1T.cpp
    SetupShockingNT.cpp
    SetupLarsen1T.cpp
    SetupLarsenTTv.cpp
    SetupProperties.cpp
)

install(FILES Solver.h DESTINATION bin/include)
install(FILES Data.h DESTINATION bin/include)
install(FILES DataShocking1T.h DESTINATION bin/include)
install(FILES DataLarsen1T.h DESTINATION bin/include)
install(FILES DataLarsenTTv.h DESTINATION bin/include)
install(FILES Mesh.h DESTINATION bin/include)
install(FILES Problem.h DESTINATION bin/include)
install(FILES Shocking1T.h DESTINATION bin/include)
install(FILES ShockingTTv.h DESTINATION bin/include)
install(FILES Larsen1T.h DESTINATION bin/include)
install(FILES LarsenTTv.h DESTINATION bin/include)
install(FILES ShockRelations.h DESTINATION bin/include)
install(FILES RankineHugoniot1T.h DESTINATION bin/include)
install(FILES SetupProblem.h DESTINATION bin/include)
install(FILES SetupShocking1T.h DESTINATION bin/include)
install(FILES SetupShockingNT.h DESTINATION bin/include)
install(FILES SetupLarsen1T.h DESTINATION bin/include)
install(FILES SetupLarsenTTv.h DESTINATION bin/include)
install(FILES SetupProperties.h DESTINATION bin/include)
install(FILES BasicOutput.h DESTINATION bin/include)
install(FILES dfdt.h DESTINATION bin/include)
install(FILES jacobian.h DESTINATION bin/include)
