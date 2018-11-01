#- Add sources for a target
#
#  add_sources(<target> <source1> [<source2> ...])
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
#
#  add_headers(<target> <source1> [<source2> ...])
#
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
