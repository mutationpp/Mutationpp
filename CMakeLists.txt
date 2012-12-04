cmake_minimum_required(VERSION 2.6)
project(mutation++)

option(BUILD_FORTRAN_WRAPPER "Generate the wrapper library for using mutation++ with Fortran" OFF)
option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)

enable_language(CXX)
if (BUILD_FORTRAN_WRAPPER)
    enable_language(Fortran)
endif()

set (CMAKE_BUILD_TYPE      "Release")
set (CMAKE_INSTALL_PREFIX  "../install")

if (CMAKE_COMPILER_IS_GNUCXX)
    #set (CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wold-style-cast")
    set (CMAKE_CXX_FLAGS "-g")
endif ()

if (BUILD_FORTRAN_WRAPPER)
    # FFLAGS depend on the compiler
    get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

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

#-- Add an Option to toggle the generation of the API documentation
if(BUILD_DOCUMENTATION)
    FIND_PACKAGE(Doxygen)
    if (NOT DOXYGEN_FOUND)
        message(FATAL_ERROR 
            "Doxygen is needed to build the documentation. Please install it correctly")
    endif()
    #-- Configure the Template Doxyfile for our specific project
    configure_file(Doxyfile.in 
                    ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
    #-- Add a custom target to run Doxygen when ever the project is built
    add_custom_target (docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
        SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
endif()

# Descend into the src directory to build all targets and libraries
add_subdirectory(src)

# Install the data directory in the install folder
install(FILES data/thermo/nasa7.dat DESTINATION include/mutation++_data/thermo)
install(FILES data/thermo/nasa9.dat DESTINATION include/mutation++_data/thermo)
install(FILES data/thermo/elements.xml DESTINATION include/mutation++_data/thermo)
install(FILES data/thermo/species.xml DESTINATION include/mutation++_data/thermo)
install(FILES data/transport/heavy.dat DESTINATION include/mutation++_data/transport)
FILE(GLOB mechanisms data/mechanisms/*.xml)
install(FILES ${mechanisms} DESTINATION include/mutation++_data/mechanisms)
FILE(GLOB mixtures data/mixtures/*.xml)
install(FILES ${mixtures} DESTINATION include/mutation++_data/mixtures)

