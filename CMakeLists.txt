################################################################################
#
# High level cmake configuration file for the Mutation++ project and associated
# helper tools and binaries.  All of the various configuration options are
# defined here.
#
# author: J.B. Scoggins (scoggins@vki.ac.be)
#
################################################################################
cmake_minimum_required(VERSION 2.6)
project(mutation++)


enable_language(CXX)
if (BUILD_FORTRAN_WRAPPER)
    enable_language(Fortran)
endif()


################################################################################
# Build types
################################################################################

# Profile
SET (CMAKE_CXX_FLAGS_PROFILE "-g -Wall -pg -O3 -DNDEBUG" CACHE STRING
    "Flags used by the C++ compiler during Profile builds."
    FORCE )
SET (CMAKE_C_FLAGS_PROFILE "-g -Wall -pedantic -pg -O3 -DNDEBUG" CACHE STRING
    "Flags used by the C compiler during Profile builds."
    FORCE )
SET (CMAKE_EXE_LINKER_FLAGS_PROFILE
    "" CACHE STRING
    "Flags used for linking binaries during Profile builds."
    FORCE )
SET (CMAKE_SHARED_LINKER_FLAGS_PROFILE
    "" CACHE STRING
    "Flags used by the shared libraries linker during Profile builds."
    FORCE )
MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_PROFILE
    CMAKE_C_FLAGS_PROFILE
    CMAKE_EXE_LINKER_FLAGS_PROFILE
    CMAKE_SHARED_LINKER_FLAGS_PROFILE)
    

# Update the documentation string of CMAKE_BUILD_TYPE for GUIs and set the
# default build type as Release
if (NOT CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: None Debug Release 
        RelWithDebInfo MinSizeRel Profile." FORCE)
endif ()


################################################################################
# Install prefix settings
################################################################################
SET (CMAKE_INSTALL_PREFIX "../install" CACHE STRING
    "Install path prefix, prepended onto install directories." FORCE)
MARK_AS_ADVANCED (CMAKE_INSTALL_PREFIX)

if (CMAKE_COMPILER_IS_GNUCXX)
    #set (CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wold-style-cast")
    set (CMAKE_CXX_FLAGS "-g")
endif ()


################################################################################
# Fortran wrapper options
################################################################################
option (BUILD_FORTRAN_WRAPPER 
    "Generate the wrapper library for using mutation++ with Fortran" OFF)

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



################################################################################
# Doxygen documentation generation
################################################################################
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
endif()


################################################################################
# Source code
################################################################################
# Descend into the src directory to build all targets and libraries
add_subdirectory(src)
add_subdirectory(mutation_tests)


################################################################################
# Install the data directory in the install folder
################################################################################
install(FILES data/thermo/nasa7.dat DESTINATION include/mutation++_data/thermo)
install(FILES data/thermo/nasa9.dat DESTINATION include/mutation++_data/thermo)
install(FILES data/thermo/elements.xml DESTINATION include/mutation++_data/thermo)
install(FILES data/thermo/species.xml DESTINATION include/mutation++_data/thermo)
install(FILES data/transport/heavy.dat DESTINATION include/mutation++_data/transport)
FILE(GLOB mechanisms data/mechanisms/*.xml)
install(FILES ${mechanisms} DESTINATION include/mutation++_data/mechanisms)
FILE(GLOB mixtures data/mixtures/*.xml)
install(FILES ${mixtures} DESTINATION include/mutation++_data/mixtures)

