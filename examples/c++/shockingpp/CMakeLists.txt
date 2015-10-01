# shocking++ code written by Pierre Schrooyen

# Basic setup of project and minimum required CMake version
cmake_minimum_required(VERSION 2.6)
project(shockingpp)
enable_language(CXX)

# Mutation++ directories
set(MPP_INSTALL_DIR "" CACHE PATH "Mutation++ installation directory")
include_directories(${MPP_INSTALL_DIR}/include/mutation++)
link_directories(${MPP_INSTALL_DIR}/lib)

# Boost Library
find_package( Boost COMPONENTS filesystem system REQUIRED )

# Descend into the src directory to build all targets and libraries
add_subdirectory(src)

# Main executable
get_property(shockingpp_SRCS GLOBAL PROPERTY shockingpp_SRCS)
set_source_files_properties(${shockingpp_SRCS} PROPERTIES COMPILE_FLAGS -cpp)
add_executable(shockingpp ${shockingpp_SRCS})

target_link_libraries(shockingpp -lmutation++ ${Boost_LIBRARIES})

# Installation
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS shockingpp DESTINATION bin)

