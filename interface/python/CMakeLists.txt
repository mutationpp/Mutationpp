cmake_minimum_required(VERSION 2.8.12)
project(mutationpp)

set(MPP_INSTALL_DIR  ../../install)
include_directories(${MPP_INSTALL_DIR}/include/mutation++)
include_directories(${MPP_INSTALL_DIR}/include/)
link_directories(${MPP_INSTALL_DIR}/lib)

include_directories(../../thirdparty/eigen/)

find_package(pybind11 REQUIRED)

pybind11_add_module(mutationpp
    src/mutationpp_python.cpp
    src/init_MixtureOptions.cpp
)
target_link_libraries(mutationpp PRIVATE mutation++.dylib)