cmake_minimum_required(VERSION 2.8.12)
project(mutationpp_python)

set(MPP_INSTALL_DIR  ../../install)
include_directories(${MPP_INSTALL_DIR}/include/mutation++)
include_directories(${MPP_INSTALL_DIR}/include/)
link_directories(${MPP_INSTALL_DIR}/lib)

include_directories(../../thirdparty/eigen/)
add_subdirectory(pybind11)

pybind11_add_module(mutationpp_python
    src/mutationpp_python.cpp
    src/init_MixtureOptions.cpp
)
target_link_libraries(mutationpp_python PRIVATE mutation++.dylib)