cmake_minimum_required(VERSION 2.6)

add_library(mutation++_fortran SHARED cwrapper.cpp cwrapper_interface.f90)
target_link_libraries(mutation++_fortran mutation++)
install(TARGETS mutation++_fortran DESTINATION lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mutationpp.mod DESTINATION include/mutation++)

