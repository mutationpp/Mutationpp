cmake_minimum_required(VERSION 2.6)

add_library(mutation++_fortran SHARED cwrapper.cpp cwrapper_interface.f90)
target_link_libraries(mutation++_fortran mutation++)
install(TARGETS mutation++_fortran DESTINATION lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mutationpp.mod DESTINATION include/mutation++)

if (BUILD_TEST_CODES)
    add_executable(wrapper_test wrapper_test.f90)
    target_link_libraries(wrapper_test mutation++_fortran)
    set_target_properties(wrapper_test PROPERTIES LINKER_LANGUAGE Fortran)
    install(TARGETS wrapper_test DESTINATION bin)
endif()

