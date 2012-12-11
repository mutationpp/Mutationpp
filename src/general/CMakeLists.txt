cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    MixtureOptions.cpp
)

# Build the checkmix code
add_executable(checkmix checkmix.cpp)
target_link_libraries(checkmix mutation++)
install(TARGETS checkmix DESTINATION bin)

# Build the mppequil code
add_executable(mppequil mppequil.cpp)
target_link_libraries(mppequil mutation++)
install(TARGETS mppequil DESTINATION bin)

#add_executable(bprime bprime.cpp)
#target_link_libraries(bprime mutation++)
#install(TARGETS bprime DESTINATION bin)

# Build the fortran wrapper library
if (BUILD_FORTRAN_WRAPPER)
    add_library(mutation++_fortran STATIC cwrapper.cpp cwrapper_interface.f90)
    target_link_libraries(mutation++_fortran mutation++)
    install(TARGETS mutation++_fortran DESTINATION lib)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mutationpp.mod DESTINATION include/mutation++)
endif()

# Build and install the test codes
option(BUILD_TEST_CODES "Build test executables" ON)
if (BUILD_TEST_CODES)
    # test
    add_executable(test test.cpp)
    target_link_libraries(test mutation++)
    install(TARGETS test DESTINATION bin)
    
    # profiling
    #add_executable(profiling profiling.cpp)
    #target_link_libraries(profiling mutation++ rt)
    #install(TARGETS profiling DESTINATION bin)
    
    if (BUILD_FORTRAN_WRAPPER)
        add_executable(wrapper_test wrapper_test.f90)
        target_link_libraries(wrapper_test mutation++_fortran)
        set_target_properties(wrapper_test PROPERTIES LINKER_LANGUAGE Fortran)
        install(TARGETS wrapper_test DESTINATION bin)
    endif()
endif()

# Install the header files
install(FILES mutation++.h DESTINATION include/mutation++)
install(FILES Mixture.h DESTINATION include/mutation++)
install(FILES MixtureOptions.h DESTINATION include/mutation++)
install(FILES Constants.h DESTINATION include/mutation++)



