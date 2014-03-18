cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    Mixture.cpp
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

add_executable(bprime bprime.cpp)
target_link_libraries(bprime mutation++)
install(TARGETS bprime DESTINATION bin)

add_executable(reduce reduce.cpp)
target_link_libraries(reduce mutation++)
install(TARGETS reduce DESTINATION bin)

# Build and install the test codes
#option(BUILD_TEST_CODES "Build test executables" OFF)
#if (BUILD_TEST_CODES)
#    
#    # profiling
#    #add_executable(profiling profiling.cpp)
#    #target_link_libraries(profiling mutation++ rt)
#    #install(TARGETS profiling DESTINATION bin)
#    
#    add_executable(minres test_minres.cpp)
#    #target_link_libraries(minres mutation++)
#    install(TARGETS minres DESTINATION bin)
#
#    add_executable(mpgfc test_mpgfc.cpp)
#    target_link_libraries(mpgfc mutation++)
#    install(TARGETS mpgfc DESTINATION bin)
#endif()

# Install the header files
install(FILES mutation++.h DESTINATION include/mutation++)
install(FILES Mixture.h DESTINATION include/mutation++)
install(FILES MixtureOptions.h DESTINATION include/mutation++)
install(FILES Constants.h DESTINATION include/mutation++)



