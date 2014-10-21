cmake_minimum_required(VERSION 2.6)

find_package( Boost COMPONENTS unit_test_framework filesystem system REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)

add_executable(stefan_maxwell test_stefan_maxwell.cpp)
target_link_libraries(stefan_maxwell mutation++ ${Boost_LIBRARIES})
add_executable(wdot test_wdot.cpp)
target_link_libraries(wdot mutation++ ${Boost_LIBRARIES})

# Define the tests to be run
add_test(stefan-maxwell_air5 stefan_maxwell air5)
add_test(stefan-maxwell_air11 stefan_maxwell air11)
add_test(wdot_air5  wdot air5)
add_test(wdot_air11 wdot air11)

