 cmake_minimum_required(VERSION 2.6)
project(Testings)
execute_process(COMMAND echo "THE GCC OF THIS MACHINE IS ${MY_GCC_VERSION}" )

set(BUILD_SHARED_LIBS ON)

add_definitions(-DTEST_DATA_DIRECTORY=${CMAKE_CURRENT_SOURCE_DIR})
#find_package( Boost "1.46.1" COMPONENTS unit_test_framework filesystem system REQUIRED)
find_package( Boost COMPONENTS unit_test_framework filesystem system REQUIRED)
add_executable(test_thermo test_thermo.cpp ThermoTests.h MppTestFixture.h )					
target_link_libraries( test_thermo   mutation++ ${Boost_LIBRARIES} )  
add_test(NAME test_thermo  COMMAND $<TARGET_FILE:test_thermo> )
