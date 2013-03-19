# cmake_minimum_required(VERSION 2.6)
 
# project(Testings)

# find_package( Boost "1.46.1" COMPONENTS unit_test_framework REQUIRED)
 
# add_executable(test_class_Position test_class_Position.cpp)
add_executable(test test_thermo.cpp MppTestFixture.hpp)

target_link_libraries( test_class_readFile /home/didi/mutation++/trunk/install/include/mutation++ 
			/home/didi/mutation++/trunk/install/lib/libmutation++.so ${Boost_LIBRARIES} ) 
 
add_test(NAME Thermodynamics COMMAND $<TARGET_FILE:test_class_readFile> )
