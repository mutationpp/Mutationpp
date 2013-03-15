# cmake_minimum_required(VERSION 2.6)
 
# project(Testings)

# find_package( Boost "1.46.1" COMPONENTS unit_test_framework REQUIRED)
 
# add_executable(test_class_Position test_class_Position.cpp)
add_executable(test_class_readFile test_class_readFile.cpp readFile.cpp readFile.hpp)

target_link_libraries( test_class_readFile ${Boost_LIBRARIES} ) 
 
add_test(NAME test_main_name_displayed COMMAND $<TARGET_FILE:test_class_readFile> )
