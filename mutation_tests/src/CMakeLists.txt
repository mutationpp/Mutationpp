cmake_minimum_required(VERSION 2.6)
 
project(myTests)
ENABLE_TESTING()

set(BUILD_SHARED_LIBS ON)

find_package( Boost "1.46.1" COMPONENTS unit_test_framework filesystem system REQUIRED)

#add_subdirectory(src)
add_subdirectory(test_src)
