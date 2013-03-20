# cmake_minimum_required(VERSION 2.6)
 
# project(Testings)

# find_package( Boost "1.46.1" COMPONENTS unit_test_framework REQUIRED)
 
# add_executable(test_class_Position test_class_Position.cpp)
add_executable(test_thermo_class test_thermo.cpp  MppTestFixture.h)
add_executable(mixture test_thermo.cpp  MppTestFixture.h)

target_link_libraries( test_thermo_class  
			/home/didi/mutation++/trunk/install/lib/libmutation++.so ${Boost_LIBRARIES} ) 
target_link_libraries( mixture  
			/home/didi/mutation++/trunk/install/lib/libmutation++.so ${Boost_LIBRARIES} ) 
 
add_test(NAME test_thermo_class  COMMAND $<TARGET_FILE:test_thermo_class> )

set_property(TARGET Mixture PROPERTY LABELS Mixture)
add_test(NAME mixture  COMMAND $<TARGET_FILE:mixture> )
set_property(TEST mixture PROPERTY LABELS Mixture)

