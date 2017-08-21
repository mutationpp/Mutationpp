#
# Copyright 2017 von Karman Institute for Fluid Dynamics (VKI)
#
# This file is part of MUlticomponent Thermodynamic And Transport
# properties for IONized gases in C++ (Mutation++) software package.
#
# Mutation++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Mutation++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with Mutation++.  If not, see
# <http://www.gnu.org/licenses/>.
#
# 
# This file originated from the Catch library contributions as an unlicensed 
# file.  It has been modified in order to better suite the needs of Mutation++.
#
#==================================================================================================#
#  supported macros                                                                                #
#    - TEST_CASE                                                                                   #
#                                                                                                  #
#  Usage                                                                                           #
# 1. make sure this module is in the path or add this otherwise:                                   #
#    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake.modules/")              #
# 2. make sure that you've enabled testing option for the project by the call:                     #
#    enable_testing()                                                                              #
# 3. add the lines to the script for testing target (sample CMakeLists.txt):                       #
#        project(testing_target)                                                                   #
#        set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")          #
#        enable_testing()                                                                          #
#                                                                                                  #
#        find_path(CATCH_INCLUDE_DIR "catch.hpp")                                                  #
#        include_directories(${INCLUDE_DIRECTORIES} ${CATCH_INCLUDE_DIR})                          #
#                                                                                                  #
#        file(GLOB SOURCE_FILES "*.cpp")                                                           #
#        add_executable(${PROJECT_NAME} ${SOURCE_FILES})                                           #
#                                                                                                  #
#        include(ParseAndAddCatchTests)                                                            #
#        ParseAndAddCatchTests(${PROJECT_NAME})                                                    #
#==================================================================================================#

cmake_minimum_required(VERSION 2.8.8)

# This removes the contents between
#  - block comments (i.e. /* ... */) 
#  - full line comments (i.e. // ... ) 
# contents have been read into '${CppCode}'.
# !keep partial line comments
function(RemoveComments CppCode)
  string(ASCII 2 CMakeBeginBlockComment)
  string(ASCII 3 CMakeEndBlockComment)
  string(REGEX REPLACE "/\\*" "${CMakeBeginBlockComment}" ${CppCode} "${${CppCode}}")
  string(REGEX REPLACE "\\*/" "${CMakeEndBlockComment}" ${CppCode} "${${CppCode}}")
  string(REGEX REPLACE "${CMakeBeginBlockComment}[^${CMakeEndBlockComment}]*${CMakeEndBlockComment}" "" ${CppCode} "${${CppCode}}")
  string(REGEX REPLACE "\n[ \t]*//+[^\n]+" "\n" ${CppCode} "${${CppCode}}")

  set(${CppCode} "${${CppCode}}" PARENT_SCOPE)
endfunction()

# Worker function
function(ParseFile SourceFile TestTarget)
    if(NOT EXISTS ${SourceFile})
        return()
    endif()
    file(STRINGS ${SourceFile} Contents NEWLINE_CONSUME)

    # Remove block and fullline comments
    RemoveComments(Contents)

    # Find definition of test names
    string(REGEX MATCHALL "[ \t]*TEST_CASE[^(]*\\([^\"]*\"([^\"]+)\"[^,]*,[^\"]*\"([^\"]+)\"[^)]*\\)" Tests "${Contents}")

    foreach(TestName ${Tests})
        # Strip newlines
        string(REGEX REPLACE "\\\\\n|\n" "" TestName "${TestName}")

        # Get string parts of test definition
        string(REGEX MATCHALL "\"+([^\\^\"]|\\\\\")+\"+" TestStrings "${TestName}")

        # Strip wrapping quotation marks
        string(REGEX REPLACE "^\"(.*)\"$" "\\1" TestStrings "${TestStrings}")
        string(REPLACE "\";\"" ";" TestStrings "${TestStrings}")

        # Validate that a test name and tags have been provided
        list(LENGTH TestStrings TestStringsLength)
        if(TestStringsLength GREATER 2 OR TestStringsLength LESS 1)
            message(FATAL_ERROR "You must provide a valid test name and tags for all tests in ${SourceFile}")
        endif()

        # Assign name and tags
        list(GET TestStrings 0 Name)
        set(CTestName "${Name}")

        if(TestStringsLength EQUAL 2)
            list(GET TestStrings 1 Tags)
            string(TOLOWER "${Tags}" Tags)
            string(REPLACE "]" ";" Tags "${Tags}")
            string(REPLACE "[" "" Tags "${Tags}")
        endif()

        # Add the test and set its properties
        add_test(NAME "\"${CTestName}\"" COMMAND ${TestTarget} ${Name} ${AdditionalCatchParameters})
        set_tests_properties("\"${CTestName}\"" PROPERTIES FAIL_REGULAR_EXPRESSION "No tests ran"
                                                LABELS "${Tags}")

    endforeach()
endfunction()

# entry point
function(ParseAndAddCatchTests TestTarget)
    get_target_property(SourceFiles ${TestTarget} SOURCES)
    foreach(SourceFile ${SourceFiles})
        ParseFile(${SourceFile} ${TestTarget})
    endforeach()
endfunction()
