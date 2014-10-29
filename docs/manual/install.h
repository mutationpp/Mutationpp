/**
@page installation Installation

@section prerequisits Prerequisits

Mutation++ is designed from the ground up to be self sufficient and does not
require any additional packages to use.  All that is required is a C++ compiler.

However, if you would like to build this reference manual, then you will need the
[Doxygen](http://www.stack.nl/~dimitri/doxygen/download.html) tool installed.  See
the section on [building the documentation](@ref build_docs) below for more instructions.

If you would like to run the unit tests which are packaged with Mutation++, then
you will need to have the [Boost](http://www.boost.org) libraries installed.  See
the section on [running the unit tests](@ref unit_tests) below.

@section basic_install Install on Mac OS X or Unix/Linux

Once the library is downloaded then use the following commands to install.

     cd mutation++; mkdir build; cd build;
     cmake ..; make install

Once the code is compiled and installed, you must add the following environment 
variables to your .bashrc (Linux) or .bash_profile (Mac) file located in your 
home directory

__Linux__

     export MPP_DATA_DIRECTORY="path_to_mutation++_directory"/data
     export PATH="path_to_mutation++_directory"/install/bin:$PATH
     export LD_LIBRARY_PATH="path_to_mutation++_directory"/install/lib:$LD_LIBRARY_PATH

__Mac__

     export MPP_DATA_DIRECTORY="path_to_mutation++_directory"/data
     export PATH="path_to_mutation++_directory"/install/bin:$PATH
     export DYLD_LIBRARY_PATH="path_to_mutation++_directory"/install/lib:$DYLD_LIBRARY_PATH

Note that you must change the "path_to_mutation++_directory" to the correct path 
on your machine.

@section test_installation Test Installation

To see if the installation is successful, use the [checkmix](@ref checkmix) command
and make sure that mixture information is output for the given mixture.

     checkmix air11

This should print out information about the species and reactions in the air11 mixture. 
You can also now use the [mppequil](@ref mppequil) program for calculating equilibrium 
properties of a given mixture.

@subsection unit_tests Running the Unit Tests

@section build_docs Compiling the Reference Docs

*/