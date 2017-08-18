/**
@page installation Installation

@tableofcontents

@section prerequisites Prerequisites

Mutation++ is designed from the ground up to be self sufficient and does not
require any additional packages to use.  All that is required is a C++ compiler, a
Python interpreter (shipped with every UNIX operating system installation) to use
the automated install process for the third-party dependencies and a working installation
of CMake (make sure that command line support is enabled if your computer is
running macOS).

However, if you would like to build this reference manual, then you will need the
[Doxygen](http://www.stack.nl/~dimitri/doxygen/download.html) tool installed.  See
the section on [building the documentation](@ref build_docs) below for more instructions.

@section basic_install Install on Mac OS X or Unix/Linux

First, download the library and unpack the source to a directory where you have write access.

Then, you must add the following environment
variables to your .bashrc (Linux) or .bash_profile (Mac) file located in your
home directory

__Linux__

     export MPP_DIRECTORY="path_to_mutation++_directory"
     export MPP_DATA_DIRECTORY=$MPP_DIRECTORY/data
     export PATH=$MPP_DIRECTORY/install/bin:$PATH
     export LD_LIBRARY_PATH=$MPP_DIRECTORY/install/lib:$LD_LIBRARY_PATH

__Mac__

     export MPP_DIRECTORY="path_to_mutation++_directory"
     export MPP_DATA_DIRECTORY=$MPP_DIRECTORY/data
     export PATH=$MPP_DIRECTORY/install/bin:$PATH
     export DYLD_LIBRARY_PATH=$MPP_DIRECTORY/install/lib:$DYLD_LIBRARY_PATH

Note that you must change the `"path_to_mutation++_directory"` to the correct path
on your machine.

Once the environment variables are set, you must install third-party dependencies.
This can be done manually or by invoking the install script from the root of the
Mutation++ source directory:

    cd $MPP_DIRECTORY
    python scripts/install_dependencies.py

Once the dependencies are installed then use the following commands from the root
of the Mutation++ repository to install.

     mkdir build
     cd build
     cmake ..
     make -j N install

where `N` is the number of CPU units to use for the build process (e.g. 4).

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
