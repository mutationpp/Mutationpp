# Installation Guide

## Prerequisites

Compiling and installing Mutation++ and its dependencies requires:

* A C++ compiler [compatible with C++11 revision of the
  language](https://en.cppreference.com/w/cpp/compiler_support)
* A working installation of a fairly recent CMake (>= 3.5)

On macOS, make sure that CMake is installed with command line support
activated. You can use [homebrew](http://brew.sh) or
[macports](http://macports.org) as package managers on macOS. 

## Download

The source code and databases may be downloaded from the Git repository using
the following command on the command line, which should be run from the parent
directory in which you wish to place the Mutation++ library.

```
git clone https://github.com/mutationpp/Mutationpp.git
```

## Install
To begin, add the following environment variables to your .bashrc (Linux) or
.bash_profile (Mac) file located in your home directory

### Dependencies
`mutation++` depends on [`Eigen`](https://eigen.tuxfamily.org) for linear
algebra purposes. We ship them in `thirdparty`but you are free to use another
version or a version installed in the system. Check
[dependencies]('dependencies.md') for more information.

### Linux
```
export MPP_DIRECTORY=path_to_mutation++_directory
export MPP_DATA_DIRECTORY=$MPP_DIRECTORY/data
export PATH=$MPP_DIRECTORY/install/bin:$PATH
export LD_LIBRARY_PATH=$MPP_DIRECTORY/install/lib:$LD_LIBRARY_PATH
```

### Mac
```
export MPP_DIRECTORY=path_to_mutation++_directory
export MPP_DATA_DIRECTORY=$MPP_DIRECTORY/data
export PATH=$MPP_DIRECTORY/install/bin:$PATH
export DYLD_LIBRARY_PATH=$MPP_DIRECTORY/install/lib:$DYLD_LIBRARY_PATH
```

Note that you must change the `path_to_mutation++_directory` to the correct
path on your machine where you cloned the Git repository.


```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=$(realpath ../install) ..
make -j N install
```

where `N` is the number of CPU units to use for the build process (e.g. 4).
This cmake command will install `mutation++` into a directory named `instal` in
the root directory of the git repository (i.e. the directory that contains the
file `README.md`). 
You're free to choose whatever path works better for you.

## Test
A simply way to check that the installation process was successful is to try
the [checkmix](checkmix.md#top) command. For a deeper description of the test
suite check [Testing](testing.md)

```
checkmix air11
```

This should print information about the species and reactions in the air11
mixture.  You can also now use the [mppequil](mppequil.md#top) program for
calculating equilibrium properties of a given mixture.  For more information
about `mppequil`, simply type

```
mppequil --help
```

which will display a help message.  To modify or create new mixtures, go to the
`data/mixtures` folder in your mutation++ directory and open or create the
appropriate file.  Mixture files may also be placed in the local directory in
which you are calling any Mutation++ tool.
