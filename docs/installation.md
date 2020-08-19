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
### Dependencies
`mutation++` depends on [`Eigen`](https://eigen.tuxfamily.org) for linear
algebra purposes. We ship them in `thirdparty`but you are free to use another
version or a version installed in the system. Check
[dependencies](dependencies.md) for more information.

### Build
In order to build the library, create a `build` directory, than run the `cmake`
command inside of it (Note: The commands listed below installs `Mutationpp` 
inside the parent directory containing the code using the `CMAKE_INSTALL_PREFIX`
variable. If you don't specify it, the code will be installed in the default
prefix `/usr/local`. 

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=$(realpath ../install) ..
make -j N install
```

where `N` is the number of CPU units to use for the build process (e.g. 4).
This cmake command will install `mutation++` into the directory specified by `CMAKE_INSTALL_PREFIX`.
You're free to choose whatever path works better for you.

### Configuration
**Linux**<br>
Some Linux distributions install libraries in `${CMAKE_INSTALL_PREFIX}/lib64` 
instead of `${CMAKE_INSTALL_PREFIX}/lib`. For this reason, after you compiled
and installed the code, check if you have a `lib` or `lib64` `LIBDIR`. Next, add the following lines
to your `.bashrc` file in your home directory.

```
export MPP_DIRECTORY=path_to_mutation++_directory
export MPP_DATA_DIRECTORY=$MPP_DIRECTORY/data
export PATH=$MPP_DIRECTORY/install/bin:$PATH
export LD_LIBRARY_PATH=$MPP_DIRECTORY/install/<LIBDIR>:$LD_LIBRARY_PATH  # Change <LIBDIR> with `lib` or `lib64`
```

**MacOS**<br>
Add the following lines to your `.bash_profile` file in your home directory.
```
export MPP_DIRECTORY=path_to_mutation++_directory
export MPP_DATA_DIRECTORY=$MPP_DIRECTORY/data
export PATH=$MPP_DIRECTORY/install/bin:$PATH
export DYLD_LIBRARY_PATH=$MPP_DIRECTORY/install/lib:$DYLD_LIBRARY_PATH
```

*Note that you must change the `path_to_mutation++_directory` to the correct
path on your machine where you cloned the Git repository* **for both Linux and MacOS installs**.


## Test
A simple way to check that the installation process was successful is to try
the [checkmix](checkmix.md#top) command. For a deeper description of the test
suite check [Testing](testing.md)

```
checkmix air_11
```

This should print information about the species and reactions in the air_11
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


## Build API Docs
If you are interested in building the API documentation for Mutation++ locally, then be sure to have [Doxygen](https://www.doxygen.nl/index.html) installed on your computer and run

```
make docs
```

from the `build` directory.  Once done, open the `docs/html/index.html` file using your favorite web browser to see a full listing of all the classes and methods available in the library.
