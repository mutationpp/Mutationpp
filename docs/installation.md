# Installation Guide

## Prerequisites

Compiling and installing Mutation++ and its dependencies requires:

* A C++ compiler
* A Python interpreter (2.7)
* A working installation of CMake

On macOS, make sure that CMake is installed with command line support activated.

## Download

The source code and databases may be downloaded from the Git repository using the following command on the command line, which should be run from the parent directory in which you wish to place the Mutation++ library.

```
git clone https://github.com/mutationpp/Mutationpp.git
```

## Install
To begin, add the following environment variables to your .bashrc (Linux) or .bash_profile (Mac) file located in your home directory

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

Note that you must change the `path_to_mutation++_directory` to the correct path
on your machine.

Once the environment variables are set, you must install third-party dependencies.
This can be done manually or by invoking the install script from the root of the
Mutation++ source directory:

```
cd $MPP_DIRECTORY
python scripts/install_dependencies.py
```

Once the dependencies are installed then use the following commands from the root
of the Mutation++ repository to install.

```
mkdir build
cd build
cmake ..
make -j N install
```

where `N` is the number of CPU units to use for the build process (e.g. 4).

## Test
A simply way to check that the installation process was successful is to try the [checkmix](checkmix.md#top) command. 

```
checkmix air11
```

This should print information about the species and reactions in the air11 mixture.  You can also now use the [mppequil](mppequil.md#top) program for calculating equilibrium properties of a given mixture.  For more information about `mppequil`, simply type

```
mppequil --help
```

which will display a help message.  To modify or create new mixtures, go to the ```data/mixtures``` folder in your mutation++ directory and open or create the appropriate file.  Mixture files may also be placed in the local directory in which you are calling any Mutation++ tool.
