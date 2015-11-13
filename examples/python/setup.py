#!/usr/bin/env python

from setuptools import setup
import shutil

# Fetch library in the CMake build directory
shutil.copy('../../build/src/python/libpyMutation.so', 'pyMutation/')

setup(

    # Library basic information
    name='pyMutation',
    version='0.0.1',

    # List of files to include (use find_packages for recursive search)
    packages=['pyMutation'],
    include_package_data=True,
    package_data={'pyMutation': ['*.so']},

    # Prevent the library package from being compressed
    zip_safe=False,

)
