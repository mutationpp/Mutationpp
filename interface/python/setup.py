import re
import os
import sys
import subprocess
import platform

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

setup(
    name='mutationpp',
    author='von Karman Institute for Fluid Dynamics',
    version='0.0.1',    
    package_dir={'': '${CMAKE_CURRENT_SOURCE_DIR}'},
    author_email='mutationpp@vki.ac.be',
    description='Python bindings for Mutation++',
    ext_modules=[CMakeExtension('mutationpp')],
    zip_safe=False,
)
