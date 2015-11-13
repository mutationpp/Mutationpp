# pyMutation: Python bindings for Mutation++

**DISCLAIMER**: pyMutation is in experimental development stage.

This file is written in [GitHub Flavored Markdown](https://help.github.com/articles/github-flavored-markdown/), it can be comfortably viewed using the [Atom text editor](https://atom.io/) and its Markdown viewer plugin.

## Package contents

The `examples/python/` directory contains:

* the pyMutation setup script and module directory
* the pyMutation tests and examples suite

## Build

pyMutation is built using the CMake build system. It uses the [Boost.Python](http://www.boost.org/doc/libs/1_59_0/libs/python/doc/) Python bindings and the [Boost.NumPy](https://github.com/ndarray/Boost.NumPy) third-party library.

When configuring the build system, the `BUILD_PYTHON_WRAPPER` flag has to be set to `ON`. A convenient way to configure CMake variables is the use the `ccmake` command instead of `cmake`.

Once CMake is configured, build Mutation++ as usual. The pyMutation shared library object will be  written to the build file tree. It can then be copied to the pyMutation module directory by executing

    python setup.py build

pyMutation is then read for use: just copy the `examples/python/pyMutation/` directory to the location of the Python project where it should be used and import the module by including the following code at the beginning of your Python script:

    import pyMutation

## Install

Once pyMutation is built, it can be installed to the system very conveniently. `cd` to the `examples/python/` directory and execute the `setup.py` script:

    python setup.py install

This will install pyMutation to the `site-packages` directory. **It is strongly advised to use [virtual environments](http://docs.python-guide.org/en/latest/dev/virtualenvs/) to manage Python packages.** If you are not using a virtual environment, pyMutation will be installed system-wide and the aforementioned command will probably require root privileges to be executed. Otherwise, pyMutation will be installed in the current virtual environment.

Once this is done, pyMutation can be accessed from anywhere: it is no longer necessary to have the module in the same directory as the calling script.

Note that pyMutation is not installed alongside the rest of the Mutation++ library when issuing the `make install` command.

## Uninstall

If installed using the previous procedure, pyMutation can be uninstalled using the PIP uninstaller:

    pip uninstall pyMutation
