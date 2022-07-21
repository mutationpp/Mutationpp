import re
import sys
from pathlib import Path

from setuptools_scm import get_version

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

ROOT_CMAKELISTS = "./CMakeLists.txt"

DESCRIPTION = """An open-source library providing thermodynamic, transport,\
chemistry, and energy transfer properties associated with subsonic to\
hypersonic flows."""

version_python = get_version()
version_cmake = ".".join(version_python.split(".")[:3])


setup(
    name="mutationpp",
    version=version_python,
    description=DESCRIPTION,
    long_description=DESCRIPTION,
    author="James B. Scoggins",
    license="LGPL3",
    package_dir={"": "interface/python"},
    packages=["mutationpp"],
    extras_require={"test": ["numpy", "pytest"], "release": ["autopub"]},
    cmake_install_dir="interface/python/mutationpp",
    cmake_args=[f"-DMUTATIONPP_VERSION_FROM_PYTHON={version_cmake}"],
)
