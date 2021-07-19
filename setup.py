import re
import sys

from pathlib import Path

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


def get_version_from_cmake(root_cmakelists):
    """A helper function to parse the root CMakeLists.txt and retrieve the
    version from there in order to have a single point holding the version
    info, easier to update"""

    pattern = re.compile(r" *VERSION *(\d+\.\d+\.\d+)")

    root_file = Path(root_cmakelists)

    with open(root_file, "r") as f:
        for line in f:
            match = pattern.search(line)
            if match is not None:
                return match.group(1)

    raise RuntimeError(
        "Couldn't parse CMakeLists.txt file to find a version statement"
    )


setup(
    name="mutationpp",
    version=get_version_from_cmake(ROOT_CMAKELISTS),
    description=DESCRIPTION,
    long_description=DESCRIPTION,
    author="James B. Scoggins",
    license="LGPL3",
    package_dir={"": "interface/python"},
    packages=["mutationpp"],
    extras_require={
        "test": ["numpy", "pytest"],
    },
    cmake_install_dir="interface/python/mutationpp",
)
