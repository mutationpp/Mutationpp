#!/usr/bin/python

from __future__ import print_function
import argparse
import multiprocessing
import os
import subprocess

from utils import colors, reindent, makedirs_exist_ok


class Dependency(object):
    """Default installation method for packages using CMake as their build system"""

    def __init__(self, name, build_cmd_seq, build_dir=None, source_dir=None):
        self.name = name
        self.source_dir = os.path.join(MPP_THIRDPARTY_DIRECTORY, self.name) if source_dir is None else source_dir
        self.build_dir = os.path.join(self.source_dir, "build") if build_dir is None else build_dir
        self.build_cmd_seq = build_cmd_seq

    def install(self):

        print(colors.fg.red + "Installing " + self.name + colors.reset)
        os.chdir(self.source_dir)
        makedirs_exist_ok(self.build_dir)
        os.chdir(self.build_dir)

        for cmd in self.build_cmd_seq:
            print(colors.fg.blue + " ".join(cmd) + colors.reset)
            try:
                subprocess.call(cmd)
            except OSError as e:
                print(colors.fg.red +
                      "Command call raised OSError: {} was not found".format(cmd[0]) +
                      colors.reset)
                raise e


# Global variables
# TODO: check if environment variables are defined
MPP_THIRDPARTY_DIRECTORY = os.path.join(os.environ["MPP_DIRECTORY"], "thirdparty")
MPP_INSTALL_DIRECTORY = os.path.join(os.environ["MPP_DIRECTORY"], "install")
NUM_PROC = multiprocessing.cpu_count()


def generate_deps_dict(cmd_dict):
    return \
        {
            "eigen": Dependency(
                name="eigen",
                build_cmd_seq=[[cmd_dict["cmake"],
                                "-DCMAKE_INSTALL_PREFIX={}".format(MPP_INSTALL_DIRECTORY),
                                ".."],
                               [cmd_dict["make"],
                                "-j", str(NUM_PROC),
                                "install"]]),
            "catch": Dependency(
                name="catch",
                build_cmd_seq=[[cmd_dict["cmake"],
                                "-DCMAKE_INSTALL_PREFIX={}".format(MPP_INSTALL_DIRECTORY),
                                ".."],
                               [cmd_dict["make"],
                                "-j", str(NUM_PROC),
                                "install"]]),
        }


def main(args):

    cmd_dict = {"cmake": args.cmd_cmake,
                "make": args.cmd_make}

    deps_dict = generate_deps_dict(cmd_dict)
    deps_list = sorted(deps_dict.keys())

    if args.list:
        print("Dependencies available for installation:\n" +
              reindent("\n".join(deps_list), 4))
        return

    if not args.deps:  # Means that no positional argument is provided
        for name in deps_list:
            deps_dict[name].install()
    else:
        for name in args.deps:
            deps_dict[name].install()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="""Install third-party dependencies for the Mutation++ software package.""")
    parser.add_argument("-l", "--list", action="store_true",
                        help="show the list of available dependencies")
    parser.add_argument('deps', nargs='*', default=None,
                        help="list of dependencies to install (if none is specified, install all dependencies)")
    parser.add_argument("--cmd_cmake", default="cmake",
                        help="command to use to call cmake (defaults to cmake)")
    parser.add_argument("--cmd_make", default="make",
                        help="command to use to call make (defaults to make)")
    args = parser.parse_args()

    main(args)
