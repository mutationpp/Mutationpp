# Copyright 2017-2018 von Karman Institute for Fluid Dynamics (VKI)
#
# This file is part of MUlticomponent Thermodynamic And Transport
# properties for IONized gases in C++ (Mutation++) software package.
#
# Mutation++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Mutation++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with Mutation++.  If not, see
# <http://www.gnu.org/licenses/>.


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


def generate_deps_dict(cmd_dict):
    return \
        {
            "eigen": Dependency(
                name="eigen",
                build_cmd_seq=[[cmd_dict["cmake"],
                                "-DCMAKE_INSTALL_PREFIX={}".format(MPP_INSTALL_DIRECTORY),
                                ".."],
                               [cmd_dict["make"],
                                "-j", str(args.jobs),
                                "install"]]),
            "catch": Dependency(
                name="catch",
                build_cmd_seq=[[cmd_dict["cmake"],
                                "-DCMAKE_INSTALL_PREFIX={}".format(MPP_INSTALL_DIRECTORY),
                                ".."],
                               [cmd_dict["make"],
                                "-j", str(args.jobs),
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
    parser.add_argument('deps', nargs='*', default=None,
                        help="list of dependencies to install (if none is specified, install all dependencies)")
    parser.add_argument("-l", "--list", action="store_true",
                        help="show the list of available dependencies")
    parser.add_argument("-j", "--jobs", default=multiprocessing.cpu_count(),
                        help="allow JOBS jobs at once (defaults to the number of CPU units)")
    parser.add_argument("--cmd_cmake", default="cmake",
                        help="command to use to call cmake (defaults to cmake)")
    parser.add_argument("--cmd_make", default="make",
                        help="command to use to call make (defaults to make)")
    args = parser.parse_args()

    main(args)
