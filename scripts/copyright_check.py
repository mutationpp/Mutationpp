from __future__ import print_function
# Copyright 2017-2020 von Karman Institute for Fluid Dynamics (VKI)
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


# TODO: single-file check


import os
import argparse
import fnmatch
import re
import datetime
import shutil
from utils import colors


# Constants defining defaults for the copyright information and the copying statement
YEAR = datetime.date.today().year
NAME = "von Karman Institute for Fluid Dynamics (VKI)"
STATEMENT_TEXT = \
"""This file is part of MUlticomponent Thermodynamic And Transport
properties for IONized gases in C++ (Mutation++) software package.

Mutation++ is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

Mutation++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with Mutation++.  If not, see
<http://www.gnu.org/licenses/>."""

# Constants defining strings used for verbose reports
CPINFO_CORRECTLY_FORMATTED_STRING = \
    {True: colors.fg.green + "correctly formatted" + colors.reset,
     False: colors.fg.red + "incorrect or missing" + colors.reset}
CPINFO_DATE_CORRECT_STRING = \
    {True: colors.fg.green + "correct" + colors.reset,
     False: colors.fg.red + "incorrect" + colors.reset}
STATEMENT_CORRECT_STRING = \
    {True: colors.fg.green + "correct" + colors.reset,
     False: colors.fg.red + "incorrect" + colors.reset}


class FileInfo(object):
    def __init__(self, filename, statement_text):
        """Default constructor.

        :param filename: Path to the file being checked.
        :param statement_text: Text for the copying statement."""

        self.filename = filename
        self.statement_text = statement_text
        self.statement_correct = None
        self.cpinfo_correctly_formatted = None
        self.cpinfo_date_correct = None
        self.cpinfo_line = None
        self.cpinfo_line_corrected = None

    def check_cpinfo_correctly_formatted(self, name=NAME, force=False):
        """Checks if the file contains a copyright information line and if it is properly formatted,
        i.e. if it has the format `Copyright <year> <name>` or `Copyright <year1>-<year2> <name>`.
        If the update has already been run, it is not run again (unless forced).

        :param name: Author name to be looked up. Defaults to the NAME global variable.
        :param force: Force update."""

        if self.cpinfo_correctly_formatted is None or force:
            # Check copyright info format
            with open(self.filename, 'r') as source_file:
                for source_line in source_file:
                    # Check format
                    search_obj = re.search("Copyright ([0-9]{4})-?([0-9]{4})? " + re.escape(name),
                                           source_line)
                    if search_obj is not None:
                        self.cpinfo_correctly_formatted = True
                        self.cpinfo_line = search_obj.group(0)
                        return

            # If no match was found, it means that no correctly formatted copyright info was found
            self.cpinfo_correctly_formatted = False

    def check_cpinfo_date_correct(self, year=YEAR, force=False):
        """If the file contains a correctly formatted copyright information line, checks if the
        date is up-to-date. If not, propose a correction. This method automatically ensures that
        the copyright information line format has been checked for and forces a check if necessary.
        If the update has already been run, it is not run again (unless forced).

        :param year: Current year. Defaults to the current year.
        :param force: Force update. This will not force self.check_cpinfo_correctly_formatted()."""

        # First, make sure that a copyright info line was searched for
        self.check_cpinfo_correctly_formatted()

        if self.cpinfo_date_correct is None or force:
            if self.cpinfo_correctly_formatted:
                # Check date
                if self.cpinfo_correctly_formatted:
                    # Search and store year info
                    search_obj = re.search("(?P<first>[0-9]{4})-?(?P<second>[0-9]{4})?",
                                           self.cpinfo_line)

                    # Check if year is up-to-date; if not, propose a correction
                    if search_obj is not None:
                        first = int(search_obj.groupdict()['first'])
                        second = int(search_obj.groupdict()['second']) \
                            if search_obj.groupdict()['second'] is not None \
                            else None

                        if second is not None:
                            if second == year:
                                self.cpinfo_date_correct = True
                                self.cpinfo_line_corrected = self.cpinfo_line
                            else:
                                self.cpinfo_date_correct = False
                                self.cpinfo_line_corrected = re.sub(str(second), str(year), self.cpinfo_line)
                        else:
                            if first == year:
                                self.cpinfo_date_correct = True
                                self.cpinfo_line_corrected = self.cpinfo_line
                            else:
                                self.cpinfo_date_correct = False
                                self.cpinfo_line_corrected = re.sub(str(first), "{}-{}".format(first, year), self.cpinfo_line)

    def check_statement_correct(self, force=False):
        """Check that the copying statement is correctly formatted. This method work with every comment types,
        but it requires the paging (newline positioning) and the case to be correct. If not, false
        negatives or positives are possible.
        If the update has already been run, it is not run again (unless forced).

        :param force: Force update."""

        if self.statement_correct is None or force:
            # Check copying statement (line by line because of comment symbols)
            statement_line_list = [re.escape(line.strip()) for line in self.statement_text.split("\n")]

            # Build the regex pattern
            pattern = r"\s+.*".join(statement_line_list)

            # Match it with file contents
            with open(self.filename, 'r') as source_file:
                if re.search(pattern, source_file.read()):
                    self.statement_correct = True
                else:
                    self.statement_correct = False

    def check(self):
        self.check_cpinfo_correctly_formatted()
        self.check_cpinfo_date_correct()
        self.check_statement_correct()

        if self.filename \
                and self.cpinfo_correctly_formatted \
                and self.cpinfo_date_correct \
                and self.statement_correct:
            return True
        else:
            return False

    def print_report(self, verbose=False):
        """Report if the file has problems.

        :param verbose: Activate verbose mode and print a detailed report. Lines after
                        the 1st are indented with 4 whitespaces."""

        report_string = colors.bold + self.filename + colors.reset

        if self.check():
            report_string += \
                ": " + colors.fg.green + "OK" + colors.reset
        else:
            report_string += \
                ": " + colors.fg.red + "NOT OK" + colors.reset

        if verbose:
            report_string += \
                "\nCopyright info {}".format(
                    CPINFO_CORRECTLY_FORMATTED_STRING[self.cpinfo_correctly_formatted])

            if self.cpinfo_correctly_formatted:
                report_string += \
                    "\nCopyright info date " + \
                    "{}".format(CPINFO_DATE_CORRECT_STRING[self.cpinfo_date_correct])

                if not self.cpinfo_date_correct:
                    report_string += \
                        "\nCurrent copyright info: " + \
                        colors.fg.red + \
                        "{}".format(self.cpinfo_line) + \
                        colors.reset
                    report_string += \
                        "\nCorrected copyright info: " + \
                        colors.fg.green + \
                        "{}".format(self.cpinfo_line_corrected) + \
                        colors.reset

            report_string += \
                "\nCopying statement " + \
                "{}".format(STATEMENT_CORRECT_STRING[self.statement_correct])

        print(report_string.replace("\n", "\n    "))

    def update_cpinfo_date(self, backup=False):
        """Update copyright info line if the date is incorrect.

        :param backup: Create backup of the file if modified. If a backup already exists, it is overwritten."""

        self.check_cpinfo_date_correct()

        if self.cpinfo_correctly_formatted and not self.cpinfo_date_correct:
            # Search for the copyright info line and replace it with the correct one
            source_data = None
            with open(self.filename, 'r') as source_file:
                source_data = source_file.read()

            source_data = source_data.replace(self.cpinfo_line, self.cpinfo_line_corrected)

            if backup:
                shutil.copy(self.filename, self.filename + ".bak")

            with open(self.filename, 'w') as source_file:
                source_file.write(source_data)


def main(args):
    path_list = args.path.split()
    ext_list = args.extensions.split()

    statement_text = STATEMENT_TEXT
    if args.statement is not None:
        with open(args.statement, 'r') as statement_file:
            statement_text = statement_file.read()

    # List files to be checked
    matches = set([])
    for path in path_list:
        for ext in ext_list:
            for root, dirnames, filenames in os.walk(path):
                for filename in fnmatch.filter(filenames, '*{}'.format(ext)):
                    matches.add(os.path.join(root, filename))

    # Check copyright information and correct files if requested
    for filename in matches:
            source_file_info = FileInfo(filename, statement_text)
            source_file_info.check()

            if args.faulty_only:
                if not source_file_info.check():
                    source_file_info.print_report(args.verbose)
            else:
                source_file_info.print_report(args.verbose)

            if args.update \
                    and source_file_info.cpinfo_correctly_formatted \
                    and not source_file_info.cpinfo_date_correct:
                source_file_info.update_cpinfo_date(backup=args.backup)


if __name__ == "__main__":
    default_extensions_list = ['.h', '.cpp', '.py', '.f90', 'CMakeLists.txt']
    default_paths_list = ['scripts', 'src', 'tests']

    parser = argparse.ArgumentParser(
        description="""Check if the copyright information and headers are present and
                       appropriately formatted. The program reports the files with
                       (1) missing, improperly formatted or outdated copyright information or
                       (2) incorrect or missing statement of copying permission. A default
                       copying statement text is provided, and it can be overridden with the
                       --statement option.""")

    parser.add_argument("-b", "--backup", action="store_true",
                        help="Create backup if files are updated.")
    parser.add_argument("-e", "--extensions", default=" ".join(default_extensions_list),
                        help="File extension to match. Multiple extensions can be "
                             "specified using quotes."
                             "Defaults to \"{}\"".format(" ".join(default_extensions_list)))
    parser.add_argument("-f", "--faulty-only", action="store_true",
                        help="Only show faulty files when reporting.")
    parser.add_argument("-p", "--path", default=" ".join(default_paths_list),
                        help="List of paths where source files "
                             "should be inspected. Multiple paths can be "
                             "specified using quotes. "
                             "Defaults to \"{}\".".format(" ".join(default_paths_list)))
    parser.add_argument("-s", "--statement", default=None,
                        help="File from which the copying statement text is loaded.")
    parser.add_argument("-u", "--update", action="store_true",
                        help="Update copyright information if it is properly formatted.")
    parser.add_argument("-v", "--verbose", action='store_true',
                        help="Print detailed report (verbose mode).")
    args = parser.parse_args()

    main(args)
