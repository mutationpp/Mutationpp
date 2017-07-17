import os
import argparse
import fnmatch
import re

YEAR = "2014"
NAME = "von Karman Institute for Fluid Dynamics (VKI)"

CPINFO_PRESENT_STRING = {True: "present", False: "not present"}
CPINFO_CORRECT_STRING = {True: "correct", False: "not correct"}


class FileInfo(object):
    def __init__(self, filename):
        self.filename = filename
        self.cpinfo_present = False
        self.cpinfo_correct = False
        self.cpinfo_line = ''

    def analyze(self, year=YEAR, name=NAME):
        with open(self.filename, 'r') as source_file:

            for source_line in source_file:
                if re.match("(.*)Copyright(.*)", source_line):
                    self.cpinfo_present = True
                    self.cpinfo_line = source_line
                    self.cpinfo_correct = \
                        bool(re.match("(.*)Copyright\ {}\ {}(.*)".format(
                            re.escape(YEAR), re.escape(NAME)),
                            source_line))

    def check(self):
        if self.filename and self.cpinfo_present and self.cpinfo_correct:
            return True
        else:
            return False

    def print_report(self, verbose=0):
        """Build a string reporting the status of the file for which
        info is stored. Lines after the 1st are indented with 4
        whitespaces."""
        report_string = self.filename

        if verbose:
            report_string += \
                "\n" + "Copyright info {}".format(
                    CPINFO_PRESENT_STRING[self.cpinfo_present])

            if self.cpinfo_present:
                report_string += \
                    "\n" + "Copyright info {}".format(
                        CPINFO_CORRECT_STRING[self.cpinfo_correct])
                report_string += \
                    "\n" + "{}".format(
                        self.cpinfo_line.strip(" *").strip())

        print report_string.replace("\n", "\n    ")


def main(args):
    path_list = args.path.split() if args.path is not None else ['src', 'tests']

    ext_list = ['.h', '.cpp']

    # List files to be checked
    matches = []
    for path in path_list:
        for ext in ext_list:
            for root, dirnames, filenames in os.walk(path):
                for filename in fnmatch.filter(filenames, '*{}'.format(ext)):
                    matches.append(os.path.join(root, filename))

    # Check copyright information
    for filename in matches:
            source_file_info = FileInfo(filename)
            source_file_info.analyze()
            if not source_file_info.check():
                source_file_info.print_report(args.verbose)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path",
                        help="List of paths where source files "
                             "should be inspected (multiple paths can be "
                             "specified using quotes.")
    parser.add_argument("-v", "--verbose", action='store_true',
                        help="Print detailed report (verbose mode.")
    args = parser.parse_args()

    main(args)
