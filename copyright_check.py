import os
import argparse
import fnmatch
import re


class colors:
    """Colors class providing easy access to ANSI escape sequences.

     The generic style formatters bold, underline, reverse, etc. are defined.
     A special reset formatter resets any applied style or coloring.

    Background and foreground colors are nested in dedicated bg and fg classes.
    """
    reset = '\033[0m'
    bold = '\033[01m'
    disable = '\033[02m'
    underline = '\033[04m'
    reverse = '\033[07m'
    strikethrough = '\033[09m'
    invisible = '\033[08m'

    class fg:
        black = '\033[30m'
        red = '\033[31m'
        green = '\033[32m'
        orange = '\033[33m'
        blue = '\033[34m'
        purple = '\033[35m'
        cyan = '\033[36m'
        lightgrey = '\033[37m'
        darkgrey = '\033[90m'
        lightred =' \033[91m'
        lightgreen = '\033[92m'
        yellow = '\033[93m'
        lightblue = '\033[94m'
        pink = '\033[95m'
        lightcyan = '\033[96m'

    class bg:
        black = '\033[40m'
        red = '\033[41m'
        green = '\033[42m'
        orange = '\033[43m'
        blue = '\033[44m'
        purple = '\033[45m'
        cyan = '\033[46m'
        lightgrey = '\033[47m'


# Constants defining defaults for the copyright information note
YEAR = "2014"
NAME = "von Karman Institute for Fluid Dynamics (VKI)"

# Constants defining strings used for verbose reports
CPINFO_PRESENT_STRING = \
    {True: colors.fg.green + "present" + colors.reset,
     False: colors.fg.red + "missing" + colors.reset}
CPINFO_CORRECTLY_FORMATTED_STRING = \
    {True: colors.fg.green + "correctly formatted" + colors.reset,
     False: colors.fg.red + "improperly formatted" + colors.reset}
CPINFO_UPTODATE_STRING = \
    {True: colors.fg.green + "up-to-date" + colors.reset,
     False: colors.fg.red + "outdated" + colors.reset}
STATEMENT_GOOD_STRING = \
    {True: colors.fg.green + "correct" + colors.reset,
     False: colors.fg.red + "incorrect" + colors.reset}


class FileInfo(object):
    def __init__(self, filename):
        self.filename = filename
        self.statement_good = False
        self.cpinfo_present = False
        self.cpinfo_correctly_formatted = False
        self.cpinfo_uptodate = False
        self.cpinfo_line = ''

    def check_copyright(self, year=YEAR, name=NAME):
        # Check copyright info
        with open(self.filename, 'r') as source_file:
            for source_line in source_file:
                # Check presence
                # TODO: this presence check may yield many false positives and should be deleted
                if re.match("(.*)(C|c)opyright(.*)", source_line):
                    self.cpinfo_present = True
                    self.cpinfo_line = source_line

                    # Check format
                    self.cpinfo_correctly_formatted = \
                        bool(re.match("(.*)Copyright [0-9]{4}(-[0-9]{4})* (.*)",
                                      source_line))

                    # Check date and author
                    self.cpinfo_uptodate = \
                        bool(
                            re.match(
                                "(.*)Copyright {} {}(.*)".format(re.escape(YEAR),
                                                                 re.escape(NAME)),
                                source_line))

    def check_statement(self):
        # Check copying statement (line by line because of comment symbols)
        # Load statement lines
        statement_line_list = []
        with open('copying_statement.txt', 'r') as statement_file:
            statement_line_list = [line.strip() for line in statement_file.readlines()]

        with open(self.filename, 'r') as source_file:
            i_statement_line = 0

            for source_line in source_file:
                source_line = source_line.strip(" * ").strip()

                statement_line = statement_line_list[i_statement_line]
                if statement_line == source_line:
                    i_statement_line += 1

                if i_statement_line >= len(statement_line_list):
                    self.statement_good = True
                    break

    def check(self):
        if      self.filename \
            and self.cpinfo_present \
            and self.cpinfo_correctly_formatted \
            and self.cpinfo_uptodate:
            return True
        else:
            return False

    def print_report(self, verbose=0):
        """Build a string reporting the status of the file for which
        info is stored. Lines after the 1st are indented with 4
        whitespaces."""
        report_string = colors.bold + self.filename + colors.reset

        if verbose:
            report_string += \
                "\n" + "Copying statement {}".format(
                    STATEMENT_GOOD_STRING[self.statement_good])

            report_string += \
                "\n" + "Copyright info {}".format(
                    CPINFO_PRESENT_STRING[self.cpinfo_present])

            if self.cpinfo_present:
                report_string += \
                    ", {}".format(
                        CPINFO_CORRECTLY_FORMATTED_STRING[self.cpinfo_correctly_formatted])

                if self.cpinfo_correctly_formatted:
                    report_string += \
                        ", {}".format(
                            CPINFO_UPTODATE_STRING[self.cpinfo_uptodate])

                report_string += \
                    "\n" + "{}".format(
                        self.cpinfo_line.strip(" * ").strip())

        print report_string.replace("\n", "\n    ")


def main(args):
    path_list = args.path.split()
    ext_list = args.extensions.split()

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

            source_file_info.check_copyright()
            source_file_info.check_statement()

            if not source_file_info.check():
                source_file_info.print_report(args.verbose)


if __name__ == "__main__":
    default_extensions_list = ['.h', '.cpp']
    default_paths_list = ['src', 'tests']

    parser = argparse.ArgumentParser(
        description="""Check if the copyright information and headers are present and 
                       appropriately formatted. By default, the program lists the names 
                       of files with (1) missing, improperly formatted or outdated copyright 
                       information or (2) incorrect or missing statement of copying 
                       permission.""")

    parser.add_argument("-e", "--extensions", default=" ".join(default_extensions_list),
                        help="File extension to match. Multiple extensions can be "
                             "specified using quotes. "
                             "Defaults to \"{}\"".format(" ".join(default_extensions_list)))
    parser.add_argument("-p", "--path", default=" ".join(default_paths_list),
                        help="List of paths where source files "
                             "should be inspected. Multiple paths can be "
                             "specified using quotes. "
                             "Defaults to \"{}\".".format(" ".join(default_paths_list)))
    parser.add_argument("-u", "--update", action="store_true",
                        help="Update copyright information if it is properly formatted.")
    parser.add_argument("-v", "--verbose", action='store_true',
                        help="Print detailed report (verbose mode).")
    args = parser.parse_args()

    main(args)
