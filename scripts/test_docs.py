#!/usr/bin/env python

import argparse
import collections
import fnmatch
import os
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
        lightred = ' \033[91m'
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

class TestCase:
    def __init__(self, name, tags, comment, body, file_name, line):  
        self.name = name.strip()
        self.tags = re.findall(r'\[(.*?)\]', tags)
        self.comment = comment.strip()
        self.body = body
        self.file_name = os.path.basename(file_name.strip())
        self.line = line

    def print_case(self):
        print '------------'
        print 'file: ' + self.file_name
        print 'line: ' + str(self.line)
        if self.comment != None:
            print self.comment
        print 'TEST_CASE'
        print '('
        print '    "' + self.name + '",'
        print '    "' + " ".join(self.tags) + '"'
        print ')'
        print '{'
        print '    ' + self.body
        print '}'


def find_enclosed_text(string, start = 0, open = '{', close = '}'):
    pos = start = string.find(open, start) + 1
    level = 1
    while True:
        if string[pos] == close:
            level = level - 1
            if level == 0:
                break
        else:
            if string[pos] == open:
                level = level + 1
        pos = pos + 1
        
    return string[start:pos]
    
def extract_test_cases(file_name):
    """ Opens the given file and extracts all of the Catch TEST_CASEs in the file, if they
        exist.  
        
        :param file_name: The name of the file."""
    
    with open(file_name, 'r') as f:
        text = f.read()
    
    regex = re.compile(r'(?:/\*((?:\*(?!/)|[^*])*)\*/\n?)?' + r'[ \t]*TEST_CASE[^(]*\([^"]*"([^"]+)"[^,]*,[^"]*"([^"]+)"[^)]*\)')
    cases = []
    for m in regex.finditer(text):
        comment = m.group(1)
        if comment == None:
            comment = ""
        cases.append(
            TestCase(
                m.group(2), 
                m.group(3), 
                comment, 
                find_enclosed_text(text, m.end()), 
                file_name,
                text.count('\n',0,m.start()) + comment.count('\n') + 1
            )
        )
    
    return cases


def generate_output(filename, cases):
    """ Generates the output file given the test cases. """
    
    # Sort the tests
    cases.sort(key=lambda c: c.name.lower())
    
    # Create a dictionary of tests linked to tags
    tags = collections.defaultdict(list)
    for i, c in enumerate(cases):
        for t in c.tags:
            tags[t].append((i, c.name))
    
    with open(filename, 'w') as f:
        f.write('/**\n')
        f.write('\\page tests Test Cases\n')
        
        # List of tags
        f.write('\\section Tags\n')
        for i, t in enumerate(sorted(tags.keys(), key=lambda s: s.lower())):
            f.write('- \\subpage test_tag_' + str(i) + ' "' + t + '"\n')
        f.write('\n')
    
        # Global list of tests
        f.write('\\section Tests\n')
        for i, c in enumerate(cases):
            f.write('- [' + c.name + '](\\ref test_case_' + str(i) + ')\n')
        f.write('*/\n')
        f.write('\n')
        
        # Single page for each tag
        for i, t in enumerate(sorted(tags.keys(), key=lambda s: s.lower())):
            f.write('/**\n')
            f.write('\\page test_tag_' + str(i) + ' ' + t + '\n')
            for val in tags[t]:
                f.write('- \\subpage test_case_' + str(val[0]) + ' "' + val[1] +'"\n')
            f.write('*/\n')
            f.write('\n')
        
        # Single page for each test
        for i, c in enumerate(cases):
            f.write('/**\n')
            f.write('\\page test_case_' + str(i) + ' Test Case: ' + c.name + '\n')
            f.write(c.file_name + ':' + str(c.line) + '\n')
            if c.comment != "":
                f.write('\\section Description\n')
                f.write(c.comment + '\n')
            if c.body != "":
                f.write('\\section Code\n')
                f.write('\\code{.cpp}\n')
                f.write(c.body + '\n')
                f.write('\\endcode\n')
            f.write('*/\n')
            f.write('\n')
        
        


def main(args):
    """ Run the main tasks of the script. """
    
    path_list = args.path.split()
    ext_list = args.extensions.split()

    # Get files to extract test cases from
    matches = set([])
    for path in path_list:
        for ext in ext_list:
            for root, dirnames, filenames in os.walk(path):
                for filename in fnmatch.filter(filenames, '*{}'.format(ext)):
                    matches.add(os.path.join(root, filename))

    # Extract the test case info from each file
    all_cases = []
    for filename in matches:
        cases = extract_test_cases(filename)
        all_cases.extend(cases)
    
        if args.verbose:
            print colors.fg.green + filename + colors.reset
            for c in cases:
                print '    line ' + str(c.line) + ': ' + c.name
                print '    ' + ' '.join(c.tags)
                
    # Generate the output
    generate_output(args.output, all_cases)



#
# Main entry point.  Parse the command line options.
#
if __name__ == "__main__":
    default_extensions_list = ['.cpp']
    default_paths_list = ['tests']
    default_output = 'docs/manual/tests.h'

    parser = argparse.ArgumentParser(
        description="""Parse the Catch TEST_CASESs from all of the testing source files 
                       and generate the API documentation to be read by Doxygen.""")

    parser.add_argument("-e", "--extensions", default=" ".join(default_extensions_list),
                        help="File extension to match. Multiple extensions can be "
                             "specified using quotes. "
                             "Defaults to \"{}\"".format(" ".join(default_extensions_list)))
    parser.add_argument("-p", "--path", default=" ".join(default_paths_list),
                        help="List of paths where source files "
                             "should be inspected. Multiple paths can be "
                             "specified using quotes. "
                             "Defaults to \"{}\".".format(" ".join(default_paths_list)))
    parser.add_argument("-o", "--output", default=default_output,
                        help="File where the output is generated. "
                             "Defaults to \"{}\".".format(default_output))
    parser.add_argument("-v", "--verbose", action='store_true',
                        help="Print detailed report (verbose mode).")
    args = parser.parse_args()

    main(args)
    
    
    
    