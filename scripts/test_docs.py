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


import argparse
import collections
import fnmatch
import os
import re
from utils import colors


class TestCase(object):
    """ Simple class representing a single test case. 
    """
    
    def __init__(self, name, tags, comment, body, file_name, line):
        """ Default constructor.
        
        :param name: Name of the test case.
        :param tags: Tags this test case belongs to.
        :param comment: Comment associated with the test case.
        :param body: The body (or code) implementing the test case.
        :param file_name: The name of the file this test case was defined in.
        :param line: The line number this test case was defined at. 
        """
         
        self.name = name.strip()
        self.tags = re.findall(r'\[(.*?)\]', tags)
        self.comment = comment.strip()
        self.body = body
        self.file_name = os.path.basename(file_name.strip())
        self.line = line

    def print_case(self):
        """ Simply prints test case parameters to the standard output.  Useful for 
        debugging. 
        """
        
        print('------------')
        print('file: ' + self.file_name)
        print('line: ' + str(self.line))
        if self.comment is None:
            print(self.comment)
        print('TEST_CASE')
        print('(')
        print('    "' + self.name + '",')
        print('    "' + " ".join(self.tags) + '"')
        print(')')
        print('{')
        print('    ' + self.body)
        print('}')
    
    def generate_doxygen_page(self, link_name):
        """ Generates a doxygen style manual page for this test case.
        
        :param page_name: The name of the link associated with this page. 
        """
        
        # Header
        doxygen_page = '/**\n'
        doxygen_page += '\\page ' + link_name + ' Test Case: ' + self.name + '\n'
        doxygen_page += self.file_name + ':' + str(self.line) + '\n'
        
        # Comment
        if self.comment != "":
            doxygen_page += '\\section Description\n'
            doxygen_page += self.comment + '\n'
        
        # Body
        if self.body != "":
            doxygen_page += '\\section Code\n'
            doxygen_page += '\\code{.cpp}\n'
            doxygen_page += self.body + '\n'
            doxygen_page += '\\endcode\n'
        
        # Tail
        doxygen_page += '*/\n'
        
        return doxygen_page


def get_nested_block(string, start_index = 0, open_char = '{', close_char = '}'):
    """ Extracts the first text block in the string which is between the open and close 
    characters, accounting for nesting.
    
    :param string: The string to search.
    :param start: Starting location in the string.
    :param open: The character which opens the nested block.
    :param close: The character which closes the nested block. 
    """
    
    pos = start_index = string.find(open_char, start_index) + 1
    level = 1
    while True:
        if string[pos] == close_char:
            level = level - 1
            if level == 0:
                break
        else:
            if string[pos] == open_char:
                level = level + 1
        pos = pos + 1
        
    return string[start_index:pos]
    
def extract_test_cases(file_name):
    """ Opens the given file and extracts all of the Catch TEST_CASEs in the file, if they
    exist.  A list of the test cases in the file is returned.
        
    :param file_name: The name of the file.
    """
    
    with open(file_name, 'r') as f:
        text = f.read()
    
    regex = re.compile(
        # Comment
        r'(?:/\*((?:\*(?!/)|[^*])*)\*/\n?)?' + 
        # TEST_CASE structure
        r'[ \t]*TEST_CASE[^(]*\([^"]*"([^"]+)"[^,]*,[^"]*"([^"]+)"[^)]*\)')
        
    test_cases = []
    for m in regex.finditer(text):
        comment = m.group(1)
        if comment is None:
            comment = ""
        test_cases.append(
            TestCase(
                m.group(2), 
                m.group(3), 
                comment, 
                get_nested_block(text, m.end()), 
                file_name,
                text.count('\n',0,m.start()) + comment.count('\n') + 1
            )
        )
    
    return test_cases


def generate_doxygen_pages(file_name, test_cases):
    """ Generates the doxygen formatted manual pages related to the tests. 
    
    :param file_name: Name of the file where the output is generated.
    :param cases: List of TestCase objects. 
    """
    
    # Sort the tests
    test_cases.sort(key=lambda c: c.name.lower())
    
    # Create a dictionary of tests linked to tags
    tags = collections.defaultdict(list)
    for i, c in enumerate(test_cases):
        for t in c.tags:
            tags[t].append((i, c.name))
    
    with open(file_name, 'w') as f:
        f.write('/**\n')
        f.write('\\page tests Test Cases\n')
        
        # List of tags
        f.write('\\section Tags\n')
        for i, t in enumerate(sorted(list(tags.keys()), key=lambda s: s.lower())):
            f.write('- \\subpage test_tag_' + str(i) + ' "' + t + '"\n')
        f.write('\n')
    
        # Global list of tests
        f.write('\\section Tests\n')
        for i, c in enumerate(test_cases):
            f.write('- [' + c.name + '](\\ref test_case_' + str(i) + ')\n')
        f.write('*/\n')
        f.write('\n')
        
        # Single page for each tag
        for i, t in enumerate(sorted(list(tags.keys()), key=lambda s: s.lower())):
            f.write('/**\n')
            f.write('\\page test_tag_' + str(i) + ' ' + t + '\n')
            for val in tags[t]:
                f.write('- \\subpage test_case_' + str(val[0]) + ' "' + val[1] +'"\n')
            f.write('*/\n')
            f.write('\n')
        
        # Single page for each test
        for i, c in enumerate(test_cases):
            f.write(c.generate_doxygen_page('test_case_' + str(i)))
            f.write('\n')
        
        
def main(args):
    """ Run the main tasks of the script. 
    """
    
    path_list = args.path.split()
    ext_list = args.extensions.split()

    # Get files to extract test cases from
    matches = set([])
    for path in path_list:
        for ext in ext_list:
            for root, dirnames, file_names in os.walk(path):
                for file_name in fnmatch.filter(file_names, '*{}'.format(ext)):
                    matches.add(os.path.join(root, file_name))

    # Extract the test case info from each file
    all_cases = []
    for file_name in matches:
        cases = extract_test_cases(file_name)
        all_cases.extend(cases)
    
        if args.verbose:
            print(colors.fg.green + file_name + colors.reset)
            for c in cases:
                print('    line ' + str(c.line) + ': ' + c.name)
                print('    ' + ' '.join(c.tags))
                
    # Generate the output
    generate_doxygen_pages(args.output, all_cases)


# Entry point for the program.  Parse command line options and run script.
if __name__ == "__main__":    
    default_extensions_list = ['.cpp']
    default_paths_list = ['tests']
    default_output = 'docs/tests.h'

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
    
    
    
    