#!/usr/bin/env python3

import re
import argparse
import sys 
import errno
import shutil
import os
from os import listdir
from os.path import isfile, join

__author__ = "Katharina J. Hoff"
__copyright__ = "Copyright 2020. All rights reserved."
__license__ = "Artistic License"
__version__ = "1.0.0"
__credits__= "Anica Hoppe"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "production"

parser = argparse.ArgumentParser(
    description='BUSCO training genes have local gene IDs This script changes IDs so that ' +
                'they are non-redundant over different gff files and concatenates them into ' +
                'a new file; only CDS entries are transferred.')
parser.add_argument('-d', '--input_directory', required=True, type=str,
                    help='Directory that contains gff files produced by BUSCO')
parser.add_argument('-o', '--outfile', required=True, type=str,
                    help='Output file with non-redundant gene IDs.')
args = parser.parse_args()

''' Check whether args.input_directory exists and is readable '''

if not os.path.isdir(args.input_directory):
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': specified argument for --input_directory ' +
          "(" + args.input_directory + ") is not a directory!")
    exit(1)
elif not os.access(args.input_directory, os.R_OK):
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': specified argument for --input_directory ' +
          "(" + args.input_directory + ") is not readable!")
    exit(1)

''' Retrieve all files in input_directory, careful, all files are used '''

input_files = [f for f in listdir(args.input_directory) if isfile(join(args.input_directory, f))]

try:
    with open(args.outfile, 'w') as out_handle:
        for file in input_files:
            bid = re.search(r'([^.]+)\.gff', file).group(1)
            try:
                with open(join(args.input_directory, file), "r") as gff_handle:
                    for line in gff_handle:
                        if re.search(r'\tCDS\t', line):
                            new_line = re.sub(r'g(\d+)', r'g\1_' + re.escape(bid), line)
                            out_handle.write(new_line)
            except IOError:
                print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': Failed to open file ' + file + ' for reading!')
                quit(1)
except IOError:
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': Failed to open file ' + args.outfile + ' for writing!')
    quit(1)
