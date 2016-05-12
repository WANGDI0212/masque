#!/usr/bin/env python
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Extract the annotation in case of blast on imomi database."""

from __future__ import print_function
import os
import sys
import argparse
import csv

__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2014, INRA"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@jouy.inra.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


#===================
# parameters
#===================
def get_arguments():
    """Extract program options
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     usage="{0} -h [options] [arg]"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='annotation_list_file',
                        type=isfile, required=True, nargs='+',
                        help='Input several annotation result.')
    parser.add_argument('-o', dest='output_file',
                        type=str, default=os.curdir + os.sep + 'output.txt',
                        help='Output annotation.')
    return parser.parse_args()


def load_annoation(annotation_file, annotation_dict):
    """
    """
    try:
        with open(annotation_file, "rt") as annotation:
            annotation_reader = csv.reader(annotation, delimiter="\t")
            header = annotation_reader.next()
            for line in annotation_reader:
                #print(line)
                line = filter(None, line)
                if line[0] in annotation_dict:
                    if len(line[1:]) > len(annotation_dict[line[0]]):
                        #print("before")
                        #print(annotation_dict[line[0]])
                        #print(line[1:])
                        #print("after")
                        assert(len(line[1:]) <= 7)
                        annotation_dict[line[0]] = line[1:]
                else:
                    annotation_dict[line[0]] = line[1:]
    except IOError:
        sys.exit("Error cannot open {0}".format(annotation_file))
    return header, annotation_dict


def write_result(header, annotation_dict, output_file):
    """
    """
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter='\t')
            output_writer.writerow(header)
            for otu in annotation_dict:
                output_writer.writerow([otu] + annotation_dict[otu])
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


#===================
# MAIN
#===================
def main():
    """Main program
    """
    args = get_arguments()
    annotation_dict = {}
    for annotation_file in args.annotation_list_file:
        print(annotation_file)
        header, annotation_dict = load_annoation(annotation_file, annotation_dict)
        print(len(header))
    write_result(header, annotation_dict, args.output_file)
    
    
if __name__ == "__main__":
    main()
# END