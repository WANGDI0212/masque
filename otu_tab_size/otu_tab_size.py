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
from __future__ import print_function
import os
import sys
import argparse
import csv

__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2015, Institut Pasteur"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@pasteur.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      Arguments:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.add_argument('-i', dest='input_file', type=isfile, required=True,
                        help='Path to the count matrix file.')
    parser.add_argument('-g', dest='gene_file', type=isfile, required=True,
                        help='Path to the database fasta file.')
    parser.add_argument('-o', dest='output_file', type=str,
                            default=os.curdir + os.sep + "output.txt",
                        help='Output file.')
    args = parser.parse_args()
    return args


def parse_fasta(gene_file):
    """Compute gene length
    """
    gene_length = {}
    header = ''
    seq = ''
    try:
        with open(gene_file, "rt") as gene:
            for line in gene:
                if line.startswith(">"):
                    size = len(seq)
                    if size > 0 and header != '':
                        gene_length[header] = size
                        seq = ''
                    header = line[1:].strip()
                else:
                    seq += line.strip()
            size = len(seq)
            if size > 0 and header != '':
                gene_length[header] = size
            assert(len(gene_length) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(gene_file))
    except AssertionError:
        sys.exit("No gene parsed from {0}".format(gene_file))
    return gene_length


def write_result(input_file, gene_length, output_file):
    """Write the new count matrix with gene length
    """
    try:
        with open(input_file, "rt") as input_table:
            with open(output_file, "wt") as output_table:
                input_reader = csv.reader(input_table, delimiter='\t')
                output_writer = csv.writer(output_table, delimiter='\t')
                output_writer.writerow(input_reader.next() + ['size'])
                for line in input_reader:
                    if line[0] not in gene_length:
                        sys.exit("The sequence {0} is not present in the "
                                 "database".format(line[0]))
                    output_writer.writerow(line + [gene_length[line[0]]])
    except IOError:
        sys.exit("Error cannot open {0}".format(input_file))


def main():
    """Main program
    """
    args = get_arguments()
    gene_length = parse_fasta(args.gene_file)
    write_result(args.input_file, gene_length, args.output_file)

if __name__ == '__main__':
    main()