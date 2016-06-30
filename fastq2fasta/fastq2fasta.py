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


def getArguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile, required=True,
                        help='Path to the fastq file.')
    parser.add_argument('-s', dest='sample_name', type=str, default="",
                        help='Add sample to the reads for usearch pipeline.')
    parser.add_argument('-o', dest='output_file', type=str, default=None,
                        help='Output file.')
    args = parser.parse_args()
    return args


def convert_fastq_fasta(fastq_file, sample_name, output_file):
    """
    """
    if not output_file:
        output = sys.stdout
    else:
        output = open(output_file, "wt")
    try:
        with open(fastq_file, "rt") as fastq:
            for line in fastq:
                header = line[1:].split(" ")[0]
                line = fastq.next()
                if sample_name:
                    print(">{0};barcodelabel={2}\n{1}".format(
                            header[1:].replace("\n", ""), line.replace("\n", ""), sample_name),
                            file=output)
                else:
                    print(">{0}\n{1}".format(header[1:], line.replace("\n", "")),
                          file=output)
                fastq.next()
                fastq.next()
    except IOError:
        sys.exit("Error cannot open {0}".format(fastq_file))
    if output_file:
        output.close()


def main():
    """Main program
    """
    args = getArguments()
    convert_fastq_fasta(args.fastq_file, args.sample_name, args.output_file)


if __name__ == '__main__':
    main()
