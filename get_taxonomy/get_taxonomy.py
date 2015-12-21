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
__copyright__ = "Copyright 2014, INRA"
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


def isdir(path):
    """Check if path is an existing file.
      Arguments:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file.".format(path)
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
    parser.add_argument('-i', dest='input_file', type=isfile, required=True,
                        help='Path to the input file.')
    parser.add_argument('-d', dest='database_file', type=isfile, required=True,
                        help='Path to the database file.')
    parser.add_argument('-dtype', dest='database_type', type=str,
                        default="silva", choices=["silva", "rdp",
                                                  "greengenes", "unite"],
                        help='Database format (default = silva).')
    parser.add_argument('-t', dest='taxonomy_file', type=isfile,
                        help='Path to the taxonomy file (Greengenes only).')
    parser.add_argument('-o', dest='output_file', type=str, default=None,
                        help='Output file.')
    parser.add_argument('-r', dest='results', type=isdir,
                        default=os.curdir + os.sep,
                        help='Path to result directory.')
    args = parser.parse_args()
    return args


def load_vsearch(input_file):
    """Load assignation provided with vsearch
    """
    vsearch_dict = {}
    try:
        with open(input_file, "rt") as input_data:
            input_reader = csv.reader(input_data, delimiter='\t')
            for line in input_reader:
                annotation = line[1].strip()
                if annotation in vsearch_dict:
                    vsearch_dict[annotation] += [[line[0], float(line[2])]]
                else:
                    vsearch_dict[annotation] = [[line[0], float(line[2])]]
        assert(len(vsearch_dict) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(input_file))
    except AssertionError:
        sys.exit("Nothing read from {0}".format(input_file))
    return vsearch_dict


def parse_rdp(header, vsearch_dict, annotation_dict):
    """Parse RDP annotation
    """
    identified = 0
    tax = header.strip().split(" ")[0][1:]
    lineage = header.strip().split("=")[1]
    if tax in vsearch_dict:
        lineage = lineage.split(";")[2:]
        lineage = [lineage[i].replace("\"","") for i in xrange(0, len(lineage), 2)]
        annotation_dict[tax] = ";".join(lineage)
        identified = 1
    return annotation_dict, identified


def parse_unite(header, vsearch_dict, annotation_dict):
    """Parse unite annotation
    """
    identified = 0
    tax = header.strip().split(" ")[0][1:]
    lineage = header.strip().split(" ")[1]
    if tax in vsearch_dict:
        lineage = lineage.split(";")
        lineage = [lineage[i].split("_")[2] for i in xrange(0, len(lineage))]
        annotation_dict[tax] = ";".join(lineage)
        identified = 1
    return annotation_dict, identified


def parse_silva(header, vsearch_dict, annotation_dict):
    """Parse silva database annotations
    """
    identified = 0
    tax = header.strip().split(" ")
    if tax[0][1:].strip() in vsearch_dict:
        annotation_dict[tax[0][1:].strip()] = tax[1]
        identified = 1
    return annotation_dict, identified


def load_taxonomy_gg(taxonomy_file, vsearch_dict):
    """Load greengenes taxonomy file
    """
    annotation_dict = {}
    try:
        with open(taxonomy_file, "rt") as taxonomy:
            taxonomy_reader = csv.reader(taxonomy, delimiter="\t")
            for line in taxonomy_reader:
                annotation_dict[line[0]] = "".join([annot[3:]
                                                    for annot in line[1].split(" ")])
    except IOError:
        sys.exit("Error cannot open {0}".format(taxonomy_file))
    return annotation_dict


def load_taxonomy(database_file, vsearch_dict, database_type):
    """Load rdp and silva annotations
    """
    annotation_dict = {}
    if database_type == "rdp":
        parse_result = parse_rdp
    elif database_type == "unite":
        parse_result = parse_unite
    else:
        parse_result = parse_silva
    nb_id = len(vsearch_dict)
    count_identified = 0
    #print(vsearch_dict)
    try:
        with open(database_file, "rt") as database:
            for line in database:
                if line.startswith(">"):
                    annotation_dict, count = parse_result(line, vsearch_dict,
                                                      annotation_dict)
                    count_identified += count
                    if count_identified == nb_id:
                        break
    except IOError:
        sys.exit("Error cannot open {0}".format(database_file))
    return annotation_dict


def write_tax_table(vsearch_dict, annotation_dict, output_file):
    """
    """
    # Identity threshold :
    # Uniting the classification of cultured and uncultured bacteria and archaea using 16S rRNA gene sequences
    # Pablo Yarza,       Pelin Yilmaz,   Elmar Pruesse,  Frank Oliver Glöckner,  Wolfgang Ludwig,        Karl-Heinz Schleifer,   William B. Whitman,     Jean Euzéby,    Rudolf Amann    & Ramon Rosselló-Móra
    # Nature Reviews Microbiology 12, 635–645 (2014) doi:10.1038/nrmicro3330
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter='\t')
            output_writer.writerow(["OTU", "Kingdom", "Phylum", "Class",
                                    "Order", "Family", "Genus", "Specie"])
            for tax in vsearch_dict:
                for OTU in vsearch_dict[tax]:
                    taxonomy = annotation_dict[tax].split(";")
                    # Genus and the rest
                    if OTU[1] >= 94.5:
                        pass
                        #taxonomy = taxonomy[:-1]
                    # Family
                    elif OTU[1] < 94.5 and OTU[1] >= 86.5:
                        taxonomy = taxonomy[:-2]
                    # Order
                    elif OTU[1] < 86.5 and OTU[1] >= 82.0:
                        taxonomy = taxonomy[:-3]
                    # Class
                    elif OTU[1] < 82.0 and OTU[1] >= 78.5:
                        taxonomy = taxonomy[:-4]
                    # Phylum
                    elif OTU[1] < 78.5 and OTU[1] >= 75.0:
                        taxonomy = taxonomy[:-5]
                    else:
                        sys.exit("Strange id is to low for {0[0]} : {0[1]} \%".format(OTU))
                    output_writer.writerow([OTU[0]] + taxonomy)
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


def main():
    """Main program
    """
    args = getArguments()
    vsearch_dict = load_vsearch(args.input_file)
    if args.database_type == "greengenes" and args.taxonomy_file:
        annotation_dict = load_taxonomy_gg(args.taxonomy_file, vsearch_dict)
    elif args.database_type == "greengenes" and not args.taxonomy_file:
        sys.exit("Please indicate the greengenes taxonomy file")
    else:
        annotation_dict = load_taxonomy(args.database_file, vsearch_dict,
                                        args.database_type)
    write_tax_table(vsearch_dict, annotation_dict, args.output_file)


if __name__ == '__main__':
    main()
