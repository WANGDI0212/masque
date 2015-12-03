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
                        help='Path to the representant swarm file.')
    parser.add_argument('-c', dest='clustering_file', type=isfile,
                        required=True,
                        help='Path to swarm clustering file.')
    parser.add_argument('-u', dest='uclust_file', type=isfile,
                        help='Path to uclust file.')
    parser.add_argument('-o', dest='output_file', type=str, default=None,
                        help='Output file.')
    parser.add_argument('-oc', dest='output_clustering_file', type=str,
                        default=None, help='Output clustering file.')
    parser.add_argument('-ou', dest='output_uclust_file', type=str,
                        default=None, help='Output uclust file.')
    parser.add_argument('-r', dest='results', type=isdir,
                        default=os.curdir + os.sep,
                        help='Path to result directory.')
    args = parser.parse_args()
    return args


def get_cluster(clustering_file):
    """
    """
    cluster_dict = {}
    try:
        with open(clustering_file, "rt") as clustering:
            clustering_reader = csv.reader(clustering, delimiter=' ')
            for i,cluster in enumerate(clustering_reader):
                cluster_dict[";".join(cluster[0].split(";")[:-2])] = (
                    ["OTU_{0}".format(i + 1),
                     [";".join(value.split(";")[:-2]) for value in cluster]])
        assert(len(cluster_dict) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(clustering_file))
    except AssertionError:
        sys.exit("Error no element read from {0}".format(clustering_file))
    return cluster_dict


def set_cluster(cluster_dict, output_clustering_file):
    """
    """
    try:
        with open(output_clustering_file, "wt") as output_clustering:
            output_clustering_writer = csv.writer(output_clustering,
                                                  delimiter='\t')
            output_clustering_writer.writerow(["OTU", "OTU_representant",
                                               "OTU_cluster"])
            for cluster in cluster_dict:
                output_clustering_writer.writerow([cluster_dict[cluster][0],
                                                   cluster, " ".join(cluster_dict[cluster][1])])
    except IOError:
        sys.exit("Error cannot open {0}".format(output_clustering_file))


def convert_swarm_fasta(input_file, cluster_dict, output_file):
    """
    """
    clust_header = ""
    sequence = ""
    if not output_file:
        output = sys.stdout
    else:
        output = open(output_file, "wt")
    try:
        with open(input_file, "rt") as fasta:
            for line in fasta:
                if line.startswith(">"):
                    if clust_header:
                        print(">{0}\n{1}".format(
                                clust_header,
                                sequence.upper().replace("\n","")),
                              file=output)
                    header = ";".join(line[1:].replace("\n","").split(";")[:-2])
                    clust_header = cluster_dict[header][0]
                    sequence = ""
                elif len(line) > 0 and header:
                    sequence += line
            if clust_header and sequence:
                print(">{0}\n{1}".format(clust_header,
                                         sequence.upper().replace("\n","")),
                        file=output)
    except IOError:
        sys.exit("Error cannot open {0}".format(input_file))


def convert_swarm_uclust(uclust_file, cluster_dict, output_uclust_file):
    """
    """
    if not output_uclust_file:
        output = sys.stdout
    else:
        output = open(output_uclust_file, "wt")
    try:
        with open(uclust_file, "rt") as uclust:
            uclust_reader = csv.reader(uclust, delimiter="\t")
            for line in uclust_reader:
                if line[9] != "*":
                    line[9] = cluster_dict[";".join(line[9].split(";")[:-2])][0]
                print("\t".join(line), file=output)
    except IOError:
        sys.exit("Error cannot open {0}".format(uclust_file))


def main():
    """Main program
    """
    args = getArguments()
    # Read clustering
    cluster_dict = get_cluster(args.clustering_file)
    # Assign OTU number to each cluster
    set_cluster(cluster_dict, args.output_clustering_file)
    #print(cluster_dict)
    # Reassign fasta file
    convert_swarm_fasta(args.input_file, cluster_dict, args.output_file)
    # Reassign uclust file
    if args.uclust_file:
        convert_swarm_uclust(args.uclust_file, cluster_dict, args.output_uclust_file)


if __name__ == '__main__':
    main()