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

"""Extract 16S information."""
from __future__ import print_function
import argparse
import sys
import os
import csv
import glob
import re

__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2015, Institut Pasteur"
__credits__ = ["Amine Ghozlane"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@pasteur.fr"
__status__ = "Developpement"


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        # If false then its a list
        if not isinstance(values, basestring):
            out = []
            for val in values:
                if os.path.isfile(val):
                    out += [os.path.abspath(os.path.expanduser(val))]
                elif os.path.isdir(val):
                    out += [os.path.abspath(os.path.expanduser(val)) + os.sep]
                else:
                    out += [val]
            setattr(namespace, self.dest, out)
        # Value is a string
        else:
            if os.path.isfile(values):
                setattr(namespace, self.dest,
                        os.path.abspath(os.path.expanduser(values)))
            elif os.path.isdir(values):
                setattr(namespace, self.dest,
                        os.path.abspath(os.path.expanduser(values)) + os.sep)


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
            msg = "{0} is a file".format(path)
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
    parser.set_defaults(results=".{0}".format(os.sep))
    parser.add_argument('-d', dest='data_dir', type=isdir, action=FullPaths,
                        required=True, help='Result directory of MASQUE.')
    parser.add_argument('-r', dest='raw_reads_dir', type=isdir,
                        action=FullPaths, help='Raw reads directory.')
    parser.add_argument('-a', dest='amplicon_dir', type=isdir, action=FullPaths,
                        help='Amplicon directory (if raw reads not available).')
    parser.add_argument('-p', dest='paired_reads', action='store_true',
                        default=False, help='Paired reads (Default:False).')
    parser.add_argument('-o1', dest='output_file1', type=str,
                        default="masque_build_process.tsv",
                        help='Output file about the OTU building.')
    parser.add_argument('-o2', dest='output_file2', type=str,
                        default="masque_annotation_process.tsv",
                        help='Output file about the OTU annotation.')
    return parser.parse_args()


def check_file(screen):
    """Get the list of files in the directory
    """
    list_file = []
    try:
        list_file = glob.glob(screen)
        assert(len(list_file) > 0)
    except AssertionError:
        print("No file found in {0}".format(screen), file=sys.stderr)
        list_file = [""]
    return list_file


def parse_fastq(fastq_file):
    """Get length of each read
    """
    seq_len_tab = []
    try:
        with open(fastq_file, "rt") as fastq:
            for line in fastq:
                # Get the sequence
                seq_len_tab.append(len(fastq.next()))
                # Pass separator
                fastq.next()
                # Pass quality
                fastq.next()
    except IOError:
        sys.exit("Error cannot open fastq file : {0}".format(fastq_file))
    return seq_len_tab


def parse_fasta(fasta_file, tag="size=.+"):
    """Parse fasta sequence
    """
    regex_name = re.compile(r"barcodelabel=(\S+);" + tag)
    header_dict = {}
    header = ""
    sequence = ""
    seq_len_tab = []
    try:
        with open(fasta_file, "rt") as fast:
            for line in fast:
                if line.startswith(">"):
                    if len(header) > 0:
                        seq_len_tab.append(len(sequence))
                        sequence = ""
                    regex_match = regex_name.search(line)
                    if regex_match:
                        header = regex_match.group(1)
                        #print(header)
                        if header in header_dict:
                            header_dict[header] += 1
                        else:
                            header_dict[header] = 1
                else:
                    sequence += line.replace("\n", "").replace("\r", "")
            seq_len_tab.append(len(sequence))
    except IOError:
        sys.exit("Error cannot open {0}".format(fasta_file))
    return header_dict, seq_len_tab

def get_size_info(seq_len_tab):
    """Get number, mean length and median_length
    """
    return [len(seq_len_tab), sum(seq_len_tab)/len(seq_len_tab),
            sorted(seq_len_tab)[len(seq_len_tab)//2]]

def get_reads_data(sample_read, list_reads, paired, tag):
    """Count simple information on fastq sets
    """
    if paired:
        for i in xrange(len(list_reads[0])):
            # Get sample name
            name = os.path.splitext(os.path.basename(list_reads[0][i]))[0]
            name = name.replace("-R1","").replace("_R1_001","")
            name = name.replace("_alien_f_filt","")
            seq_len_tab_fwd = parse_fastq(list_reads[0][i])
            seq_len_tab_rev = parse_fastq(list_reads[1][i])
            #print(name)
            if name in sample_read:
                #print("here")
                sample_read[name].update({tag+"_fwd":get_size_info(seq_len_tab_fwd)})
                sample_read[name].update({tag+"_rev":get_size_info(seq_len_tab_rev)})
            else:
                #print("here2")
                sample_read[name] = {tag+"_fwd":get_size_info(seq_len_tab_fwd)}
                sample_read[name].update({tag+"_rev":get_size_info(seq_len_tab_rev)})
    else:
        for sample in list_reads:
            name = os.path.splitext(os.path.basename(sample))[0]
            name = name.replace("_alien_filt","")
            seq_len_tab = parse_fastq(sample)
            if name in sample_read:
                sample_read[name].update({tag:get_size_info(seq_len_tab)})
            else:
                #print(seq_len_tab)
                sample_read[name] = {tag:get_size_info(seq_len_tab)}
    return sample_read


def parse_log(log_file, regex, tag_alien):
    """Parse log data
    """
    data = []
    try:
        with open(log_file, "rt") as log:
            for line in log:
                regex_match = regex.match(line)
                if regex_match:
                    data += [int(i.replace(",","").replace("%","")) 
                            for i in regex_match.groups()]
            if tag_alien > 0 and len(data) == 0:
                data = [0] * tag_alien
            else:
                assert(len(data) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(log_file))
    except AssertionError:
        sys.exit("No data parsed from {0}".format(log_file))
    return data


def get_log(sample_read, list_file, soft, tag, paired=False):
    """
    """
    tag_alien = 0
    if soft == "alientrimmer":
        if paired :
            regex = re.compile(r".+\s+(\S+)\s+trimmed\s+\(fwd:\s+(\S+)\s+rev:\s+(\S+)\)"
                               "\s+(\S+)\s+removed\s+\(fwd:\s+(\S+)\s+rev:\s+(\S+)\)")
            tag_alien = 6
        else:
            regex = re.compile(r".+\s+(\S+)\s+trimmed\s+(\S+)\s+removed")
            tag_alien = 2
    elif soft == "flash":
        regex = re.compile(r"\S+\s+\S+\s+pairs:\s+(\S+)")
    else:
        #regex = re.compile(r"\s+(\S+)\s+\((\S+)\)\s+aligned.+1\s+time")
        regex = re.compile(r"\s+(\S+)\s+.+\s+aligned\s+(?:exactly\s+|\>)1\s+time")
    for sample in list_file:
        name = os.path.basename(sample).replace("log_{0}_".format(soft),"").replace(tag+".txt", "")
        if name in sample_read:
            sample_read[name].update({soft+tag:parse_log(sample, regex,
                                                         tag_alien)})
        else:
            print("Sample not identified : {0}".format(name), file=sys.stderr)
            #raise KeyError("Sample not identified : {0}".format(name))
    return sample_read


def update_step(sample_read, header_dict, step):
    """Update sample_read information
    """
    for sample in sample_read:
        if sample in header_dict:
            sample_read[sample].update({step:header_dict[sample]})
        else:
            sample_read[sample].update({step:0})
    return sample_read


def parse_otu_table(sample_read, otu_table_file):
    """Count number and length of OTU
    """
    try:
        with open(otu_table_file, "rt") as otu_table:
            otu_reader = csv.reader(otu_table, delimiter='\t')
            header_name = otu_reader.next()[1:]
            header = [0] * len(header_name)
            for line in otu_reader:
                # ignore ligne name
                line = line[1:]
                for i in xrange(len(header)):
                    header[i] += int(line[i])
        assert(len(header) > 0)
        for i in xrange(len(header_name)):
            sample_read[header_name[i]].update({"mapped":header[i]})
    except IOError:
        sys.exit("Error cannot open {0}".format(otu_table_file))
    except AssertionError:
        sys.exit("Error nothing read from {0}".format(otu_table_file))
    return sample_read



def write_sample_result(sample_read, output_file, paired):
    """Write otu building process data
    """
    if paired:
        header = ["sample", "Raw_reads_fwd", "Raw_mean_length_fwd",
                  "Raw_median_length_fwd", "Raw_reads_rev",
                  "Raw_mean_length_rev", "Raw_median_length_rev",
                  "Trimmed", "Trimmed_fwd", "Trimmed_rev", "Removed",
                  "Removed_fwd", "Removed_rev", "mapping_human_1_time",
                  "mapping_human_>1_time", "mapping_phiX_1_time",
                  "mapping_phiX_>1_time", "Filtered_reads_fwd",
                  "Filtered_mean_length_fwd", "Filtered_median_length_fwd",
                  "Filtered_reads_rev", "Filtered_mean_length_rev",
                  "Filtered_median_length_rev", "Combined pairs",
                  "Uncombined pairs", "Selected_dereplication",
                  "Selected_singleton", "Selected_chimera", "Selected_otu",
                  "Mapped_reads", "Mapping_percent_combined"]
    else:
        header = ["sample", "Raw_reads", "Raw_mean_length", "Raw_median_length",
                  "Trimmed", "Removed", "mapping_human_1_time",
                  "mapping_human_>1_time", "mapping_phiX_1_time",
                  "mapping_phiX_>1_time", "Filtered_reads",
                  "Filtered_mean_length", "Filtered_median_length",
                  "Selected_dereplication", "Selected_singleton",
                  "Selected_chimera", "Selected_otu",
                  "Mapped_reads", "Mapping_percent_proc"]
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter='\t')
            output_writer.writerow(header)
            for sample in sample_read:
                if paired:
                    #print(sample)
                    #print(sample_read[sample])
                    flash_data = sample_read[sample]['flash'][1:]
                    read_data = (sample_read[sample]['raw_fwd'] +
                                 sample_read[sample]['raw_rev'])
                    proc_data = (sample_read[sample]['proc_fwd'] +
                                 sample_read[sample]['proc_rev'])
                    mapping_perc = round(float(sample_read[sample]['mapped'])/
                                    float(sample_read[sample]['flash'][1])*100.0, 2)
                else:
                    flash_data = []
                    read_data = sample_read[sample]['raw']
                    proc_data = sample_read[sample]['proc']
                    mapping_perc = round(float(sample_read[sample]['mapped'])/
                                    float(sample_read[sample]['proc'][0])*100.0, 2)
                output_writer.writerow(
                    [sample] + read_data +
                    sample_read[sample]['alientrimmer']+
                    sample_read[sample]['mapping_1'] +
                    sample_read[sample]['mapping_2'] +
                    proc_data + flash_data +
                    [sample_read[sample]['dereplication'],
                    sample_read[sample]['singleton'],
                    sample_read[sample]['chimera'],
                    sample_read[sample]['otu'],
                    sample_read[sample]['mapped'],
                    mapping_perc])
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))
    #except KeyError:
    #    sys.exit("Sample {0} failed to be processed".format(sample))


def write_otu_annotation(global_data, output_file, tag):
    """Write global data
    """
    list_step = ["extended", "dereplication", "singleton", "chimera", "otu"]
    header = ["Type", "Count", "Mean_length", "Median_length"]
    rownames = ["Amplicon", "Dereplication", "Singleton_removed",
                "Chimera_removed", "OTU"]
    tag_names = [step  + "_annotation" for step in tag]
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter='\t')
            output_writer.writerow(header)
            for i in xrange(len(list_step)):
                output_writer.writerow([rownames[i]] + global_data[list_step[i]])
            for i in xrange(len(tag)):
                output_writer.writerow([tag_names[i], global_data[tag[i]]])
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    args = get_arguments()
    sample_read = {}
    ## Get raw data
    if args.raw_reads_dir:
        if args.paired_reads:
            list_r1 = check_file(args.raw_reads_dir + "*R1*.f*q")
            list_file = ([list_r1] +
                        [[r.replace("R1", "R2") for r in list_r1]])
        else:
            list_file = check_file(args.raw_reads_dir + "*.f*q")
        #print(list_file)
        sample_read = get_reads_data(sample_read, list_file, args.paired_reads,
                                     "raw")
    elif args.amplicon_dir:
        sample_read = get_reads_data(sample_read,
                                     check_file(args.raw_reads_dir + "*.f*q"),
                                     args.paired_reads, "raw")
    else:
        sys.exit("Error no raw reads or amplicon dir data")
    ## Get processed reads information
    #if not args.processed_reads_dir:
    processed_reads_dir = args.data_dir + "reads" + os.sep
    if os.path.isdir(processed_reads_dir):
        if args.paired_reads:
            list_r1 = check_file(processed_reads_dir + "*alien_f_filt.fastq")
            list_file = ([list_r1] +
                        [[r.replace("alien_f_filt", "alien_r_filt")
                          for r in list_r1]])
        else:
            list_file = check_file(processed_reads_dir + "*alien_filt.f*q")
        sample_read = get_reads_data(sample_read, list_file, args.paired_reads,
                                     "proc")
    else:
        print("Read directory is missing",file=sys.stderr)
    ## Get log data
    log_dir = args.data_dir + "log" + os.sep
    if os.path.isdir(log_dir):
        # Get log alientrimmer
        sample_read = get_log(sample_read, 
                              check_file(log_dir + "log_alientrimmer*.txt"),
                              "alientrimmer", "", args.paired_reads)
        # Get log flash
        if args.paired_reads:
            sample_read = get_log(sample_read, 
                                check_file(log_dir + "log_flash*.txt"),
                                "flash", "")
        # Get log mapping
        sample_read = get_log(sample_read,
                              check_file(log_dir + "log_mapping*_1.txt"),
                              "mapping", "_1")
        sample_read = get_log(sample_read,
                              check_file(log_dir + "log_mapping*_2.txt"),
                              "mapping", "_2")
    # Get dereplication data
    header_dict, seq_len_tab = parse_fasta(check_file(args.data_dir +
                                                      "*_extendedFrags.fasta")[0],
                                           "")
    global_data = {"extended":[
                    len(seq_len_tab), sum(seq_len_tab)/len(seq_len_tab),
                    sorted(seq_len_tab)[len(seq_len_tab)//2]]}
    # Get dereplication data
    header_dict, seq_len_tab = parse_fasta(check_file(args.data_dir +
                                                      "*_drep.fasta")[0])
    sample_read = update_step(sample_read, header_dict, "dereplication")
    global_data.update({"dereplication":[
                        len(seq_len_tab),
                        sum(seq_len_tab)/len(seq_len_tab),
                        sorted(seq_len_tab)[len(seq_len_tab)//2]]})
    ## Get singleton data
    header_dict, seq_len_tab = parse_fasta(check_file(args.data_dir +
                                                     "*_sorted.fasta")[0])
    sample_read = update_step(sample_read, header_dict, "singleton")
    global_data.update({"singleton":[
                    len(seq_len_tab),
                    sum(seq_len_tab)/len(seq_len_tab),
                    sorted(seq_len_tab)[len(seq_len_tab)//2]]})
    ## Get Chimera data
    header_dict, seq_len_tab = parse_fasta(check_file(args.data_dir +
                                                     "*_nochim.fasta")[0])
    sample_read = update_step(sample_read, header_dict, "chimera")
    global_data.update({"chimera":[
                    len(seq_len_tab),
                    sum(seq_len_tab)/len(seq_len_tab),
                    sorted(seq_len_tab)[len(seq_len_tab)//2]]})
    ## Get otu data
    header_dict, seq_len_tab = parse_fasta(check_file(args.data_dir +
                                                     "*_otu_compl.fasta")[0])
    #print(header_dict)
    sample_read = update_step(sample_read, header_dict, "otu")
    #print(sample_read)
    global_data.update({"otu":[
                    len(seq_len_tab),
                    sum(seq_len_tab)/len(seq_len_tab),
                    sorted(seq_len_tab)[len(seq_len_tab)//2]]})
    ## Get otu table
    sample_read = parse_otu_table(sample_read, check_file(args.data_dir +
                                                          "*_otu_table.tsv")[0])
    # annotation
    tag = ["silva", "greengenes", "unite", "findley", "rdp"]
    annotation_files = [check_file(args.data_dir + "*_silva_annotation_eval*.tsv")[0],
                        check_file(args.data_dir + "*_greengenes_annotation_*.tsv")[0],
                        check_file(args.data_dir + "*_unite_annotation_*.tsv")[0],
                        check_file(args.data_dir + "*_findley_annotation_*.tsv")[0],
                        check_file(args.data_dir + "*_rdp.tsv")[0]]
    tag_present = []
    for i in xrange(len(annotation_files)):
        if os.path.isfile(annotation_files[i]):
            try:
                with open(annotation_files[i], "rt") as annotation:
                    global_data.update({tag[i]:sum(1 for line in annotation) - 1})
                tag_present += [1]
            except IOError:
                sys.exit("Error cannot open {0}".format(annotation_files[i]))
        else:
            tag_present += [0]
    # write result
    write_sample_result(sample_read, args.output_file1, args.paired_reads)
    write_otu_annotation(global_data, args.output_file2,
                         [tag[i] for i in xrange(len(tag)) if tag_present[i]])

if __name__ == '__main__':
    main()
