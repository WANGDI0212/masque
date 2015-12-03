#!/bin/bash
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


function say_parameters {
    # echo blue
    echo -e "\e[34m# $1\e[0m" >&2
}

function say {
    # echo green
    echo -e "\033[1;32m* $1\033[0m" >&2
}

function error {
    # echo red
    echo -e  "\e[31m* $1\e[0m" >&2
}

function check_log {
    # check if log file is not empty
    if [ -s $1 ]
    then
        error "$1 is not empty !"
        exit 1
    fi
}

function check_integer {
    if [[ $1 != [0-9]* ]]
    then
        error "\"$1\" is not an integer value"
        exit 1
    fi
}

function check_file {
    # Check if result is well produced
    if [ ! -f $1 ] && [ ! -s $1 ]
    then
        error "File \"$1\" does not exist or is empty !"
        exit 1
    fi
}

function check_dir {
    # Check if directory doesnt exist
    if [ ! -d $1 ]
    then
        mkdir $1
        if [ ! -d $1 ]
        then 
            error "The program cannot create the directory \"$1\""
            exit 1
        fi
    fi
}

display_help() {
    if [ "$1" -eq "0" ]
    then
        echo """$0 -1 <read_R1.fastq> -2 <read_R2.fastq> -s <sample_name> -o </path/to/result/directory/>"""
    else
        display_parameters
    fi
    exit
}

display_parameters() {
   # Display the parameters of the analysis
   say_parameters "Sample fastq [-1] [-2] :"
   echo $input1 >&2
   echo $input2 >&2
}

function timer()
{
   if [[ $# -eq 0 ]]; then
         echo $(date '+%s')
   else
      local  stime=$1
      etime=$(date '+%s')
      if [[ -z "$stime" ]]; then stime=$etime; fi
      dt=$((etime - stime))
      ds=$((dt % 60))
      dm=$(((dt / 60) % 60))
      dh=$((dt / 3600))
      printf '%d:%02d:%02d' $dh $dm $ds
  fi
}

SCRIPTPATH=$(dirname "${BASH_SOURCE[0]}")

#############
# Databases #
#############
# ChimeraSlayer reference database
#http://drive5.com/uchime/uchime_download.html
gold="$SCRIPTPATH/databases/gold.fa"
# Alien sequences
alienseq="$SCRIPTPATH/databases/alienTrimmerPF8contaminants.fasta"
# RDP
#http://rdp.cme.msu.edu/misc/resources.jsp
rdp="$SCRIPTPATH/databases/rdp_11_4.fa"
# Silva
#http://www.arb-silva.de/no_cache/download/archive/release_123/Exports/
silva="$SCRIPTPATH/databases/SILVA_123_SSURef_Nr99_tax_silva.fasta"

#######################
# Assembly Parameters #
#######################
maxoverlap=200
minoverlap=50
NbProc=$(grep -c ^processor /proc/cpuinfo)

############
# Programs #
############
# AlienTrimmer
alientrimmer="$SCRIPTPATH/AlienTrimmer_0.4.0/src/AlienTrimmer.jar"
# Fastq2fasta
fastq2fasta="$SCRIPTPATH/fastq2fasta/fastq2fasta.py"
# Fastqc
fastqc="$SCRIPTPATH/FastQC/fastqc"
# FLASH
flash="$SCRIPTPATH/FLASH-1.2.11/flash" #$(which flash)
# usearch
usearch="$SCRIPTPATH/usearch8.1.1756_i86linux32"
#usearch -makeudb_utax 16s_ref.fa -output 16s_ref.udb -report 16s_report.txt
# vsearch
vsearch="$SCRIPTPATH/vsearch-1.1.3-linux-x86_64"

########
# Main #
########
# Execute getopt on the arguments passed to this program, identified by the special character $@
PARSED_OPTIONS=$(getopt -n "$0"  -o hs:1:2:o:r:n: --long "help,sample_name:,input1:,input2:,output:,NbProc:,SamplePrefix:"  -- "$@")

#Check arguments
if [ $# -eq 0 ]
then
    display_help 0
fi

#Bad arguments, something has gone wrong with the getopt command.
if [ $? -ne 0 ];
then
    display_help 1
    exit 1
fi

# A little magic, necessary when using getopt.
eval set -- "$PARSED_OPTIONS"

# Get Cmd line arguments depending on options
while true;
do
  case "$1" in
    -h|--help)
        display_help 0
        shift;;
    -s|--sample_name)
        sample_name=$2
        shift 2;;
    -1|--input1) 
        check_file $2
        input1=$2
        shift 2;;
    -2|--input2)
        check_file $2
        input2=$2
        shift 2;;
    -r|--SamplePrefix)
        SamplePrefix=$2
        shift 2;;
    -o|--output)
        resultDir=$2
        logDir=$resultDir/log/
        errorlogDir=$resultDir/error_log/
        check_dir $resultDir
        check_dir $logDir
        check_dir $errorlogDir
        shift 2;;
    -n|--NbProc)
        check_integer $2
        NbProc=$2
        shift 2;;
    --)
      shift
      break;;
  esac
done


if [ "$resultDir" = "" ]
then
    error "Please indicate the output directory."
    exit 1
fi

# Check sample
if [ -f "$input1" ] && [ -f "$input2" ]
then
    if [ "$sample_name" = ""  ]
    then
      error "Please provide a sample name [-s]"
      exit 1
    fi
    filename=$(basename "$input1")
    extension=".${filename##*.}"
    if [ "$extension" != ".fastq" ] && [ "$extension" != ".fq" ]
    then
        error "The input file should be a fastq file."
        display_help
    fi
    SamplePath1=$(dirname $input1)
    SamplePath2=$(dirname $input2)
    if [ "$SamplePath1" != "$SamplePath2" ]
    then
        error "The sample must be in the same directory (for the moment)."
        exit 1
    else
        SamplePath="$SamplePath1"
    fi
    SampleName=$sample_name
else
    error "No input file !"
    exit 1
fi

# display parameters
display_parameters

# Start timer
say "Start analysis"
wall_time=$(timer)

if [ -f "$input1" ] && [ -f "$input2" ] && [ ! -f "${resultDir}/${SampleName}_alien_f.fastq" ]
then
    say "Triming reads with Alientrimmer"
    start_time=$(timer)
    java -jar $alientrimmer -if $input1 -ir $input2 -of ${resultDir}/${SampleName}_alien_f.fastq -or ${resultDir}/${SampleName}_alien_r.fastq -os ${resultDir}/${SampleName}_alien_s.fastq -c $alienseq  > ${logDir}/log_alientrimmer_${SampleName}.txt 2> ${errorlogDir}/error_log_alientrimmer_${SampleName}.txt
    check_file ${resultDir}/${SampleName}_alien_f.fastq
    check_file ${resultDir}/${SampleName}_alien_r.fastq
    check_log ${errorlogDir}/error_log_alientrimmer_${SampleName}.txt
    say "Elapsed time to trim with alientrimmer : $(timer $start_time)"
fi

if [ -f "${resultDir}/${SampleName}_alien_f.fastq" ] && [ -f "${resultDir}/${SampleName}_alien_r.fastq" ] && [ ! -f "${resultDir}/${SampleName}.extendedFrags.fastq" ]
then
    say "Merging paired reads with FLASH"
    start_time=$(timer)
    $flash ${resultDir}/${SampleName}_alien_f.fastq ${resultDir}/${SampleName}_alien_r.fastq -M $maxoverlap -m $minoverlap -d $resultDir/ -o $sample_name -t $NbProc  > ${logDir}/log_flash_${SampleName}.txt
    check_file ${resultDir}/${SampleName}.extendedFrags.fastq
    say "Elapsed time to merging paired reads : $(timer $start_time)"
fi

if [ -f "${resultDir}/${SampleName}.extendedFrags.fastq" ] && [ ! -f "${resultDir}/${SampleName}.extendedFrags_fastqc.html" ] 
then
    say "Quality control with Fastqc"
    start_time=$(timer)
    $fastqc ${resultDir}/${SampleName}.extendedFrags.fastq --nogroup -q -t $NbProc 2> ${errorlogDir}/error_log_fastqc_${SampleName}.txt
    check_file ${resultDir}/${SampleName}.extendedFrags_fastqc.html
    check_log ${errorlogDir}/error_log_fastqc_${SampleName}.txt
    say "Elapsed time to quality control: $(timer $start_time)"
fi


# FLASH THEN ALIENTRIMMER
# if [ -f "$input1" ] && [ -f "$input2" ] && [ ! -f "${resultDir}/${SampleName}.extendedFrags.fastq" ]
# then
#   say "Merging paired reads"
#   start_time=$(timer)
#   $flash $input1 $input2 -M $maxoverlap -m $minoverlap -d $resultDir/ -o $sample_name -t $NbProc  > ${logDir}/log_flash_${SampleName}.txt
#   check_file ${resultDir}/${SampleName}.extendedFrags.fastq
#   say "Elapsed time to merging paired reads : $(timer $start_time)"
# fi
# 
# if [ -f "${resultDir}/${SampleName}.extendedFrags.fastq" ] && [ ! -f ${resultDir}/${SampleName}_alien.fastq ]
# then
#     say "Triming reads with alientrimmer"
#     start_time=$(timer)
#     java -jar $alientrimmer -i ${resultDir}/${SampleName}.extendedFrags.fastq -o ${resultDir}/${SampleName}_alien.fastq -c $alienseq  > ${logDir}/log_alientrimmer_${SampleName}.txt 2> ${errorlogDir}/error_log_alientrimmer_${SampleName}.txt
#     check_file ${resultDir}/${SampleName}_alien.fastq
#     check_log ${errorlogDir}/error_log_alientrimmer_${SampleName}.txt
#     say "Elapsed time to trim with alientrimmer : $(timer $start_time)"
# fi

#if [ -f "${resultDir}/${SampleName}_alien.fastq" ]



if [ -f "${resultDir}/${SampleName}.extendedFrags.fastq" ] && [ ! -f ${resultDir}/${SampleName}_extendedFrags.fasta ]
then
    say "Convert fastq to fasta..."
    start_time=$(timer)
    $fastq2fasta -i ${resultDir}/${SampleName}.extendedFrags.fastq -o ${resultDir}/${SampleName}_extendedFrags.fasta  2> ${errorlogDir}/error_log_fastq2fasta_${SampleName}.txt
    check_file ${resultDir}/${SampleName}_extendedFrags.fasta
    check_log ${errorlogDir}/error_log_fastq2fasta_${SampleName}.txt
    say "Elapsed time to convert fastq to fasta : $(timer $start_time)"
fi

if [ -f "${resultDir}/${SampleName}_extendedFrags.fasta" ] && [ ! -f "${resultDir}/${SampleName}_drep.fasta" ]
then
     say "Dereplication"
     start_time=$(timer)
     #$usearch -derep_fulllength ${resultDir}/${SampleName}.extendedFrags.fasta -fastaout ${resultDir}/${SampleName}_drep.fasta -sizeout 
     # -minseqlength 64
     $vsearch --derep_fulllength ${resultDir}/${SampleName}_extendedFrags.fasta -output ${resultDir}/${SampleName}_drep.fasta -sizeout
     check_file ${resultDir}/${SampleName}_drep.fasta
     say "Elapsed time to dereplicate: $(timer $start_time)"
fi

if [ -f "${resultDir}/${SampleName}_drep.fasta" ] && [ ! -f "${resultDir}/${SampleName}_sorted.fasta" ]
then
     say "Abundance sort and discard singletons"
     start_time=$(timer)
     #$usearch -sortbysize ${resultDir}/${SampleName}_drep.fasta -fastaout ${resultDir}/${SampleName}_sorted.fasta -minsize 4
     $vsearch -sortbysize ${resultDir}/${SampleName}_drep.fasta -output ${resultDir}/${SampleName}_sorted.fasta -minsize 4
 > ${logDir}/log_search_sort_${SampleName}.txt 2>&1
     check_file ${resultDir}/${SampleName}_sorted.fasta
     say "Elapsed time to sort: $(timer $start_time)"
fi


if [ -f "${resultDir}/${SampleName}_sorted.fasta" ] && [ ! -f "${resultDir}/${SampleName}_otu.fasta" ]
then
     say "OTU clustering"
     start_time=$(timer)
     $usearch -cluster_otus ${resultDir}/${SampleName}_sorted.fasta -otus ${resultDir}/${SampleName}_otu.fasta -uparseout ${resultDir}/${SampleName}_uparse.txt -relabel OTU_ -sizein -sizeout 
     #$vsearch --cluster_size ${resultDir}/${SampleName}_sorted.fasta --id 0.97 --centroids ${resultDir}/${SampleName}_otu.fasta --relabel OTU_ --sizein --sizeout
     check_file ${resultDir}/${SampleName}_otu.fasta 
     say "Elapsed time to OTU clustering: $(timer $start_time)"
fi

if [ -f "${resultDir}/${SampleName}_otu.fasta" ] &&  [ ! -f "${resultDir}/${SampleName}_otu_nochim.fasta" ]
then
     say "Chimera filtering using reference database"
     start_time=$(timer)
     #$usearch -uchime_ref ${resultDir}/${SampleName}_otu.fasta -db $gold -strand plus -nonchimeras ${resultDir}/${SampleName}_otu_nochim.fasta
     $vsearch --uchime_ref ${resultDir}/${SampleName}_otu.fasta --db $gold --strand plus --nonchimeras ${resultDir}/${SampleName}_otu_nochim.fasta
     check_file ${resultDir}/${SampleName}_otu_nochim.fasta 
     say "Elapsed time to filter chimera: $(timer $start_time)"
fi

if [ -f "${resultDir}/${SampleName}_otu_nochim.fasta" ] &&  [ ! -f "${resultDir}/${SampleName}_map.txt" ]
then
    say "Map reads back to OTUs"
    start_time=$(timer)
    #$usearch -usearch_global ${resultDir}/${SampleName}_extendedFrags.fasta -db ${resultDir}/${SampleName}_otu_nochim.fasta -strand plus -id 0.97 -uc ${resultDir}/${SampleName}_map.txt
    $vsearch -usearch_global ${resultDir}/${SampleName}_extendedFrags.fasta -db ${resultDir}/${SampleName}_otu_nochim.fasta --strand plus --id 0.97 -uc ${resultDir}/${SampleName}_map.txt
    check_file ${resultDir}/${SampleName}_map.txt
    say "Elapsed time to map reads: $(timer $start_time)"
fi

if [ -f "${resultDir}/${SampleName}_otu_nochim.fasta" ]
then
    if [ ! -f "${resultDir}/${SampleName}_vs_rdp.txt" ]
    then
        say "Assign taxonomy against rdp"
        start_time=$(timer)
        #$usearch -utax ${resultDir}/${SampleName}_otu_nochim.fasta -db $rdp -strand both -taxconfs rdp_16s_short.tc -utaxout ${resultDir}/${SampleName}_otu_tax.txt -utax_cutoff 0.8
        #--usearch_global
        $vsearch  --cluster_fast ${resultDir}/${SampleName}_otu_nochim.fasta --db $rdp --id 0.9 --blast6out ${resultDir}/${SampleName}_vs_rdp.txt 
        say "Elapsed time to assign taxonomy against rdp: $(timer $start_time)"
    fi
    if [ ! -f "${resultDir}/${SampleName}_vs_silva.txt" ]
    then
        say "Assign taxonomy against silva"
        start_time=$(timer)
        # --usearch_global
        $vsearch --cluster_fast ${resultDir}/${SampleName}_otu_nochim.fasta --db $silva --id 0.9 --blast6out ${resultDir}/${SampleName}_vs_silva.txt
        say "Elapsed time to assign taxonomy against silva: $(timer $start_time)"
    fi
fi

say "16S analysis is done. Elapsed time: $(timer $wall_time)"