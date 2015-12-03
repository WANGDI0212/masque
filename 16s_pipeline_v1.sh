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
# ------------------------------------------------------------------
# Author: Amine Ghozlane (amine.ghozlane@pasteur.fr)
# Title:  16S pipeline
# Description : De novo 16S pipeline assignation
# ------------------------------------------------------------------

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
        echo """$0 -i </path/to/input/directory/> -o </path/to/result/directory/>
        case V1-V3 regions of 16S rRNA - $0 -i </path/to/input/directory/> -o </path/to/result/directory/> --minoverlap --maxoverlap"""
    else
        display_parameters
    fi
    exit
}

display_parameters() {
   # Display the parameters of the analysis
   say_parameters "Sample input [-i]:"
   echo $input_dir >&2
   say_parameters "Result output [-o]:"
   echo $resultDir >&2
   say_parameters "Merge reads"
   echo "Minoverlap= $minoverlap" >&2
   echo "Maxoverlap= $maxoverlap" >&2
   say_parameters "OTU taxonomy threshold (SILVA, Greengenes)"
   echo "Identity threshold= $identity_threshold" >&2
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
# rdp="$SCRIPTPATH/databases/rdp_11_4.fa"
# Silva
#http://www.arb-silva.de/no_cache/download/archive/release_123/Exports/
silva="$SCRIPTPATH/databases/SILVA_123_SSURef_Nr99_tax_silva.fasta"
# Greengenes
#ftp://greengenes.microbio.me/greengenes_release/gg_13_5/
greengenes="$SCRIPTPATH/databases/gg_13_5.fasta"

#######################
# Assembly Parameters #
#######################
maxoverlap=200
minoverlap=50
NbProc=$(grep -c ^processor /proc/cpuinfo)
# evalueTaxAnnot="1E-5"
# maxTargetSeqs=1
identity_threshold=0.8

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
# get_taxonomy
get_taxonomy="$SCRIPTPATH/get_taxonomy/get_taxonomy.py"
# rdp classifier
rdp_classifier="$SCRIPTPATH/rdp_classifier_2.11/dist/classifier.jar"
# swarm
swarm="$SCRIPTPATH/swarm/bin/swarm"
# swarm2vsearch
swarm2vsearch="$SCRIPTPATH/swarm2vsearch/swarm2vsearch.py"
# uc2otutab
uc2otutab="$SCRIPTPATH/usearch_python_scripts/uc2otutab.py"
# usearch
usearch="$SCRIPTPATH/usearch8.1.1756_i86linux32"
#usearch -makeudb_utax 16s_ref.fa -output 16s_ref.udb -report 16s_report.txt
# vsearch
#vsearch="$SCRIPTPATH/vsearch-1.4.0-linux-x86_64/bin/vsearch"
vsearch="$SCRIPTPATH/vsearch_bin/bin/vsearch"

########
# Main #
########
# Execute getopt on the arguments passed to this program, identified by the special character $@
PARSED_OPTIONS=$(getopt -n "$0"  -o hi:o:r:n: --long "help,input_dir:,output:,NbProc:,maxoverlap:,minoverlap:,identity_threshold:"  -- "$@")

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
    -i|--input_dir)
        input_dir=$2
        shift 2;;
    -o|--output)
        resultDir=$2
        readsDir=$resultDir/reads/
        logDir=$resultDir/log/
        errorlogDir=$resultDir/error_log/
        check_dir $resultDir
        check_dir $logDir
        check_dir $errorlogDir
        check_dir $readsDir
        shift 2;;
    -n|--NbProc)
        check_integer $2
        NbProc=$2
        shift 2;;
    --maxoverlap)
        check_integer $2
        maxoverlap=$2
        shift 2;;
    --minoverlap)
        check_integer $2
        minoverlap=$2
        shift 2;;
    --identity_threshold)
        identity_threshold=$2
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

if [ -d $input_dir ]
then
    ProjectName=$(basename "$input_dir")
else
    error "No input directory !"
    exit 1
fi


# display parameters
display_parameters

# Start timer
say "Start analysis"
wall_time=$(timer)


list_product_fa=""
nb_samples=$(ls $input_dir/*R1*.fastq -1 |wc -l)
num_sample=0
for r1_file in $(ls $input_dir/*R1*.fastq)
do
    let "num_sample=$num_sample+1"
    input1=$r1_file
    input2=$(echo $r1_file|sed "s:R1:R2:g")
    check_file $input1
    check_file $input2
    # Get the sample name
    SampleName=$(basename $input1 |sed "s:_L001:@:g"|cut -f 1 -d"@")
    
    if [ "$SampleName" = "" ]
    then
        error "Error failed to parse the sample name:"
        echo "sample name=$SampleName"
        exit 1
    else
        list_product_fa+="${resultDir}/reads/${SampleName}_extendedFrags.fasta "
    fi
    if [ -f "$input1" ] && [ -f "$input2" ] && [ ! -f "${readsDir}/${SampleName}_alien_f.fastq" ]
    then
        say "$num_sample/$nb_samples - Triming reads with Alientrimmer"
        start_time=$(timer)
        java -jar $alientrimmer -if $input1 -ir $input2 -of ${readsDir}/${SampleName}_alien_f.fastq -or ${readsDir}/${SampleName}_alien_r.fastq -os ${readsDir}/${SampleName}_alien_s.fastq -c $alienseq  > ${logDir}/log_alientrimmer_${SampleName}.txt 2> ${errorlogDir}/error_log_alientrimmer_${SampleName}.txt
        check_file ${readsDir}/${SampleName}_alien_f.fastq
        check_file ${readsDir}/${SampleName}_alien_r.fastq
        check_log ${errorlogDir}/error_log_alientrimmer_${SampleName}.txt
        say "$num_sample/$nb_samples - Elapsed time with Alientrimmer : $(timer $start_time)"
    fi

    if [ -f "${readsDir}/${SampleName}_alien_f.fastq" ] && [ -f "${readsDir}/${SampleName}_alien_r.fastq" ] && [ ! -f "${readsDir}/${SampleName}.extendedFrags.fastq" ]
    then
        say "$num_sample/$nb_samples - Merging paired reads with FLASH"
        start_time=$(timer)
        $flash ${readsDir}/${SampleName}_alien_f.fastq ${readsDir}/${SampleName}_alien_r.fastq -M $maxoverlap -m $minoverlap -d $readsDir/ -o $SampleName -t $NbProc  > ${logDir}/log_flash_${SampleName}.txt
        check_file ${readsDir}/${SampleName}.extendedFrags.fastq
        say "$num_sample/$nb_samples - Elapsed time with FLASH : $(timer $start_time)"
    fi

    if [ -f "${readsDir}/${SampleName}.extendedFrags.fastq" ] && [ ! -f "${readsDir}/${SampleName}.extendedFrags_fastqc.html" ] 
    then
        say "$num_sample/$nb_samples - Quality control with Fastqc"
        start_time=$(timer)
        $fastqc ${readsDir}/${SampleName}.extendedFrags.fastq --nogroup -q -t $NbProc 2> ${errorlogDir}/error_log_fastqc_${SampleName}.txt
        check_file ${readsDir}/${SampleName}.extendedFrags_fastqc.html
        check_log ${errorlogDir}/error_log_fastqc_${SampleName}.txt
        say "$num_sample/$nb_samples - Elapsed time with Fastqc: $(timer $start_time)"
    fi

    if [ -f "${readsDir}/${SampleName}.extendedFrags.fastq" ] && [ ! -f ${readsDir}/${SampleName}_extendedFrags.fasta ]
    then
        say "$num_sample/$nb_samples - Convert fastq to fasta with fastq2fasta"
        start_time=$(timer)
        $fastq2fasta -i ${readsDir}/${SampleName}.extendedFrags.fastq -o ${readsDir}/${SampleName}_extendedFrags.fasta -s ${SampleName}  2> ${errorlogDir}/error_log_fastq2fasta_${SampleName}.txt
        check_file ${readsDir}/${SampleName}_extendedFrags.fasta
        check_log ${errorlogDir}/error_log_fastq2fasta_${SampleName}.txt
        say "$num_sample/$nb_samples - Elapsed time with fastq2fasta : $(timer $start_time)"
    fi
done




# Combine all files
if [ ! -f ${resultDir}/${ProjectName}_extendedFrags.fasta ]
then
    say "Combine fasta files"
    start_time=$(timer)
    cat $list_product_fa > ${resultDir}/${ProjectName}_extendedFrags.fasta
    say "Elapsed time to combine fasta files : $(timer $start_time)"
fi

if [ -f "${resultDir}/${ProjectName}_extendedFrags.fasta" ] && [ ! -f "${resultDir}/${ProjectName}_drep.fasta" ]
then
     say "Dereplication"
     start_time=$(timer)
     #$usearch -derep_fulllength ${resultDir}/${ProjectName}.extendedFrags.fasta -fastaout ${resultDir}/${ProjectName}_drep.fasta -sizeout 
     # -minseqlength 64
     $vsearch --derep_fulllength ${resultDir}/${ProjectName}_extendedFrags.fasta -output ${resultDir}/${ProjectName}_drep.fasta -sizeout
     check_file ${resultDir}/${ProjectName}_drep.fasta
     say "Elapsed time to dereplicate: $(timer $start_time)"
fi

if [ -f "${resultDir}/${ProjectName}_drep.fasta" ] && [ ! -f "${resultDir}/${ProjectName}_sorted.fasta" ]
then
     say "Abundance sort and discard singletons"
     start_time=$(timer)
     #$usearch -sortbysize ${resultDir}/${ProjectName}_drep.fasta -fastaout ${resultDir}/${ProjectName}_sorted.fasta -minsize 4
     $vsearch -sortbysize ${resultDir}/${ProjectName}_drep.fasta -output ${resultDir}/${ProjectName}_sorted.fasta -minsize 4
 > ${logDir}/log_search_sort_${ProjectName}.txt 2>&1
     check_file ${resultDir}/${ProjectName}_sorted.fasta
     say "Elapsed time to sort: $(timer $start_time)"
fi


if [ -f "${resultDir}/${ProjectName}_sorted.fasta" ] && [ ! -f "${resultDir}/${ProjectName}_otu.fasta" ]
then
     say "OTU clustering"
     start_time=$(timer)
     #$usearch -cluster_otus ${resultDir}/${ProjectName}_sorted.fasta -otus ${resultDir}/${ProjectName}_otu.fasta -uparseout ${resultDir}/${ProjectName}_uparse.txt -relabel OTU_ -sizein #-sizeout 
     #$vsearch --cluster_size ${resultDir}/${ProjectName}_sorted.fasta --id 0.97 --centroids ${resultDir}/${ProjectName}_otu.fasta --relabel OTU_ --sizein #--sizeout
#      echo "$swarm -t $NbProc -f -z -w ${resultDir}/${ProjectName}_representant_swarm.fasta -o ${resultDir}/${ProjectName}_swarm_clustering.txt -s ${resultDir}/${ProjectName}_swarm_stats.txt ${resultDir}/${ProjectName}_sorted.fasta"
     $swarm -t $NbProc -f -z -w ${resultDir}/${ProjectName}_representant_swarm.fasta -o ${resultDir}/${ProjectName}_swarm_clustering.txt -s ${resultDir}/${ProjectName}_swarm_stats.txt ${resultDir}/${ProjectName}_sorted.fasta
#      echo "python $swarm2vsearch -i ${resultDir}/${ProjectName}_representant_swarm.fasta  -c ${resultDir}/${ProjectName}_swarm_clustering.txt -o ${resultDir}/${ProjectName}_otu.fasta -oc ${resultDir}/${resultDir}/${ProjectName}_otu_swarm_clustering.txt"
     python $swarm2vsearch -i ${resultDir}/${ProjectName}_representant_swarm.fasta  -c ${resultDir}/${ProjectName}_swarm_clustering.txt -o ${resultDir}/${ProjectName}_otu.fasta -oc ${resultDir}/${ProjectName}_otu_swarm_clustering.txt
     check_file ${resultDir}/${ProjectName}_otu.fasta
     say "Elapsed time to OTU clustering: $(timer $start_time)"
fi

if [ -f "${resultDir}/${ProjectName}_otu.fasta" ] &&  [ ! -f "${resultDir}/${ProjectName}_otu_nochim.fasta" ]
then
     say "Chimera filtering using reference database"
     start_time=$(timer)
     #$usearch -uchime_ref ${resultDir}/${ProjectName}_otu.fasta -db $gold -strand plus -nonchimeras ${resultDir}/${ProjectName}_otu_nochim.fasta
     $vsearch --uchime_ref ${resultDir}/${ProjectName}_otu.fasta --db $gold --strand plus --nonchimeras ${resultDir}/${ProjectName}_otu_nochim.fasta
     check_file ${resultDir}/${ProjectName}_otu_nochim.fasta 
     say "Elapsed time to filter chimera: $(timer $start_time)"
fi

if [ -f "${resultDir}/${ProjectName}_otu_nochim.fasta" ] &&  [ ! -f "${resultDir}/${ProjectName}_map.txt" ]
then
    say "Map reads back to OTUs with vsearch"
    start_time=$(timer)
    #$usearch -usearch_global ${resultDir}/${ProjectName}_extendedFrags.fasta -db ${resultDir}/${ProjectName}_otu_nochim.fasta -strand plus -id 0.97 -uc ${resultDir}/${ProjectName}_map.txt
    $vsearch -usearch_global ${resultDir}/${ProjectName}_extendedFrags.fasta -db ${resultDir}/${ProjectName}_otu_nochim.fasta --strand plus --id 0.97 -uc ${resultDir}/${ProjectName}_map.txt
    check_file ${resultDir}/${ProjectName}_map.txt
    say "Elapsed time to map reads: $(timer $start_time)"
fi

if [ -f ${resultDir}/${ProjectName}_map.txt ] && [ ! -f ${resultDir}/${ProjectName}_otu_table.txt ]
then
    say "Build OTUs table"
    start_time=$(timer)
    python $uc2otutab ${resultDir}/${ProjectName}_map.txt > ${resultDir}/${ProjectName}_otu_table.txt
    say "Elapsed time to build OTUs table: $(timer $start_time)"
fi

if [ -f "${resultDir}/${ProjectName}_otu_nochim.fasta" ]
then
    if [ ! -f "${resultDir}/${ProjectName}_vs_rdp.txt" ]
    then
        say "Assign taxonomy with rdp_classifier"
        start_time=$(timer)
        #$usearch -utax ${resultDir}/${ProjectName}_otu_nochim.fasta -db $rdp -strand both -taxconfs rdp_16s_short.tc -utaxout ${resultDir}/${ProjectName}_otu_tax_rdp.txt -utax_cutoff 0.8
        #$vsearch  --usearch_global ${resultDir}/${ProjectName}_otu_nochim.fasta --db $rdp --id 0.9 --blast6out ${resultDir}/${ProjectName}_vs_rdp.txt
        java -Xmx1g -jar $rdp_classifier classify  -q ${resultDir}/${ProjectName}_otu_nochim.fasta -o  ${resultDir}/${ProjectName}_vs_rdp.txt
        check_file ${resultDir}/${ProjectName}_vs_rdp.txt
        #python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_silva.txt -d $rdp -dtype rdp -o ${resultDir}/${ProjectName}_vs_rdp_annotation.txt
        say "Elapsed time with rdp_classifier: $(timer $start_time)"
    fi
    if [ ! -f "${resultDir}/${ProjectName}_vs_silva_${identity_threshold}.txt" ]
    then
        say "Assign taxonomy against silva with vsearch"
        start_time=$(timer)
        #$usearch -utax ${resultDir}/${ProjectName}_otu_nochim.fasta -db $silva -strand both -taxconfs silva_16s_short.tc -utaxout ${resultDir}/${ProjectName}_otu_tax_silva.txt -utax_cutoff 0.8
        $vsearch --usearch_global ${resultDir}/${ProjectName}_otu_nochim.fasta --db $silva --id $identity_threshold --blast6out ${resultDir}/${ProjectName}_vs_silva_${identity_threshold}.txt
        #$blastn -query $input_gene -db $silva -outfmt 6 -evalue $evalueTaxAnnot -num_threads $NbProc -out ${resultDir}/${ProjectName}_vs_silva.txt -max_target_seqs $maxTargetSeqs -task megablast -outfmt "6 qseqid sseqid  pident qcovs evalue"
        check_file ${resultDir}/${ProjectName}_vs_silva_${identity_threshold}.txt
        say "Elapsed time with vsearch: $(timer $start_time)"
    fi
    if [ -f "${resultDir}/${ProjectName}_vs_silva_${identity_threshold}.txt" ] && [ ! -f "${resultDir}/${ProjectName}_vs_silva_annotation_${identity_threshold}.txt" ]
    then
        say "Extract silva annotation with get_taxonomy"
        start_time=$(timer)
        python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_silva_${identity_threshold}.txt -d $silva -o ${resultDir}/${ProjectName}_vs_silva_annotation_${identity_threshold}.txt
        check_file ${resultDir}/${ProjectName}_vs_silva_annotation_${identity_threshold}.txt
        say "Elapsed time with get_taxonomy: $(timer $start_time)"
    fi
#     if [ ! -f "${resultDir}/${ProjectName}_vs_greengenes_${identity_threshold}.txt" ]
#     then
#         say "Assign taxonomy against greengenes with vsearch"
#         start_time=$(timer)
#         $vsearch --usearch_global ${resultDir}/${ProjectName}_otu_nochim.fasta --db $greengenes --id $identity_threshold --blast6out ${resultDir}/${ProjectName}_vs_greengenes_${identity_threshold}.txt
#         check_file ${resultDir}/${ProjectName}_vs_greengenes_${identity_threshold}.txt
#         say "Elapsed time with vsearch: $(timer $start_time)"
#     fi
#         if [ -f "${resultDir}/${ProjectName}_vs_greengenes_${identity_threshold}.txt" ] && [ ! -f "${resultDir}/${ProjectName}_vs_greengenes_annotation_${identity_threshold}.txt" ]
#     then
#         say "Extract greengenes annotation with get_taxonomy"
#         start_time=$(timer)
#         python $get_taxonomy -i ${resultDir}/${ProjectName}_vs_greengenes_${identity_threshold}.txt -d $greengenes -o ${resultDir}/${ProjectName}_vs_greengenes_annotation_${identity_threshold}.txt -dtype greengenes
#         check_file ${resultDir}/${ProjectName}_vs_greengenes_annotation_${identity_threshold}.txt
#         say "Elapsed time with get_taxonomy: $(timer $start_time)"
#     fi
fi


say "16S analysis is done. Elapsed time: $(timer $wall_time)"