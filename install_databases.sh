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
# Title:  16S-18S-23S-28S-ITS pipeline
# Description : De novo 16S-18S-23S-28S-ITS pipeline assignation
# ------------------------------------------------------------------
function say {
    # echo green
    echo -e "\033[1;32m* $1\033[0m" >&2
}

function check_soft {
 if type $1 > /dev/null  2>&1;
 then
    echo $1
 else
    echo $2
 fi
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

#############
# Variables #
#############
SCRIPTPATH=$(dirname "${BASH_SOURCE[0]}")
databases_dir="$SCRIPTPATH/databases/"

#############
# Databases #
#############
blast_databases="$SCRIPTPATH/blast_indexing.txt"
bowtie2_databases="$SCRIPTPATH/bowtie2_indexing.txt"
md5_check="$SCRIPTPATH/md5_check.txt"
NbProc=$(grep -c ^processor /proc/cpuinfo)
url_list="$SCRIPTPATH/url_list.txt"

############
# Programs #
############
bowtie2_build=$(check_soft "bowtie2-build" "$SCRIPTPATH/bowtie2-2.2.9/bowtie2-build")
makeblastdb=$(check_soft "makeblastdb" "$SCRIPTPATH/ncbi-blast-2.4.0+/bin/makeblastdb")
makembindex=$(check_soft "makembindex" "$SCRIPTPATH/ncbi-blast-2.4.0+/bin/makembindex")

########
# Main #
########
# Start timer
say "Start databases preparation for MASQUE"
wall_time=$(timer)

# Download databases
say "Download databases: 1.6GB will be download."
start_time=$(timer)
wget -i $url_list -P $databases_dir
say "Elapsed time to download : $(timer $start_time)"

# Check md5sum
say "Check md5sum"
start_time=$(timer)
cd $SCRIPTPATH/databases/ && md5sum -c $md5_check
cd $SCRIPTPATH/
say "Elapsed time to check md5 : $(timer $start_time)"

# Decompress database
say "Decompress databases"
start_time=$(timer)
gunzip $databases_dir/*.gz
unzip -f -j '$databases_dir/*.zip' -d $databases_dir 
say "Elapsed time to decompress databases : $(timer $start_time)"

# Homo sapiens one file
say "Prepare Homo sapiens database"
start_time=$(timer)
cat $databases_dir/hs_ref_GRCh38.p7_*.fa > $databases_dir/homo_sapiens.fna
say "Elapsed time for Homo sapiens databases : $(timer $start_time)"

# Blast indexing
say "Indexing databases for blast"
start_time=$(timer)
while read fasta_file
do
    $makeblastdb -in $databases_dir/$fasta_file -dbtype nucl
    $makembindex -input $databases_dir/$fasta_file -iformat blastdb
done < $blast_databases
say "Elapsed time to index for blast : $(timer $start_time)"

# Bowtie2 indexing
say "Indexing databases for bowtie2"
start_time=$(timer)
while read fasta_file
do
    version=$($bowtie2_build --version |grep "bowtie2-build version"|cut -f 3 -d ' ')
    if [ "$version"  == "2.2.9" ]
    then
        $bowtie2_build --threads $NbProc $databases_dir/$fasta_file $databases_dir/$fasta_file
    else
        $bowtie2_build $NbProc $databases_dir/$fasta_file $databases_dir/$fasta_file
    fi
done < $bowtie2_databases
say "Elapsed time to index for bowtie2 : $(timer $start_time)"

# Cleanup
rm -f $databases_dir/*.zip $databases_dir/sh_general_release_dynamic_31.01.2016_dev.fasta $databases_dir/README.txt  $databases_dir/._README.txt $databases_dir/ITSdb.findley.taxonomy $databases_dir/._ITSdb.findley.fasta $databases_dir/._ITSdb.findley.taxonomy $databases_dir/._ITSdb_v1 $databases_dir/hs_ref_GRCh38.p7_*.fa 

say "Databases preparation is done. Elapsed time: $(timer $wall_time)"