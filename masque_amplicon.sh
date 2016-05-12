#!/bin/bash
PBS_SCRIPT=$HOME/masque_submission.sh

#Check arguments
if [ $# -ne 6  ]
then
    echo "$0 <amplicon_file> <output_dir> <project-name> <nb_cpu> <email> <queue>"
    exit
fi

mkdir -p $2

SCRIPTPATH=$(dirname "${BASH_SOURCE[0]}")

echo """#!/bin/bash
#$ -S /bin/bash
#$ -M $5
#$ -m bea
#$ -q $6
#$ -N "masque_$3"
#$ -pe thread $4
#$ -l mem_total=50G
source /local/gensoft2/adm/etc/profile.d/modules.sh
module purge
module add Python/2.7.8 FastTree/2.1.8 FLASH/1.2.11 fasta mafft/7.149 bowtie2/2.2.6 blast+/2.2.31 vsearch/1.4.1 AlienTrimmer/0.4.0 fastqc/0.11.5 rdp_classifier/2.11

/bin/bash $SCRIPTPATH/masque.sh -a $1 -o $2/ -t $4 -n $3 -b  &> $2/${3}_stat_process.txt
""">$PBS_SCRIPT
PBSID=`qsub $PBS_SCRIPT`
echo "! Soumission PBS :> JOBID = $PBSID"
