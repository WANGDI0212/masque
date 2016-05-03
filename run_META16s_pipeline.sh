#!/bin/bash
PBS_SCRIPT=$HOME/submission_16s.sh

#Check arguments
if [ $# -ne 4  ]
then
    echo "$0 <reads_dir> <output_dir> <name> <nb_cpu>"
    exit
fi

mkdir -p $2

echo """#!/bin/bash
#$ -S /bin/bash
#$ -M amine.ghozlane@pasteur.fr
#$ -m bea
##$ -q test
#$ -q hubbioit
#$ -pe thread $4
#$ -l mem_total=50G
#$ -N "META_$3"
source /local/gensoft2/adm/etc/profile.d/modules.sh
module purge
module add Python/2.7.8 FastTree/2.1.8 FLASH/1.2.11 fasta mafft/7.149 bowtie2/2.2.3 blast+/2.2.31 vsearch/1.4.1 AlienTrimmer/0.4.0 fastqc/0.11.2 rdp_classifier/2.11 swarm/2.1.1

/bin/bash $HOME/META16s_pipeline/META16s_pipeline.sh -i $1/ -o $2/ -n $3 -t $4 --minoverlap 10 --maxoverlap 550 -b &> $2/${3}_stat_process.txt
""">$PBS_SCRIPT
PBSID=`qsub $PBS_SCRIPT`
echo "! Soumission PBS :> JOBID = $PBSID"
