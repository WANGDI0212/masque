#!/bin/bash
PBS_SCRIPT=$HOME/masque_submission.sh

#Check arguments
if [ $# -ne 6  ]
then
    echo "Usage: $0 <reads_dir> <output_dir> <project-name> <nb_cpu> <email> <queue>"
    exit
fi

mkdir -p $2

readdir=$(readlink -e "$1")
outdir=$(readlink -e "$2")

SCRIPTPATH=$(dirname "${BASH_SOURCE[0]}")

echo """#!/bin/bash
#$ -S /bin/bash
#$ -M $5
#$ -m bea
#$ -q $6
#$ -pe thread $4
#$ -l mem_total=50G
#$ -N "masque_$3"
source /local/gensoft2/adm/etc/profile.d/modules.sh
module purge
export PATH=/pasteur/projets/Matrix/metagenomics/python-lib/bin:$PATH
export PYTHONPATH=/pasteur/projets/Matrix/metagenomics/python-lib/lib/python2.7/site-packages:$PYTHONPATH
module add Python/2.7.8 FastTree/2.1.8 FLASH/1.2.11 fasta mafft/7.149 bowtie2/2.2.9 blast+/2.2.40  AlienTrimmer/0.4.0 fastqc/0.11.5 rdp_classifier/2.12 BMGE/1.12 openmpi/2.0.1 IQ-TREE/1.5.1

/bin/bash $SCRIPTPATH/masque.sh -i $readdir/ -o $outdir/ -n $3 -t $4 -b &> $outdir/${3}_stat_process.txt
""">$PBS_SCRIPT
PBSID=`qsub $PBS_SCRIPT`
echo "Submission PBS :> JOBID = $PBSID"
