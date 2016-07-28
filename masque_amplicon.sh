#!/bin/bash
PBS_SCRIPT=$HOME/masque_submission.sh

#Check arguments
if [ $# -ne 6  ]
then
    echo "Usage: $0 <amplicon_file> <output_dir> <project-name> <nb_cpu> <email> <queue>"
    exit
fi

mkdir -p $2

amplicon=$(readlink -e "$1")
outdir=$(readlink -e "$2")

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
export PATH=/pasteur/projets/Matrix/metagenomics/python-lib/bin:$PATH
export PYTHONPATH=/pasteur/projets/Matrix/metagenomics/python-lib/lib/python2.7/site-packages:$PYTHONPATH
module add Python/2.7.8 FastTree/2.1.8 FLASH/1.2.11 fasta mafft/7.149 bowtie2/2.2.9 blast+/2.2.40 AlienTrimmer/0.4.0 fastqc/0.11.5 rdp_classifier/2.12

/bin/bash $SCRIPTPATH/masque.sh -a $amplicon -o $outdir/ -t $4 -n $3 -b  &> $outdir/${3}_stat_process.txt || exit 1
exit 0
""">$PBS_SCRIPT
PBSID=`qsub $PBS_SCRIPT`
echo "! Soumission PBS :> JOBID = $PBSID"
