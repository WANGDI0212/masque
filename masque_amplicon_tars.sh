#!/bin/bash

SLURM_SCRIPT=$HOME/masque_submission.sh

#Check arguments
if [ $# -lt 8 ]
then
    echo "Usage: $0 <amplicon_file> <output_dir> <contaminants> <project-name> <nb_cpu> <email> <qos> <partition> (<account>)"
    echo "contaminants: danio,human,mouse,mosquito,phi (add contaminants separated by comma)"
    echo "nb_cpu: max is 12 on tars"
    echo "qos: fast or normal or long"
    echo "partition: common or dedicated or common,dedicated with qos fast"
    exit
fi

mkdir -p $2

amplicon=$(readlink -e "$1")
outdir=$(readlink -e "$2")

account=""
if [ "$9" != "" ]
then
    account="#SBATCH --account=$9"
fi
SCRIPTPATH=$(dirname "${BASH_SOURCE[0]}")

echo """#!/bin/bash
#SBATCH --mail-user=$6
#SBATCH --mail-type=ALL
#SBATCH --qos=$7
#SBATCH -p $8
#SBATCH --cpus-per-task=$5
#SBATCH --mem=50000
#SBATCH --job-name="masque_$4"
$account
source /local/gensoft2/adm/etc/profile.d/modules.sh
module purge
export PATH=/pasteur/projets/policy01/Matrix/metagenomics/python-lib/bin:$PATH
export PYTHONPATH=/pasteur/projets/policy01/Matrix/metagenomics/python-lib/lib/python2.7/site-packages:$PYTHONPATH
module add Python/2.7.8 FastTree/2.1.8 FLASH/1.2.11 fasta mafft/7.149 bowtie2/2.2.9 blast+/2.2.40 AlienTrimmer/0.4.0 fastqc/0.11.5 rdp_classifier/2.12 BMGE/1.12 openmpi/2.0.1 IQ-TREE/1.5.1

/bin/bash $SCRIPTPATH/masque.sh -a $amplicon -o $outdir/ -t $5 -n $4 -c $3  &> $outdir/${4}_stat_process.txt || exit 1

exit 0
""">$SLURM_SCRIPT

SLURMID=`sbatch $SLURM_SCRIPT`

echo "Submission SLURM :> JOBID = $SLURMID"
