# MASQUE on the clusters of the Institut Pasteur 
[Amine Ghozlane](https://research.pasteur.fr/fr/member/amine-ghozlane/) (amine.ghozlane@pasteur.fr)

## Introduction

Targeted metagenomics are for the moment hard to be performed on a simple laptop when we have hundreds of samples with a million of reads. These constraints lead us to supercomputers to perform the calculation in a short amount of time and without memory considerations.
We describe here the use of MASQUE on the clusters of the Institut Pasteur. The Institut Pasteur possess two clusters in production bic and tars. tars is the new cluster and bic will be closed in the futur month.
Please referer to the [README file](https://github.com/aghozlane/masque) for more information on what perform MASQUE.

## Load MASQUE

```
module use /pasteur/projets/Matrix/modules
module add masque/bic or masque/tars or masque/standalone
```

* masque/bic is supposed to be used on bic cluster
* masque/tars is supposed to be used on tars cluster
* masque/standalone is available on both bic and tars cluster and provide an access to every masque parameters. You need to write your own submission script or use interactive access to use it.

## Use MASQUE

Three programs will then be available with the following options :
- masque_16S-18S <reads_dir> <output_dir> <project-name> <nb_cpu> <email> <queue>
- masque_23S-28S <reads_dir> <output_dir> <project-name> <nb_cpu> <email> <queue>
- masque_ITS <reads_dir> <output_dir> <project-name> <nb_cpu> <email> <queue>
- masque_amplicon <amplicon_file> <output_dir> <project-name> <nb_cpu> <email> <queue>
With masque/standalone, you access to masque main program. Please consider the [README file](https://github.com/aghozlane/masque) to use it properly. 

which correspond to :
- <reads_dir>: directory with the read files (*.fq or *.fastq, if paired : *_R1.fastq, *_R2.fastq or *_R1.fq, *_R2.fq)
- <output_dir>: directory where to put the results
- <project-name>: a name for the project
- <nb_cpu>: Number of cpu used, max on bic is 12.
- <email>: your email address (it will send a start|error|success message). You have to provide a pasteur mail.
- <queue>: On bic, it corresponds to your team group on the cluster. On tars, you should indicated if it is a short job : "fast" (less than 2 hours, high priority), "normal" (less than 24 hours, normal priority) or unlimited (no time limit, 5 job per user)  

Masque will submit the automatically calculation on the cluster. You don't have to write any submission script. Masque must be used on the master node that you reach automatically when you perform:

- bic: access by sending an email to informatique at pasteur.fr
`ssh user@bic.pasteur.fr`

- tars: access by sending an email to informatique at pasteur.fr, plus perform every exercise [here](https://moocs.pasteur.fr/courses/Institut_Pasteur/DSI_01/1/info)
`ssh user@tars.pasteur.fr`

The calculation status can be considered with :
- bic: `qstat`
- tars:`squeue -u username`
