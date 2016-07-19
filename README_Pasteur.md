# MASQUE
[Amine Ghozlane](https://research.pasteur.fr/fr/member/amine-ghozlane/) (amine.ghozlane@pasteur.fr)

## Introduction

The aim of this project is to provide an easy cluster interface to perform targeted metagenomic analysis.
Masque allows :

* to analyse 16S/18S/23S/28S/ITS data. It builds a count matrix, an annotation table and a phylogeny of the OTU.
* to perform to use a set of parameters already tested on serveral projects for the numerous software used to perform the clustering and the annotation.
* to perform an "uptodate" analysis considering the scientific litterature.

## Load MASQUE on Institut Pasteur cluster's

```
module use /pasteur/projets/Matrix/modules
module add masque/bic or masque/tars or masque/standalone
```

* masque/bic is supposed to be used on bic cluster
* masque/tars is supposed to be used on tars cluster
* masque/standalone is available on both bic and tars cluster and provide an access to every masque parameters. You need to write your own submission script to use masque

## Use MASQUE on Institut Pasteur cluster's

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
- <email>: your email address (it will send a start|error|success message)
- <queue>: On bic, it correspond to your team group on the cluster. On tars, you should indicated if it's a short job : "fast" (less than 2 hours, high priority), "normal" (less than 24 hours, normal priority) or unlimited (no time limit, 5 job per user)  

Masque will submit the automatically calculation on the cluster. You don't have to write any submission script. Masque must be used on the master node that you reach automatically when you perform:

- bic: access by sending an email to informatique at pasteur.fr
'''
ssh user@bic.pasteur.fr
'''
- tars: access by sending an email to informatique at pasteur.fr, plus perform every exercise [here](https://moocs.pasteur.fr/courses/Institut_Pasteur/DSI_01/1/info)
'''
ssh user@tars.pasteur.fr
'''

The calculation status can be considered with :
- bic: `qstat`
- tars:`squeue -u username`

## Results

## Results

In output_dir, you will find after calculation :

- project-name_stat_process.txt : every step progress (during calculation : tail -f project-name_stat_process.txt, at the end : less project-name_stat_process.txt)
- project-name_annotation_process.tsv : summary of the annotation process (Numbers of amplicon after dereplication, )
- project-name_build_process.tsv : summary of the otu-build process (Number reads, contaminants and OTU identified per samples...)
- project-name_otu.fasta : otu centroid sequence 
- project-name_otu_table.tsv : count table including the raw count obtained for each OTU and each sample
- project-name_vs_database-name_annotation_eval_1E-5.tsv : OTU annotation performed by blast against the several databank
- project-name_vs_rdp.tsv : OTU annotation performed by rdp.
- project-name_otu.tree : OTU phylogeny
- reads/*_fastqc.html : fastq quality after trimming/clipping

The other files correspond to intermediate results.

If some of theses files are missing, it might be worth to re-run the program with the same command line. masque check each step and can continue after the last step. It happens sometimes that one cluster node fail for no particular reason.
The project-name_otu_table.tsv and project-name_vs_database-name_annotation_eval_1E-5.tsv  can be directly analysed with [SHAMAN](http://shaman.c3bi.pasteur.fr/).

## Citation

No papers about MASQUE alone will be published, but you can cite the first publication that use this program :
- A bacteriocin from epidemic Listeria strains alters the host intestinal microbiota to favor infection. Quereda JJ, Dussurget O, Nahori MA, Ghozlane A, Volant S, Dillies MA, Regnault B, Kennedy S, Mondot S, Villoing B, Cossart P, Pizarro-Cerda J.; PNAS 2016. [paper](http://www.ncbi.nlm.nih.gov/pubmed/27140611).

## Acknowledgements

Thanks to Emna Achouri  - emna.achouri@pasteur.fr for tars support.
