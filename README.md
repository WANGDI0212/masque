# MASQUE : Metagenomic Analysis with a Quantitative pipeline 
[Amine Ghozlane](https://research.pasteur.fr/fr/member/amine-ghozlane/) (amine.ghozlane@pasteur.fr)
 
## Introduction

The aim of this project is to provide an easy cluster interface to perform targeted metagenomic analysis.
Masque allows :

* to analyse 16S/18S/23S/28S/ITS data. It builds a count matrix, an annotation table and a phylogeny of the OTU.
* to perform to use a set of parameters already tested on serveral projects for the numerous software used to perform the clustering and the annotation.
* to perform an "uptodate" analysis considering the scientific litterature.

## Installation

MASQUE comes with many binaries for Linux. It will always use your existing installed versions if they exist, but will use the included ones if that fails. You can consult the list of dependencies later in this document.

## Command line options

```
 masque.sh -i </path/to/input/directory/> -o </path/to/result/directory/>
- case high sensitive annotation: /pasteur/projets/policy01/Matrix/metagenomics/masque/masque.sh -i </path/to/input/directory/> -o </path/to/result/directory/> -b
- case 23S/28S: /pasteur/projets/policy01/Matrix/metagenomics/masque/masque.sh -i </path/to/input/directory/> -o </path/to/result/directory/> -l
- case its: /pasteur/projets/policy01/Matrix/metagenomics/masque/masque.sh -i </path/to/input/directory/> -o </path/to/result/directory/> -f
- case amplicon: /pasteur/projets/policy01/Matrix/metagenomics/masque/masque.sh -a <amplicon file> -o </path/to/result/directory/>
- All parameters:
-i      Provide </path/to/input/directory/>
-a      Provide <amplicon file>
-o      Provide </path/to/result/directory/>
-n      Indicate <project-name> (default: use the name of the input directory)
-t      Number of <thread> (default all cpu will be used)
-s      Perform OTU clustering with swarm
-b      Perform taxonomical annotation with blast (Default vsearch)
-l      Perform taxonomical annotation against LSU databases: Silva/RDP
-f      Perform taxonomical annotation against ITS databases: Unite/Findley/Underhill/RDP
--minreadlength Minimum read length take in accound in the study (Default 35nt)
--minphred      Qvalue must lie between [0-40] (Default minimum qvalue 20)
--minphredperc  Minimum allowed percentage of correctly called nucleotides [0-100] (Default 80)
--NbMismatchMapping     Maximum number of mismatch when mapping end-to-end against Human genome and Phi174 genome (Default 1 mismatch is accepted)
--maxoverlap    Minimum overlap when paired reads are considered (Default 200 nt)
--minoverlap    Maximum overlap when paired reads are considered (Default 50 nt)
--minampliconlength     Minimum amplicon length (Default 64nt)
--minotusize    Indicate minimum OTU size (Default 4)
--prefixdrep    Perform prefix dereplication (Default full length dereplication)
--chimeraslayerfiltering        Use ChimeraSlayer database for chimera filtering (Default : Perform a de novo chimera filtering)
--otudiffswarm  Number of difference accepted in an OTU with swarm (Default 1)
--evalueTaxAnnot        evalue threshold for taxonomical annotation with blast (Default evalue=1E-5)
--maxTargetSeqs Number of hit per OTU with blast (Default 1)
--identity_threshold    Identity threshold for taxonomical annotation with vsearch (Default 0.75)
```

## Results

In the output_dir, you will find after calculation :

File | Description
---|---
**project_stat_process.txt** | every step progress (during calculation : tail -f project-name_stat_process.txt, at the end : less project-name_stat_process.txt)
**project_annotation_process.tsv** | summary of the annotation process (Numbers of amplicon after dereplication, )
**project_build_process.tsv** | summary of the otu-build process (Number reads, contaminants and OTU identified per samples...)
**project_otu.fasta** | OTU centroid sequence 
**project_otu_table.tsv** | count table including the raw count obtained for each OTU and each sample
**project_vs_database_annotation_eval_val.tsv** | OTU annotation performed by blast against the several databank
**project_vs_rdp.tsv** | OTU annotation performed by rdp.
**project_otu.tree** | OTU phylogeny
**reads/*_fastqc.html** | fastq quality after trimming/clipping

The other files correspond to intermediate results.


## Dependencies

## Databases

## Citation

No papers about MASQUE alone will be published, but you can cite the first publication that use this program :
- A bacteriocin from epidemic Listeria strains alters the host intestinal microbiota to favor infection. Quereda JJ, Dussurget O, Nahori MA, Ghozlane A, Volant S, Dillies MA, Regnault B, Kennedy S, Mondot S, Villoing B, Cossart P, Pizarro-Cerda J.; PNAS 2016. [paper](http://www.ncbi.nlm.nih.gov/pubmed/27140611).

## Acknowledgements

Thanks to Emna Achouri  - emna.achouri@pasteur.fr for tars support.
