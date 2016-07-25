# MASQUE : Metagenomic Analysis with a Quantitative pipeline 
[Amine Ghozlane](https://research.pasteur.fr/fr/member/amine-ghozlane/) (amine.ghozlane@pasteur.fr) (@xealf8)

## Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Command line options](#command-line-options)
- [SGE and SLURM deployments](#sge-and-slurm-deployments)
- [Results](#results)
- [Dependencies](#dependencies)
- [Databases](#databases)
- [Test](#test)
- [Citation](#citation)
- [Acknowledgements](#acknowledgements)

 
## Introduction

The aim of this project is to provide an easy cluster interface to perform targeted metagenomic analysis.
MASQUE allows :

* to analyse 16S/18S/23S/28S/ITS data. It builds a count matrix, an annotation table and a phylogeny of the OTU.
* to perform to use a set of parameters already tested on serveral projects for the numerous software used to perform the clustering and the annotation.
* to perform an "uptodate" analysis considering the scientific litterature.

## Installation

MASQUE comes with many binaries for Linux 64 bits. It will always use your existing installed versions if they exist, but will use the included ones if that fails. 
You can consult the list of dependencies later in this document. For the correct deployment by git, install first [git-lfs](https://git-lfs.github.com/). 
```
sudo ./git-lfs-1.2.1/install.sh
git lfs install
```
Then, you can clone masque :
```
git clone https://github.com/aghozlane/masque.git
```
Only biom program need to be installed by the user :  
```
pip install biom
```
Then, install the databases as follow :  
```
/bin/bash install_databases.sh
```

## Command line options

```
 /bin/bash masque.sh -i </path/to/input/directory/> -o </path/to/result/directory/>
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

## SGE and SLURM deployments

Template scripts are provided for SGE and SLURM deployments :  

masque_16S-18S.sh  
masque_16S-18S_tars.sh  
masque_23S-28S.sh  
masque_23S-28S_tars.sh  
masque_ITS.sh  
masque_ITS_tars.sh  
masque_amplicon.sh  
masque_amplicon_tars.sh  

For users from Institut Pasteur, please consider the [README_PASTEUR](README_Pasteur.md).


## Results

In the output_dir, you will find after calculation :

File | Description
---|---
**project_stat_process.txt** | Every step progress (during calculation : tail -f project-name_stat_process.txt, at the end : less project-name_stat_process.txt)
**project_annotation_process.tsv** | Summary of the annotation process
**project_build_process.tsv** | Summary of the otu-build process (Number reads, contaminants and OTU identified per samples...)
**project_otu.fasta** | OTU centroid sequence in fasta format 
**project_otu_table.tsv** | Count table including the raw count obtained for each OTU and each sample
**project_vs_database_annotation_eval_val.tsv** | OTU annotation performed by blast against the several databank
**project_database_eval_val.biom** | Biom file including the count and the annotation
**project_vs_rdp.tsv** | OTU annotation performed by rdp.
**project_otu.tree** | OTU phylogeny
**reads/*_fastqc.html** | fastq quality after trimming/clipping

The other files correspond to intermediate results.

## Dependencies


*  __AlienTrimmer__  
    Performs the trimming and clipping of the reads  
    _Criscuolo, A., Brisse, S., AlienTrimmer: a tool to quickly and accurately trim off multiple short contaminant sequences from high-throughput sequencing reads., Genomics, 2013, 102(5), 500-506._  
*  __Biom__  
   Combine count matrix and taxonomical annotation table  
   _Daniel McDonald, et al., The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome. GigaScience, 2012, 1:7. doi:10.1186/2047-217X-1-7_  
*  __Blastn__  
   Performs the taxonomical annotation  
   _Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. "Basic local alignment search tool." J. Mol. Biol., 1990, 215:403-410._  
*  __Bowtie2__  
   Finds contaminants  
   _Langmead B, Salzberg S0,. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359._  
*  __Fastqc__  
   Checks read quality.  
   _Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data._  
*  __Fasttree__  
   Computes the phylogeny of OTU sequences.  
   _Price, M.N., Dehal, P.S., and Arkin, A.P., FastTree 2 -- Approximately Maximum-Likelihood Trees for Large Alignments. PLoS ONE, 2010, 5(3):e9490._  
*  __FLASH__  
   Merges paired reads to get amplicons.  
   _Magoc T. and Salzberg S.., FLASH: Fast length adjustment of short reads to improve genome assemblies. Bioinformatics 27:21, 2011, 2957-63._  
*  __Mafft__  
   Performs a multiple alignment of OTU sequences.  
   _Katoh, Misawa, Kuma, Miyata, MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Res., 2002, 30:3059-3066._  
*  __rdp classifier__  
   Performs taxonomical annotation for 16S, 18S, ITS.  
   _Wang, Q, Garrity G.M., Tiedje J. M., and  ColeJ. R., Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl Environ Microbiol., 2007, 73(16):5261-5267._  
*  __swarm__  
   Performs OTU clustering.  
   _Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M., Swarm v2: highly-scalable and high-resolution amplicon clustering. PeerJ, 2015._  
*  __vsearch__  
   Performs OTU clustering (Default selection).  
   _Torbjørn Rognes https://github.com/torognes/vsearch_  


## Databases

MASQUE use several databases for taxonomical annotation and data filtering as follow :

### Taxonomical annotation 

*  __FINDLEY__  
   Used for the taxonomical annotation of ITS sequences.  
   _Findley, K., et al., Topographic diversity of fungal and bacterial communities in human skin. Nature, 2013, 498(7454), 367-370._  
   http://www.mothur.org/wiki/Findley_ITS_Database  
*  __GREENGENES__  
    Used for the taxonomical annotation of 16S, 18S sequences.  
    _DeSantis, T. Z., et al. Greengenes, a chimera-checked 16S rRNA gene database and workbench compatible with ARB. Applied and environmental microbiology, 2006, 72(7), 5069-5072._  
   http://greengenes.secondgenome.com/downloads/database/13_5  
*  __SILVA LSU, SSU__  
   Used for the taxonomical annotation of 16S, 18S, 23S, 28S sequences.  
   _Pruesse, E., et al., SILVA: a comprehensive online resource for quality checked and aligned ribosomal RNA sequence data compatible with ARB. Nucleic acids research, 2007, 35(21), 7188-7196._  
   https://www.arb-silva.de/  
*  __UNDERHILL__  
   Used for the taxonomical annotation of ITS sequences.  
   _Tang J, Iliev I, Brown J, Underhill D and Funari V. Mycobiome: Approaches to Analysis of Intestinal Fungi. Journal of Immunological Methods, 2015, 421:112-21._  
   https://risccweb.csmc.edu/microbiome/thf/  
*  __UNITE__  
   Used for the taxonomical annotation of ITS sequences.  
   _Abarenkov, K., et al., The UNITE database for molecular identification of fungi–recent updates and future perspectives. New Phytologist, 2010, 186(2), 281-285._  
   https://unite.ut.ee/repository.php  

### Filtering databases

*  __AlienTrimmer__  
   Adapters sequences provided by Illumina and Life technologies.  
*  __GOLD__  
   ChimeraSlayer reference database used for chimera filtering (default mode use de novo filtering instead)  
   http://drive5.com/uchime/uchime_download.html  
*  __NCBI Homo Sapiens, PhiX174__  
   Used to search Human and PhiX174 contaminants (manipulators and Phi phage used for the sequencing)  

## Test

Samples from the mock communities are available for testing \[[The NIH HMP Working Group, 2009](http://www.ncbi.nlm.nih.gov/pubmed/19819907)\].

Label |  Community
---|---
**SRR053818** | even
**SRR072220** | even
**SRR072221** | even
**SRR072223** | staggered
**SRR072237** | staggered
**SRR072239** | staggered

Mock communities are composed of 21 species mixed in even or staggered proportions : 
<img src="test/mock/mock.png" align="center" />

You can run:
```
gunzip test/result/*.gz
/bin/bash ./masque.sh -i test/data -o test/result
```

The results can be visualized with [SHAMAN](http://shaman.c3bi.pasteur.fr/) and compared with the results obtained with the MOCK reference genome available in test/mock/.


## Citation

No papers about MASQUE alone is published for the moment, but you can cite the first publication that use this program:  
- A bacteriocin from epidemic Listeria strains alters the host intestinal microbiota to favor infection. Quereda JJ, Dussurget O, Nahori MA, Ghozlane A, Volant S, Dillies MA, Regnault B, Kennedy S, Mondot S, Villoing B, Cossart P, Pizarro-Cerda J.; PNAS 2016. [PUBMED](http://www.ncbi.nlm.nih.gov/pubmed/27140611).

## Acknowledgements

Thanks to Emna Achouri  - emna.achouri@pasteur.fr for tars support.
