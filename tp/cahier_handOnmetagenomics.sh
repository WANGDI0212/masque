# WORKDIR
WORKDIR=$HOME/masque/tp/
# DATADIR
DATADIR=$HOME/masque/test/data/

################
# 16S analysis #
################
# Even samples from MOCK communities
sample=SRR072220
sample=SRR072239
# Staggered samples from MOCK communities
sample=SRR072221
sample=SRR072223
sample=SRR072237
# Filtering quality
$HOME/masque/vsearch_bin/bin/vsearch --fastq_filter $DATADIR/${sample}.fastq --fastqout $WORKDIR/${sample}_filt.fastq --fastq_truncqual 16 --fastq_trunclen 250
# Dereplication
$HOME/masque/fastq2fasta/fastq2fasta.py -i $WORKDIR/${sample}_filt.fastq -s $sample -o $WORKDIR/${sample}_filt.fasta
$HOME/masque/vsearch_bin/bin/vsearch --derep_fulllength $WORKDIR/${sample}_filt.fasta -output $WORKDIR/${sample}_drep.fasta -sizeout
# Remove singleton
$HOME/masque/vsearch_bin/bin/vsearch -sortbysize $WORKDIR/${sample}_drep.fasta -output $WORKDIR/${sample}_nosing.fasta  -minsize 2
# Remove chimera
$HOME/masque/vsearch_bin/bin/vsearch --uchime_denovo $WORKDIR/${sample}_nosing.fasta --chimeras $WORKDIR/${sample}_chim.fasta --nonchimeras $WORKDIR/${sample}_nochim.fasta
# Clustering
$HOME/masque/vsearch_bin/bin/vsearch --cluster_size $WORKDIR/${sample}_nochim.fasta --id 0.97 --centroids $WORKDIR/${sample}_otu.fasta --sizein --relabel OTU_
# Mapping
$HOME/masque/vsearch_bin/bin/vsearch -usearch_global $WORKDIR/${sample}_filt.fasta -db $WORKDIR/${sample}_otu.fasta   --id 0.97 -uc $WORKDIR/${sample}_map.txt
# Count matrix
python $HOME/masque/usearch_python_scripts/uc2otutab.py $WORKDIR/${sample}_map.txt > $WORKDIR/${sample}_otu_table.tsv
# Annotation MOCK
$HOME/masque/vsearch_bin/bin/vsearch --usearch_global $WORKDIR/${sample}_otu.fasta --db $DATADIR/HMP_MOCK_v35_annotated.fasta --id 0.9  --blast6out $WORKDIR/${sample}_vs_mock.tsv --alnout $WORKDIR/${sample}_vs_mock_ali.txt
$HOME/masque/get_taxonomy/get_taxonomy.py -i $WORKDIR/${sample}_vs_mock.tsv -d $DATADIR/HMP_MOCK_v35_annotated.fasta -u $WORKDIR/${sample}_otu.fasta -o $WORKDIR/${sample}_annotation_mock.tsv -ob $WORKDIR/${sample}_annotation_mock_biom.txt
# Biom
biom convert -i $WORKDIR/${sample}_otu_table.tsv -o $WORKDIR/${sample}.biom --table-type="OTU table" --to-json
biom add-metadata -i $WORKDIR/${sample}.biom -o $WORKDIR/${sample}_annotated.biom --observation-metadata-fp $WORKDIR/${sample}_annotation_mock_biom.txt --observation-header id,taxonomy --sc-separated taxonomy
