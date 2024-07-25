# De novo Transcriptome Assembly of Spirodela polyrhiza

## Description

This project involves the de novo assembly and functional annotation of the transcriptome of *Spirodela polyrhiza* grown under blue light conditions, as described in Zhong et al. (2022). Each participant used different seeds to compare the results, using a higher number of reads than in previous tasks.

## Background

The type of light (red, blue, and white) with which some plants such as *Spirodela polyrhiza* are irradiated affects functions such as growth and certain physiological processes by producing responses in the transcriptome. This is of significant interest for the various applications that *S. polyrhiza* may have, especially in the area of energy sources.

## Objective

The aim of this work is to perform de novo assembly of the transcriptome of the duckweed species, *Spirodela polyrhiza*, and annotate its functions.

## Materials and Methods

### Experimental Design and Workflow

According to the reference paper (Zhong et al. 2022), the SRA accession number from which the data is obtained is SRR15716419. A seed based on the sum of the digits of a group member's DNI is used to ensure reproducibility. The workflow consists of the following steps:

1. **Subsampling**: Extract 20,000 reads from SRR15716419.sra using a specific seed.
2. **Preprocessing**: Quality control, error correction, adapter removal, and rRNA elimination.
3. **Assembly**: De novo transcriptome assembly using Trinity.
4. **Quality Check**: Evaluate the assembly quality using various tools.
5. **Functional Annotation**: Identify coding regions, predict signal peptides, and search for sequence homologies.

### Commands Used

Here are the commands used for the analysis:

```bash
# Download data
SRA=SRR15716419
prefetch -p -v $SRA
fasterq-dump -e 4 --progress --split-3 ./SRR15716419.sra

# Subsampling
seqtk sample -s 38 SRR15716419_1.fastq 20000 > sub_SRR15716419_1.fastq
seqtk sample -s 38 SRR15716419_2.fastq 20000 > sub_SRR15716419_2.fastq
sed -i -r 's/(^[\@\+]SRR\S+)/\1\/1/' sub_SRR15716419_1.fastq
sed -i -r 's/(^[\@\+]SRR\S+)/\1\/2/' sub_SRR15716419_2.fastq

# Quality Verification
fastqc --outdir=./quality_check sub_SRR15716419_1.fastq sub_SRR15716419_2.fastq

# Error Correction
perl ~/denovo/software/run_rcorrector.pl -1 sub_SRR15716419_1.fastq -2 sub_SRR15716419_2.fastq

# Adapter and Low-Quality Sequence Removal
mkdir trimmed_data
cp ~/denovo/software/TruSeq3-PE.fa .
java -jar ~/denovo/software/trimmomatic-0.39.jar PE -threads 2 -trimlog ./trimmed_data/trimmomatic.log sub_SRR15716419_1.cor.fq sub_SRR15716419_2.cor.fq ./trimmed_data/sub_SRR15716419_1P.qtrim.fq.gz ./trimmed_data/sub_SRR15716419_1U.qtrim.fq.gz ./trimmed_data/sub_SRR15716419_2P.qtrim.fq.gz ./trimmed_data/sub_SRR15716419_2U.qtrim.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:15 TRAILING:15 MINLEN:130

# rRNA Elimination
sortmerna --ref ~/denovo/denovo_DB/rRNA_db/rfam-5.8s-database-id98.fasta --ref ~/denovo/denovo_DB/rRNA_db/silva-euk-18s-id95.fasta --ref ~/denovo/denovo_DB/rRNA_db/silva-euk-28s-id98.fasta --reads ./trimmed_data/sub_SRR15716419_1P.qtrim.fq.gz --reads ./trimmed_data/sub_SRR15716419_2P.qtrim.fq.gz --aligned rRNA --fastx --other ./cleanedrRNA.fastq --threads 2 -v --out2 --workdir ./rmRNA

# Quality Check
fastqc -t 5 --outdir=./quality_check cleanedrRNA.fastq_fwd.fq.gz cleanedrRNA.fastq_rev.fq.gz

# Normalization
TRINITYPATH=/home/fran/opt/trinityrnaseq-Trinity-v2.8.4
$TRINITYPATH/util/insilico_read_normalization.pl --seqType fq --CPU 15 --JM 60G --max_cov 20 --left cleanedrRNA.fastq_fwd.fq.gz --right cleanedrRNA.fastq_rev.fq.gz --output ./normalization_reads --pairs_together --PARALLEL_STATS

# Assembly with TRINITY
mkdir trinity_out
TRINITY_DIR=/home/fran/opt/trinityrnaseq-Trinity-v2.8.4
DATPATH=$HOME/grupo5/tarea3/normalization_reads
$TRINITY_DIR/Trinity --seqType fq --no_normalize_reads --left $DATPATH/left.norm.fq --right $DATPATH/right.norm.fq --SS_lib_type RF --max_memory 2G --CPU 2 --output ./trinity_out

## Results

### Summary Statistics
- The small average size of the reads (~400) suggests low-quality assembly.
- Low mapping percentage indicates fragmented transcriptome.
- Only a few complete proteins were found in the transcripts, indicating low coverage.
- BUSCO analysis shows 98.2% missing values, confirming poor quality reads.

### Functional Annotation
- Only 20 protein sequences were found with homology to the SWISSPROT database.
- High scoring proteins belong to the V-type aminotransferase family and ATPases.
- Signal peptide prediction shows a type 1 signal peptide with high probability.

## Conclusions
This task involved the use of short reads for de novo assembly of the Spirodela polyrhiza transcriptome. The results indicate a low-quality transcriptome, with poor representation of reads and almost no genes found as per BUSCO analysis. Functional annotation showed limited results with BLAST, identifying a few highly expressed transcripts under the study conditions.

Comparing this analysis with the class example on Helianthus annuus, both showed low-quality transcriptomes due to the low number of reads used, resulting in poor assembly and functional annotation outcomes.

## References
Zhong, Y. et al., 2022. Physiological responses and transcriptome analysis of Spirodela polyrhiza under red, blue, and white light. Planta 255, 1–15. [https://doi.org/10.1007/s00425-021-03764-4](https://doi.org/10.1007/s00425-021-03764-4)
