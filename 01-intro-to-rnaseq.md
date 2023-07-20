---
source: Rmd
title: Introduction to RNA-seq
teaching: XX
exercises: XX
---



::::::::::::::::::::::::::::::::::::::: objectives

- Explain what RNA-seq is.
- Describe some of the most common design choices that have to be made before running an RNA-seq experiment.
- Provide an overview of the procedure to go from the raw data to the read count matrix that will be used for downstream analysis.
- Show some common types of results and visualizations generated in RNA-seq analyses.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- What are the different choices to consider when planning an RNA-seq experiment? 
- How does one process the raw fastq files to generate a table with read counts per gene and sample?
- Where does one find information about annotated genes for a given organism?
- Which are the most commonly used Bioconductor packages for statistical analysis of RNA-seq data?
- What are the typical steps in an RNA-seq analysis?


::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

## Contribute!

This episode is intended to introduce important concepts in RNA-seq, such as the biology of RNA-seq, typical research questions, experimental design,  data processing, QC issues and ways around them, functional analysis of the expressed genes and use of Bioconductor for RNA-seq analysis pipelines to bring everyone up to speed.  See this [link](https://docs.google.com/document/d/12hUqVo2MhgYH9IT8ShJKSZgeMLzjeUh2EKdxM_siWwM) on how to contribute


::::::::::::::::::::::::::::::::::::::::::::::::::


# RNA-seq quantification: from reads to count matrix

There is a plethroa of RNA quantification pipelines, and the most common approaches can be categorized into three main types:

1. Align reads to the genome, and count the number of reads that map to each gene.
   This is the one of simplest methods. For species for which the transcriptome is poorly annotated, this would be the preferred approach.
   Example: `STAR` alignment to GRCz11 + `Rsubread` [featureCounts](https://doi.org/10.1093%2Fnar%2Fgkz114)

2. Align reads to the transcriptome, quantify transcript expression, and summarize transcript expression into gene expression.
   This approach can produce accurate quantification results based [independent benchmarking](https://doi.org/10.1186/s13059-016-0940-1), 
   particularly for high-quality samples without DNA contamination.
   Example: RSEM quantification using `rsem-calculate-expression --star` on the GENCODE GRCh38 transcriptome + `tximport`

3. Pseudoalign reads against the transcriptome, using the corresponding genome as a decoy, quantifying transcript expression in the process, 
   and summarize the transcript-level expression into gene-level expression.
   The advantages of this approach include: computational efficiency, mitigation of the effect of DNA contamination, and GC bias correction.
   Example: `salmon quant --gcBias` + `tximport`
   
At typical sequencing read depth, gene expression quantification is often more accurate than transcript expression quantification.
However, gene expression quantification can be [improved](https://doi.org/10.12688/f1000research.7563.1)
by first quantifying transcript expression and then summarizing it to gene expression.

Other tools used in RNA-seq quantification include: TopHat2, bowtie2, kallisto, HTseq, among many others.

Choosing the appropriate RNA-seq quantification would depend on the quality of the transcriptome annotation,
the quality of the RNA-seq library preparation, the presence of contaminating sequences, among many factors.
Often time, it would be important to compare the quantification results of multiple approaches.


:::::::::::::::::::::::::::::::::::::::: keypoints

- RNA-seq is a technique of measuring the amount of RNA expressed within a cell/tissue and state at a given time.
- Many choices have to be made when planning an RNA-seq experiment, such as whether to perform poly-A selection or ribosomal depletion, whether to apply a stranded or an unstranded protocol, and whether to sequence the reads in a single-end or paired-end fashion. Each of the choices have consequences for the processing and interpretation of the data. 
- Many approaches exist for quantification of RNA-seq data. Some methods align reads to the genome and count the number of reads overlapping gene loci. Other methods map reads to the transcriptome and use a probabilistic approach to estimate the abundance of each gene or transcript. 
- Information about annotated genes can be accessed via several sources, including Ensembl, UCSC and GENCODE. 
- The most commonly used Bioconductor packages for statistical analysis of RNA-seq data are `DESeq2`, `edgeR` and `limma`. 


::::::::::::::::::::::::::::::::::::::::::::::::::


