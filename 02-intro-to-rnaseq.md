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

This episode is intended to introduce important concepts in RNA-seq, such as the biology of RNA-seq, typical research questions, experimental design,  data processing, QC issues and ways around them, functional analysis of the expressed genes and use of Bioconductor for RNA-seq analysis pipelines to bring everyone up to speed.


::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: keypoints

- RNA-seq is a technique of measuring the amount of RNA expressed within a cell/tissue and state at a given time.
- Many choices have to be made when planning an RNA-seq experiment, such as whether to perform poly-A selection or ribosomal depletion, whether to apply a stranded or an unstranded protocol, and whether to sequence the reads in a single-end or paired-end fashion. Each of the choices have consequences for the processing and interpretation of the data. 
- Many approaches exist for quantification of RNA-seq data. Some methods align reads to the genome and count the number of reads overlapping gene loci. Other methods map reads to the transcriptome and use a probabilistic approach to estimate the abundance of each gene or transcript. 
- Information about annotated genes can be accessed via several sources, including Ensembl, UCSC and GENCODE. 
- The most commonly used Bioconductor packages for statistical analysis of RNA-seq data are `DESeq2`, `edgeR` and `limma`. 


::::::::::::::::::::::::::::::::::::::::::::::::::


