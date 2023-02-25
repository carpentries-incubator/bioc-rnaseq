---
source: Rmd
title: Introduction to RNA-seq
teaching: XX
exercises: XX
---



::::::::::::::::::::::::::::::::::::::: objectives

- To explain RNA-seq
- To provide sources of RNA-seq data - and typical research questions needing this approach.
- How does the RNA-seq data look like?
- To provide an overview of data handling procedures to answer biological questions.
- To provide an overview of common quality control issues and steps for the raw data.
- Explain how gene expression levels can be estimated from raw data.
- To determine gene expression level changes in more than one state/treatment
- To provide an overview of gene expression presentation/visualisation 
- To provide an overview of gene set analysis

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- What is RNA-seq?
- How does RNA-seq data look like?
- Where can the data be obtained?
- What are the typical research questions that would require an RNA-seq solution? Gives an overview of experimental design as well.
- What would be the general workflow for a typical RNA-seq analysis? This should also state the Tools/resources required.
- What are the typical quality issues encountered in RNA-seq analysis?
- How is gene expression estimated?
- In what functional roles are the expressed genes envolved in?

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

## Contribute!

This episode is intended to introduce important concepts in RNA-seq, such as the biology of RNA-seq, typical research questions, experimental design,  data processing, QC issues and ways around them, functional analysis of the expressed genes and use of Bioconductor for RNA-seq analysis pipelines to bring everyone up to speed.


::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: keypoints

- RNA-seq is a technique of measuring the amount of RNA expressed within a cell/tissue and state at a given time
- Sequencing procedures are used to determine the order of bases in each RNA
- The number of unique fragments gives gene counts
- The counts are normalised to factor in other quality parameter such as library size, etc
- Before processing of reads - they should be checked and filetered for quality.
- A number of pipelines can be used for reproducibility of the output (web based or command line.
- Calculated expression levels can be compared across states to get differential gene expression.
- The set of genes in any of the states at a given time can tell us about the functional /physiological state of the cell/tissue /etc

::::::::::::::::::::::::::::::::::::::::::::::::::


