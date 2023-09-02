---
source: Rmd
title: Next steps
teaching: 20
exercises: 0
---

:::::::::::::::::::::::::::::::::::::: questions 

- How to go further from here?
- What other types of analyses can be done with RNA-seq data?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Get an overview of usages of RNA-seq data that are not covered in this workshop.

::::::::::::::::::::::::::::::::::::::::::::::::

## Transcript-level analyses

The analyses covered in this workshop all assumed that we have generated a matrix with read counts for each *gene*. 
Some questions, however, require expression estimates on a more fine-grained level, typically individual transcript isoforms. 
This is the case, for example, if we would like to look for differences in expression of individual isoforms, or changes in splicing patterns induced by a specific treatment. 
As already mentioned in episode 1, some quantification tools (such as Salmon, kallisto and RSEM) do in fact estimate abundances on the transcript level.
To perform transcript-level differential expression analysis, these estimates would be used directly, without aggregating them on the gene level. 
Some caution is warranted, however, as expression estimates for individual transcripts are often more noisy than after aggregation to the gene level.
This is due to the high similarity often observed between different isoforms of a gene, which leads to higher uncertainty when mapping reads to transcripts. 
For this reason, several approaches have been developed specifically for performing differential expression analysis on the transcript level [@Zhu2019-swish, @Baldoni2023-catchsalmon].

Analysis of differential splicing, or differential transcript usage, differs from the differential expression analyses mentioned so far as it is concerned with changes in the *relative* abundances of the transcripts of a gene.
Hence, considering the transcripts in isolation is no longer sufficient.
A good resource for learning more about differential transcript usage with Bioconductor is provided by @Love2018-swimming.

## De novo transcript assembly

In this lesson, we have assumed that we are working with a well-annotated transcriptome, that is, that we know the sequence of all expressed isoforms and how they are grouped together into genes. 
Depending on the organism or disease type we are interested in, this may or may not be a reasonable assumption. 
In situations where we do not believe that the annotated transcriptome is complete enough, we can use the RNA-seq data to *assemble* transcripts and create a custom annotation. 
This assembly can be performed either guided by a genome sequence (if one is available), or completely de novo. 
It should be noted that transcript assembly is a challenging task, which requires a deeply sequenced library to get the best results. 
In addition, data from more recent long-read sequencing technologies can be very helpful, as they are in principle able to sequence entire transcript molecules, thus circumventing the need for assembly.
Methods for transcript assembly represent an active area of research. 
Recent review are provided e.g. by @Raghavan2022-denovo and @Amarasinghe2020-longread.


::::::::::::::::::::::::::::::::::::: keypoints 

- RNA-seq data is very versatile and can be used for a number of different purposes. It is important, however, to carefully plan one's analyses, to make sure that enough data is available and that abundances for appropriate features (e.g., genes, transcripts, or exons) are quantified.

::::::::::::::::::::::::::::::::::::::::::::::::

