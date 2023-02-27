---
source: Rmd
title: Introduction and setup
teaching: XX
exercises: XX
---



::::::::::::::::::::::::::::::::::::::: objectives

- Install required packages.
- Download the data.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- 
::::::::::::::::::::::::::::::::::::::::::::::::::

## Download data

The data we will use in this lesson is obtained from the [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/), accession number [GSE96870](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96870).


```r
dir.create("data", showWarnings = FALSE)
download.file(
    url = "https://github.com/Bioconductor/bioconductor-teaching/blob/master/data/GSE96870/GSE96870_counts_cerebellum.csv?raw=true", 
    destfile = "data/GSE96870_counts_cerebellum.csv"
)

download.file(
    url = "https://github.com/Bioconductor/bioconductor-teaching/blob/master/data/GSE96870/GSE96870_coldata_cerebellum.csv?raw=true", 
    destfile = "data/GSE96870_coldata_cerebellum.csv"
)

download.file(
    url = "https://github.com/Bioconductor/bioconductor-teaching/blob/master/data/GSE96870/GSE96870_coldata_all.csv?raw=true", 
    destfile = "data/GSE96870_coldata_all.csv"
)

download.file(
    url = "https://github.com/Bioconductor/bioconductor-teaching/blob/master/data/GSE96870/GSE96870_rowranges.tsv?raw=true", 
    destfile = "data/GSE96870_rowranges.tsv"
)
```

:::::::::::::::::::::::::::::::::::::::: keypoints

- Key point 1

::::::::::::::::::::::::::::::::::::::::::::::::::


