---
source: Rmd
title: Setup
teaching: XX
exercises: XX
---

:::::::::::::::::::::::::::::::::::::: questions 

- 

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Perform a quick-start analysis.
- Download the data set that will be used for the remaining episodes.

::::::::::::::::::::::::::::::::::::::::::::::::

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


::::::::::::::::::::::::::::::::::::: keypoints 

- When starting to use a new Bioconductor package, it is often good practice to go through the accompanying vignette and test some of the commands on the provided example data, before applying it to one's own data.

::::::::::::::::::::::::::::::::::::::::::::::::

