---
source: Rmd
title: Importing and annotating quantified data into R
teaching: XX
exercises: XX
---



::::::::::::::::::::::::::::::::::::::: objectives

- Learn how to import the quantifications into a SummarizedExperiment object.
- Learn how to add additional gene annotations to the object.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do we get our data into R?

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

### Contribute!

This episode is intended to show how we can assemble a SummarizedExperiment
starting from individual count, rowdata and coldata files. Moreover, we will
practice adding annotations for the genes, and discuss related concepts
and things to keep in mind (annotation sources, versions, 'helper' packages
like tximeta).


::::::::::::::::::::::::::::::::::::::::::::::::::

## Read the data

### Counts


```r
counts <- read.csv("data/GSE96870_counts_cerebellum.csv", 
                   row.names = 1)
```

### Sample annotations


```r
coldata <- read.csv("data/GSE96870_coldata_cerebellum.csv",
                    row.names = 1)
```

### Gene annotations

Need to be careful - the descriptions contain both commas and ' (e.g., 5')


```r
rowranges <- read.delim("data/GSE96870_rowranges.tsv", sep = "\t", 
                        colClasses = c(ENTREZID = "character"),
                        header = TRUE, quote = "", row.names = 5)
```

Mention other ways of getting annotations, and practice querying org package.
Important to use the right annotation source/version.


```r
suppressPackageStartupMessages({
    library(org.Mm.eg.db)
})
mapIds(org.Mm.eg.db, keys = "497097", column = "SYMBOL", keytype = "ENTREZID")
```

```{.output}
'select()' returned 1:1 mapping between keys and columns
```

```{.output}
497097 
"Xkr4" 
```

Check feature types


```r
table(rowranges$gbkey)
```

```{.output}

     C_region     D_segment          exon     J_segment      misc_RNA 
           20            23          4008            94          1988 
         mRNA         ncRNA precursor_RNA          rRNA          tRNA 
        21198         12285          1187            35           413 
    V_segment 
          535 
```

## Assemble SummarizedExperiment


```r
stopifnot(rownames(rowranges) == rownames(counts),
          rownames(coldata) == colnames(counts))

se <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    rowRanges = as(rowranges, "GRanges"),
    colData = coldata
)
```

## Save SummarizedExperiment


```r
saveRDS(se, "data/GSE96870_se.rds")
```

## Session info


```r
sessionInfo()
```

```{.output}
R version 4.2.2 Patched (2022-11-10 r83330)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.1 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] org.Mm.eg.db_3.16.0         AnnotationDbi_1.60.0       
 [3] SummarizedExperiment_1.28.0 Biobase_2.58.0             
 [5] MatrixGenerics_1.10.0       matrixStats_0.63.0         
 [7] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
 [9] IRanges_2.32.0              S4Vectors_0.36.1           
[11] BiocGenerics_0.44.0        

loaded via a namespace (and not attached):
 [1] compiler_4.2.2         XVector_0.38.0         bitops_1.0-7          
 [4] tools_4.2.2            zlibbioc_1.44.0        bit_4.0.5             
 [7] RSQLite_2.3.0          evaluate_0.20          memoise_2.0.1         
[10] lattice_0.20-45        pkgconfig_2.0.3        png_0.1-8             
[13] rlang_1.0.6            Matrix_1.5-3           DelayedArray_0.24.0   
[16] DBI_1.1.3              cli_3.6.0              xfun_0.37             
[19] fastmap_1.1.0          GenomeInfoDbData_1.2.9 httr_1.4.4            
[22] knitr_1.42             Biostrings_2.66.0      vctrs_0.5.2           
[25] bit64_4.0.5            grid_4.2.2             R6_2.5.1              
[28] blob_1.2.3             KEGGREST_1.38.0        renv_0.16.0           
[31] RCurl_1.98-1.10        cachem_1.0.6           crayon_1.5.2          
```

:::::::::::::::::::::::::::::::::::::::: keypoints

- Key point 1

::::::::::::::::::::::::::::::::::::::::::::::::::


