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

- How can one import quantified gene expression data into an object suitable for downstream statistical analysis in R?
- What types of gene identifiers are typically used, and how are mappings between them done? 


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
R version 4.3.1 (2023-06-16)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] org.Mm.eg.db_3.17.0         AnnotationDbi_1.62.2       
 [3] SummarizedExperiment_1.30.2 Biobase_2.60.0             
 [5] MatrixGenerics_1.12.2       matrixStats_1.0.0          
 [7] GenomicRanges_1.52.0        GenomeInfoDb_1.36.1        
 [9] IRanges_2.34.1              S4Vectors_0.38.1           
[11] BiocGenerics_0.46.0        

loaded via a namespace (and not attached):
 [1] Matrix_1.5-4.1          bit_4.0.5               compiler_4.3.1         
 [4] BiocManager_1.30.21     renv_1.0.0              crayon_1.5.2           
 [7] blob_1.2.4              Biostrings_2.68.1       bitops_1.0-7           
[10] png_0.1-8               fastmap_1.1.1           yaml_2.3.7             
[13] lattice_0.21-8          R6_2.5.1                XVector_0.40.0         
[16] S4Arrays_1.0.4          knitr_1.43              DelayedArray_0.26.6    
[19] GenomeInfoDbData_1.2.10 DBI_1.1.3               rlang_1.1.1            
[22] KEGGREST_1.40.0         cachem_1.0.8            xfun_0.39              
[25] bit64_4.0.5             RSQLite_2.3.1           memoise_2.0.1          
[28] cli_3.6.1               zlibbioc_1.46.0         grid_4.3.1             
[31] rstudioapi_0.15.0       vctrs_0.6.3             evaluate_0.21          
[34] RCurl_1.98-1.12         httr_1.4.6              pkgconfig_2.0.3        
[37] tools_4.3.1            
```

:::::::::::::::::::::::::::::::::::::::: keypoints

- Depending on the gene expression quantification tool used, there are different ways (often distributed in Bioconductor packages) to read the output into a SummarizedExperiment or DGEList object for further processing in R.
- Stable gene identifiers such as Ensembl or Entrez IDs should preferably be used as the main identifiers throughout an RNA-seq analysis, with gene symbols added for easier interpretation.


::::::::::::::::::::::::::::::::::::::::::::::::::


