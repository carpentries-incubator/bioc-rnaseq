---
title: "Importing and annotating quantified data into R"
source: Rmd
teaching: XX
output:
  html_document:
    df_print: paged
exercises: XX
---


```{.warning}
Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'
```

::::::::::::::::::::::::::::::::::::::: objectives
-   Learn how to import the quantifications into a SummarizedExperiment object.
-   Learn how to add additional gene annotations to the object.
::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions
-   How can one import quantified gene expression data into an object suitable for downstream statistical analysis in R?
-   What types of gene identifiers are typically used, and how are mappings between them done?
::::::::::::::::::::::::::::::::::::::::::::::::::

## Load libraries

```r
suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.Mm.eg.db)
    library(hgu95av2.db)
    library(SummarizedExperiment)
})
```

## Load data

In the last episode, we used R to download 4 files from the internet and saved them on our computer. But we do not have these files loaded into R yet so that we can work with them. The original experimental design in [Blackmore et al. 2017](https://pubmed.ncbi.nlm.nih.gov/28696309/) was fairly complex: 8-12 weeks old male and female mice, with two different genetic backgrounds (wild-type and D2D), were collected at Day 0 (before influenza infection), Day 4 and Day 8 after influenza infection. From each mouse, cerebellum and spinal cord tissues were taken for RNA-seq. There were originally 4 mice per 'Sex X Day x Genotype' group, but a few were lost along the way resulting in a total of 45 samples. For this workshop, we are going to simplify the analysis by only using the 22 cerebellum samples. Expression quantification was done using STAR to align to the mouse genome and then counting reads that map to genes. In addition to the counts per gene per sample, we also need information on which sample belongs to which Sex/Time point/Replicate. And for the genes, it is helpful to have extra information called annotation.
Let's read in the data files that we downloaded in the last episode:


### Counts


```r
counts <- read.csv("data/GSE96870_counts_cerebellum.csv", 
                   row.names = 1)
dim(counts)
```

```{.output}
[1] 41786    22
```

```r
# View(counts)
```

Genes are in rows and samples are in columns, so we have counts for 41,786 genes and 22 samples. The `View()` command has been comment out for the website, but running it will open a tab in RStudio that lets us look at the data and even sort the table by a particular column. However, the viewer cannot change the data inside the `counts` object, so we can only look, not permanently sort nor edit the entries. When finished, close the viewer using the X in the tab. Looks like the rownames are gene symbols and the column names are the GEO sample IDs, which are not very informative for telling us which sample is what.

### Sample annotations

Next read in the sample annotations. Because samples are in columns in the count matrix, we will name the object `coldata`:


```r
coldata <- read.csv("data/GSE96870_coldata_cerebellum.csv",
                    row.names = 1)
dim(coldata)
```

```{.output}
[1] 22 10
```

```r
# View(coldata)
```

Now samples are in rows with the GEO sample IDs as the rownames, and we have 10 columns of information. The columns that are the most useful for this workshop are `geo_accession` (GEO sample IDs again), `sex` and `time`.

### Gene annotations
The counts only have gene symbols, which while short and somewhat recognizable to the human brain, are not always good absolute identifiers for exactly what gene was measured. For this we need additional gene annotations that were provided by the authors. The `count` and `coldata` files were in comma separated value (.csv) format, but we cannot use that for our gene annotation file because the descriptions can contain commas that would prevent a .csv file from being read in correctly. Instead the gene annotation file is in tab separated value (.tsv) format. Likewise, the descriptions can contain the single quote `'` (e.g., 5'), which by default R assumes indicates a character entry. So we have to use a more generic function `read.delim()` with extra arguments to specify that we have tab-separated data (`sep = "\t"`) with no quotes used (`quote = ""`). We also put in other arguments to specify that the first row contains our column names (`header = TRUE`), the gene symbols that should be our `row.names` are in the 5th column (`row.names = 5`), and that NCBI's species-specific gene ID (i.e., ENTREZID) should be read in as character data even though they look like numbers (`colClasses` argument). You can look up this details on available arguments by simply entering the function name starting with question mark. (e.g., `?read.delim`)


```r
rowranges <- read.delim("data/GSE96870_rowranges.tsv", 
                        sep = "\t", 
                        colClasses = c(ENTREZID = "character"),
                        header = TRUE, 
                        quote = "", 
                        row.names = 5)
dim(rowranges)
```

```{.output}
[1] 41786     7
```

```r
# View(rowranges)
```

For each of the 41,786 genes, we have the `seqnames` (e.g., chromosome number), 
`start` and `end` positions, `strand`, `ENTREZID`, gene product description 
(`product`) and the feature type (`gbkey`). These gene-level metadata are 
useful for the downstream analysis. For example, from the `gbkey` column, we
can check what types of genes and how many of them are in our dataset:


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

:::::::::::::::::::::::::::::::::::::::  challenge

## Challenge: Discuss the following points with your neighbor

- Have you used spreadsheets, in your research, courses,
  or at home?
- What kind of operations do you do in spreadsheets?
- Which ones do you think spreadsheets are good for?
- Have you accidentally done something in a spreadsheet program that made you
  frustrated or sad?
  
::::::::::::::::::::::::::::::::::::::::::::::::::

Suppose we decide we only want to analyze mRNA genes and not any of the others. We could remove these genes from the `rowranges` object, but we also have to remember to remove them correctly from the `counts` object. Likewise, we may decide one or more samples are outliers and we would have to remove them from both the `counts` columns and the `coldata` rows. You can see how this could easily lead to mis-matches between our counts and our annotations.

Instead, Bioconductor has created a specialized S4 class called a `SummarizedExperiment`. The details of a `SummarizedExperiment` were covered extensively at the end of the [Introduction to data analysis with R and Bioconductor](https://carpentries-incubator.github.io/bioc-intro/60-next-steps.html#next-steps) workshop. 
As a reminder, let's take a look at the figure below representing the anatomy of the `SummarizedExperiment` class:

<img src="https://uclouvain-cbio.github.io/WSBIM1322/figs/SE.svg" width="80%" style="display: block; margin: auto;" />

It is designed to hold any type of quantitative 'omics data (`assays`) along with linked sample annotations (`colData`) and feature annotations with (`rowRanges`) or without (`rowData`) chromosome, start and stop positions. Once these three tables are (correctly!) linked, subsetting either samples and/or features will correctly subset the `assay`, `colData` and `rowRanges`. Additionally, most Bioconductor packages are built around the same core data infrastructure so they will recognize and be able to manipulate `SummarizedExperiment` objects. Two of the most popular RNA-seq statistical analysis packages have their own extended S4 classes similar to a `SummarizedExperiment` with the additional slots for statistical results: [DESeq2](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#the-deseqdataset)
s `DESeqDataSet` and [edgeR](https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/DGEList-class)'s `DGEList`. No matter which one you end up using for statistical analysis, you can start by putting your data in a `SummarizedExperiment`. 

## Assemble SummarizedExperiment
We will create a `SummarizedExperiment` from these objects:

- The `count` object will be saved in `assays` slot    
- The `coldata` object with sample information will be stored in `colData` slot (_**sample metadata**_)    
- The `rowranges` object describing the genes will be stored in `rowRanges` slot (_**features metadata**_)     

Before we put them together, you ABSOLUTELY MUST MAKE SURE THE SAMPLES AND GENES ARE IN THE SAME ORDER! Even though we saw that `count` and `coldata` had the same number of samples and `count` and `rowranges` had the same number of genes, we never explicitly checked to see if they were in the same order. One way to check:



```r
all.equal(colnames(counts), rownames(coldata)) # samples
```

```{.output}
[1] TRUE
```

```r
all.equal(rownames(counts), rownames(rowranges)) # genes
```

```{.output}
[1] TRUE
```

```r
# If both not TRUE, you could do something like:

if(sum(colnames(counts) %in%  rownames(coldata)) == ncol(counts)) {
  tempindex <- match(colnames(counts), rownames(coldata))
  coldata <- coldata[tempindex, ]
} else {
  print("Warning: the number of samples are not the same in counts and coldata")
}
```

:::::::::::::::::::::::::::::::::::::::  challenge

If the feature (i.e., genes) in the assay (e.g., `counts`) and the gene
annotation table (e.g., `rowranges`) are different, how can we fix them? 
Write a code using `if-else` statement. 

::::::::::::::::::::::::::::::::::: solution


```r
if(sum(rownames(counts) %in%  rownames(rowranges)) == nrow(counts)) {
  tempindex <- match(rownames(counts), rownames(rowranges))
  coldata <- coldata[tempindex, ]
} else {
  print("Warning: the number of features are not the same in counts and rowranges")
}
```


:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::



Once we have verified that samples and genes are in the same order, we can 
then create our `SummarizedExperiment` object.


```r
# One final check:
stopifnot(rownames(rowranges) == rownames(counts), # features
          rownames(coldata) == colnames(counts)) # samples

se <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    rowRanges = as(rowranges, "GRanges"),
    colData = coldata
)
```


A brief recap of how to access the various data slots in a `SummarizedExperiment`:


```r
# Access the counts
head(assay(se))
```

```{.output}
             GSM2545336 GSM2545337 GSM2545338 GSM2545339 GSM2545340 GSM2545341
Xkr4               1891       2410       2159       1980       1977       1945
LOC105243853          0          0          1          4          0          0
LOC105242387        204        121        110        120        172        173
LOC105242467         12          5          5          5          2          6
Rp1                   2          2          0          3          2          1
Sox17               251        239        218        220        261        232
             GSM2545342 GSM2545343 GSM2545344 GSM2545345 GSM2545346 GSM2545347
Xkr4               1757       2235       1779       1528       1644       1585
LOC105243853          1          3          3          0          1          3
LOC105242387        177        130        131        160        180        176
LOC105242467          3          2          2          2          1          2
Rp1                   3          1          1          2          2          2
Sox17               179        296        233        271        205        230
             GSM2545348 GSM2545349 GSM2545350 GSM2545351 GSM2545352 GSM2545353
Xkr4               2275       1881       2584       1837       1890       1910
LOC105243853          1          0          0          1          1          0
LOC105242387        161        154        124        221        272        214
LOC105242467          2          4          7          1          3          1
Rp1                   3          6          5          3          5          1
Sox17               302        286        325        201        267        322
             GSM2545354 GSM2545362 GSM2545363 GSM2545380
Xkr4               1771       2315       1645       1723
LOC105243853          0          1          0          1
LOC105242387        124        189        223        251
LOC105242467          4          2          1          4
Rp1                   3          3          1          0
Sox17               273        197        310        246
```

```r
dim(assay(se))
```

```{.output}
[1] 41786    22
```

```r
# Access the sample annotations
colData(se)
```

```{.output}
DataFrame with 22 rows and 10 columns
                     title geo_accession     organism         age         sex
               <character>   <character>  <character> <character> <character>
GSM2545336 CNS_RNA-seq_10C    GSM2545336 Mus musculus     8 weeks      Female
GSM2545337 CNS_RNA-seq_11C    GSM2545337 Mus musculus     8 weeks      Female
GSM2545338 CNS_RNA-seq_12C    GSM2545338 Mus musculus     8 weeks      Female
GSM2545339 CNS_RNA-seq_13C    GSM2545339 Mus musculus     8 weeks      Female
GSM2545340 CNS_RNA-seq_14C    GSM2545340 Mus musculus     8 weeks        Male
...                    ...           ...          ...         ...         ...
GSM2545353  CNS_RNA-seq_3C    GSM2545353 Mus musculus     8 weeks      Female
GSM2545354  CNS_RNA-seq_4C    GSM2545354 Mus musculus     8 weeks        Male
GSM2545362  CNS_RNA-seq_5C    GSM2545362 Mus musculus     8 weeks      Female
GSM2545363  CNS_RNA-seq_6C    GSM2545363 Mus musculus     8 weeks        Male
GSM2545380  CNS_RNA-seq_9C    GSM2545380 Mus musculus     8 weeks      Female
             infection      strain        time      tissue     mouse
           <character> <character> <character> <character> <integer>
GSM2545336  InfluenzaA     C57BL/6        Day8  Cerebellum        14
GSM2545337 NonInfected     C57BL/6        Day0  Cerebellum         9
GSM2545338 NonInfected     C57BL/6        Day0  Cerebellum        10
GSM2545339  InfluenzaA     C57BL/6        Day4  Cerebellum        15
GSM2545340  InfluenzaA     C57BL/6        Day4  Cerebellum        18
...                ...         ...         ...         ...       ...
GSM2545353 NonInfected     C57BL/6        Day0  Cerebellum         4
GSM2545354 NonInfected     C57BL/6        Day0  Cerebellum         2
GSM2545362  InfluenzaA     C57BL/6        Day4  Cerebellum        20
GSM2545363  InfluenzaA     C57BL/6        Day4  Cerebellum        12
GSM2545380  InfluenzaA     C57BL/6        Day8  Cerebellum        19
```

```r
dim(colData(se))
```

```{.output}
[1] 22 10
```

```r
# Access the gene annotations
head(rowData(se))
```

```{.output}
DataFrame with 6 rows and 3 columns
                ENTREZID                product       gbkey
             <character>            <character> <character>
Xkr4              497097 X Kell blood group p..        mRNA
LOC105243853   105243853 uncharacterized LOC1..       ncRNA
LOC105242387   105242387 uncharacterized LOC1..       ncRNA
LOC105242467   105242467 lipoxygenase homolog..        mRNA
Rp1                19888 retinitis pigmentosa..        mRNA
Sox17              20671 SRY (sex determining..        mRNA
```

```r
dim(rowData(se))
```

```{.output}
[1] 41786     3
```



## Save SummarizedExperiment

This was a bit of code and time to create our `SummarizedExperiment` object. We will need to keep using it throughout the workshop, so it can be useful to save it as an actual single file on our computer to read it back in to R's memory if we have to shut down RStudio. To save an R-specific file we can use the `saveRDS()` function and later read it back into R using the `readRDS()` function. 


```r
saveRDS(se, "data/GSE96870_se.rds")
rm(se) # remove the object!
se <- readRDS("data/GSE96870_se.rds")
```


## Annotations
Depending on who generates your count data, you might not have a nice file of 
additional gene annotations. There may only be the count row names, which 
could be gene symbols or ENTREZIDs or another database's ID. Characteristics 
of gene annotations differ based on their annotation strategies and information 
sources. For example, RefSeq human gene models (i.e., Entrez from NCBI) are 
well supported and broadly used in various studies. The UCSC Known Genes 
dataset is based on protein data from Swiss-Prot/TrEMBL (UniProt) and the 
associated mRNA data from GenBank, and serves as a foundation for the UCSC 
Genome Browser. Ensembl genes contain both automated genome annotation and 
manual curation.

You can find more information in Bioconductor [Annotation Workshop](https://jmacdon.github.io/Bioc2022Anno/articles/AnnotationWorkshop.html)
material.

Bioconductor has many packages and functions that can help you to get additional annotation information for your genes. The available resources are covered in more detail in [Episode 6 Gene set enrichment analysis](https://carpentries-incubator.github.io/bioc-rnaseq/06-gene-set-analysis.html#gene-set-resources). 

Here, we will introduce one of the gene ID mapping functions, `mapIds`:
```
mapIds(annopkg, keys, column, keytype, ..., multiVals)
```

Where 
- *annopkg* is the annotation package        
- *keys* are the IDs that we **know**       
- *column* is the value we **want**    
- *keytype* is the type of key used    


```r
mapIds(org.Mm.eg.db, keys = "497097", column = "SYMBOL", keytype = "ENTREZID")
```

```{.output}
'select()' returned 1:1 mapping between keys and columns
```

```{.output}
497097 
"Xkr4" 
```

Different from the `select()` function, `mapIds()` function handles 1:many 
mapping between keys and columns through an additional argument, `multiVals`.
The below example demonstrate this functionality using the `hgu95av2.db` 
package, an Affymetrix Human Genome U95 Set annotation data.


```r
keys <- head(keys(hgu95av2.db, "ENTREZID"))
last <- function(x){x[[length(x)]]}

mapIds(hgu95av2.db, keys = keys, column = "ALIAS", keytype = "ENTREZID")
```

```{.output}
'select()' returned 1:many mapping between keys and columns
```

```{.output}
       10       100      1000     10000 100008586     10001 
   "AAC2"    "ADA1"   "ACOGS"    "MPPH"     "AL4"   "ARC33" 
```

```r
mapIds(hgu95av2.db, keys = keys, column = "ALIAS", keytype = "ENTREZID", multiVals = last)
```

```{.output}
'select()' returned 1:many mapping between keys and columns
```

```{.output}
       10       100      1000     10000 100008586     10001 
   "NAT2"     "ADA"    "CDH2"    "AKT3" "GAGE12F"    "MED6" 
```

```r
mapIds(hgu95av2.db, keys = keys, column = "ALIAS", keytype = "ENTREZID", multiVals = "list")
```

```{.output}
'select()' returned 1:many mapping between keys and columns
```

```{.output}
$`10`
[1] "AAC2"  "NAT-2" "PNAT"  "NAT2" 

$`100`
[1] "ADA1" "ADA" 

$`1000`
[1] "ACOGS"  "ADHD8"  "ARVD14" "CD325"  "CDHN"   "CDw325" "NCAD"   "CDH2"  

$`10000`
[1] "MPPH"         "MPPH2"        "PKB-GAMMA"    "PKBG"         "PRKBG"       
[6] "RAC-PK-gamma" "RAC-gamma"    "STK-2"        "AKT3"        

$`100008586`
[1] "AL4"     "CT4.7"   "GAGE-7"  "GAGE-7B" "GAGE-8"  "GAGE7"   "GAGE7B" 
[8] "GAGE12F"

$`10001`
[1] "ARC33"     "NY-REN-28" "MED6"     
```


:::::::::::::::::::::::::::::::::::::::::  callout

## Important!

Saving an object as an .rds file and reading it back in can be a convenient time saver for an extended analysis. However, there is no record in the object or .rds file of how it was made or manipulated. Your code to make the .rds file can become lost or corrupted, such that you cannot create the same object again from input data files. If this happens you should consider the .rds file as lost/corrupted as well. For reproducibility, the codes are more important than the actual .rds file.

::::::::::::::::::::::::::::::::::::::::::::::::::


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
 [1] hgu95av2.db_3.13.0          org.Hs.eg.db_3.17.0        
 [3] org.Mm.eg.db_3.17.0         AnnotationDbi_1.62.2       
 [5] SummarizedExperiment_1.30.2 Biobase_2.60.0             
 [7] MatrixGenerics_1.12.2       matrixStats_1.0.0          
 [9] GenomicRanges_1.52.0        GenomeInfoDb_1.36.1        
[11] IRanges_2.34.1              S4Vectors_0.38.1           
[13] BiocGenerics_0.46.0        

loaded via a namespace (and not attached):
 [1] Matrix_1.5-4.1          bit_4.0.5               highr_0.10             
 [4] compiler_4.3.1          BiocManager_1.30.21.1   renv_1.0.0             
 [7] crayon_1.5.2            blob_1.2.4              Biostrings_2.68.1      
[10] bitops_1.0-7            png_0.1-8               fastmap_1.1.1          
[13] lattice_0.21-8          R6_2.5.1                XVector_0.40.0         
[16] S4Arrays_1.0.5          knitr_1.43              DelayedArray_0.26.6    
[19] GenomeInfoDbData_1.2.10 DBI_1.1.3               rlang_1.1.1            
[22] KEGGREST_1.40.0         cachem_1.0.8            xfun_0.39              
[25] bit64_4.0.5             RSQLite_2.3.1           memoise_2.0.1          
[28] cli_3.6.1               zlibbioc_1.46.0         grid_4.3.1             
[31] vctrs_0.6.3             evaluate_0.21           abind_1.4-5            
[34] RCurl_1.98-1.12         httr_1.4.6              pkgconfig_2.0.3        
[37] tools_4.3.1            
```

::: keypoints
-   Depending on the gene expression quantification tool used, there are different ways (often distributed in Bioconductor packages) to read the output into a `SummarizedExperiment` or `DGEList` object for further processing in R.
-   Stable gene identifiers such as Ensembl or Entrez IDs should preferably be used as the main identifiers throughout an RNA-seq analysis, with gene symbols added for easier interpretation.
:::
