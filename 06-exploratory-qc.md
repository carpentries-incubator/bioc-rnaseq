---
source: Rmd
title: Exploratory analysis and quality control
teaching: XX
exercises: XX
editor_options:
  chunk_output_type: console
---



::::::::::::::::::::::::::::::::::::::: objectives

- Learn how to explore the gene expression matrix and perform common quality control steps.
- Learn how to set up an interactive application for exploratory analysis.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- Why is exploratory analysis an essential part of an RNA-seq analysis? 
- How should one preprocess the raw count matrix for exploratory analysis?  
- Are two dimensions sufficient to represent your data?

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

### Contribute!

This episode is intended to introduce various types of exploratory analysis
and QC steps taken before a formal statistical analysis is done.


::::::::::::::::::::::::::::::::::::::::::::::::::


```r
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(DESeq2)
    library(vsn)
    library(ggplot2)
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(hexbin)
})
```


```r
se <- readRDS("data/GSE96870_se.rds")
```

Exploratory analysis is crucial for quality control and to get to know our data.
It can help us detect quality problems, sample swaps and contamination, as well as give us a sense of the most salient patterns present in the data.
In this episode, we will learn about two common ways of performing exploratory analysis for RNA-seq data; namely clustering and principal component analysis (PCA).
These tools are in no way limited to (or developed for) analysis of RNA-seq data.
However, there are certain characteristics of count assays that need to be taken into account when they are applied to this type of data.


```r
se <- se[rowSums(assay(se, "counts")) > 5, ]
dds <- DESeq2::DESeqDataSet(se[, se$tissue == "Cerebellum"],
                            design = ~ sex + time)
```

```{.warning}
Warning in DESeq2::DESeqDataSet(se[, se$tissue == "Cerebellum"], design = ~sex
+ : some variables in design formula are characters, converting to factors
```

## Library size differences


```r
ggplot(data.frame(sample = colnames(dds), 
                  libSize = colSums(assay(dds, "counts"))),
       aes(x = sample, y = libSize)) + 
    geom_bar(stat = "identity") + theme_bw() + 
    labs(x = "Sample", y = "Total count") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

<img src="fig/06-exploratory-qc-rendered-unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

Differences in the total number of reads assigned to genes between samples typically occur for technical reasons.
In practice, it means that we can not simply compare the raw read count directly between samples and conclude that a sample with a higher read count also expresses the gene more strongly - the higher count may be caused by an overall higher number of reads in that sample.
In the rest of this section, we will use the term *library size* to refer to the total number of reads assigned to genes for a sample.
We need to adjust for the differences in library size between samples, to avoid drawing incorrect conclusions.
The way this is typically done for RNA-seq data can be described as a two-step procedure.
First, we estimate *size factors* - sample-specific correction factors such that if the raw counts were to be divided by these factors, the resulting values would be more comparable across samples.
Next, these size factors are incorporated into the statistical analysis of the data.
It is important to pay close attention to how this is done in practice for a given analysis method.
Sometimes the division of the counts by the size factors needs to be done explicitly by the analyst.
Other times (as we will see for the differential expression analysis) it is important that they are provided separately to the analysis tool, which will then use them appropriately in the statistical model.

With `DESeq2`, size factors are calculated using the `estimateSizeFactors()` function.
The size factors estimated by this function combines an adjustment for differences in library sizes with an adjustment for differences in the RNA composition of the samples.
The latter is important due to the compositional nature of RNA-seq data.
There is a fixed number of reads to distribute between the genes, and if a single (or a few) very highly expressed gene consume a large part of the reads, all other genes will consequently receive very low counts.


```r
dds <- estimateSizeFactors(dds)
ggplot(data.frame(libSize = colSums(assay(dds, "counts")),
                  sizeFactor = sizeFactors(dds)),
       aes(x = libSize, y = sizeFactor)) + 
    geom_point(size = 5) + theme_bw() + 
    labs(x = "Library size", y = "Size factor")
```

<img src="fig/06-exploratory-qc-rendered-unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

## Transform data

There is a rich literature on methods for exploratory analysis.
Most of these work best in situations where the variance of the input data (here, each gene) is relatively independent of the average value.
For read count data such as RNA-seq, this is not the case.
In fact, the variance increases with the average read count.


```r
meanSdPlot(assay(dds), ranks = FALSE)
```

<img src="fig/06-exploratory-qc-rendered-unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

There are two ways around this: either we develop methods specifically adapted to count data, or we adapt (transform) the count data so that the existing methods are applicable.
Both ways have been explored; however, at the moment the second approach is arguably more widely applied in practice.


```r
vsd <- DESeq2::vst(dds, blind = TRUE)
meanSdPlot(assay(vsd), ranks = FALSE)
```

<img src="fig/06-exploratory-qc-rendered-unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

## Heatmaps and clustering


```r
dst <- dist(t(assay(vsd)))
colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)
ComplexHeatmap::Heatmap(
    as.matrix(dst), 
    col = colors,
    name = "Euclidean\ndistance",
    cluster_rows = hclust(dst),
    cluster_columns = hclust(dst),
    bottom_annotation = columnAnnotation(
        sex = vsd$sex,
        time = vsd$time,
        col = list(sex = c(Female = "red", Male = "blue"),
                   time = c(Day0 = "yellow", Day4 = "forestgreen", Day8 = "purple")))
)
```

<img src="fig/06-exploratory-qc-rendered-unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

## PCA

Principal component analysis is a dimensionality reduction method, which projects the samples into a lower-dimensional space.
This lower-dimensional representation can be used for visualization, or as the input for other analysis methods.
The principal components are defined in such a way that they are orthogonal, and that the projection of the samples into the space they span contains as much variance as possible.
It is an *unsupervised* method in the sense that no external information about the samples (e.g., the treatment condition) is taken into account.
In the plot below we represent the samples in a two-dimensional principal component space.
For each of the two dimensions, we indicate the fraction of the total variance that is represented by that component.
By definition, the first principal component will always represent more of the variance than the subsequent ones.
The fraction of explained variance is a measure of how much of the 'signal' in the data that is retained when we project the samples from the original, high-dimensional space to the low-dimensional space for visualization.


```r
pcaData <- DESeq2::plotPCA(vsd, intgroup = c("sex", "time"),
                           returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sex, shape = time), size = 5) +
    theme_minimal() +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() + 
    scale_color_manual(values = c(Male = "blue", Female = "red"))
```

<img src="fig/06-exploratory-qc-rendered-unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

## Session info


```r
sessionInfo()
```

```{.output}
R version 4.2.3 (2023-03-15)
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

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] hexbin_1.28.3               RColorBrewer_1.1-3         
 [3] ComplexHeatmap_2.14.0       ggplot2_3.4.1              
 [5] vsn_3.66.0                  DESeq2_1.38.3              
 [7] SummarizedExperiment_1.28.0 Biobase_2.58.0             
 [9] MatrixGenerics_1.10.0       matrixStats_0.63.0         
[11] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
[13] IRanges_2.32.0              S4Vectors_0.36.2           
[15] BiocGenerics_0.44.0        

loaded via a namespace (and not attached):
 [1] httr_1.4.5             bit64_4.0.5            foreach_1.5.2         
 [4] highr_0.10             BiocManager_1.30.20    affy_1.76.0           
 [7] blob_1.2.4             renv_0.17.3            GenomeInfoDbData_1.2.9
[10] pillar_1.9.0           RSQLite_2.3.1          lattice_0.20-45       
[13] glue_1.6.2             limma_3.54.2           digest_0.6.31         
[16] XVector_0.38.0         colorspace_2.1-0       preprocessCore_1.60.2 
[19] Matrix_1.5-3           XML_3.99-0.14          pkgconfig_2.0.3       
[22] GetoptLong_1.0.5       zlibbioc_1.44.0        xtable_1.8-4          
[25] scales_1.2.1           affyio_1.68.0          BiocParallel_1.32.6   
[28] tibble_3.2.1           annotate_1.76.0        KEGGREST_1.38.0       
[31] farver_2.1.1           generics_0.1.3         cachem_1.0.7          
[34] withr_2.5.0            cli_3.6.1              magrittr_2.0.3        
[37] crayon_1.5.2           memoise_2.0.1          evaluate_0.20         
[40] fansi_1.0.4            doParallel_1.0.17      tools_4.2.3           
[43] GlobalOptions_0.1.2    lifecycle_1.0.3        munsell_0.5.0         
[46] locfit_1.5-9.7         cluster_2.1.4          DelayedArray_0.24.0   
[49] AnnotationDbi_1.60.2   Biostrings_2.66.0      compiler_4.2.3        
[52] rlang_1.1.0            RCurl_1.98-1.12        iterators_1.0.14      
[55] circlize_0.4.15        rjson_0.2.21           labeling_0.4.2        
[58] bitops_1.0-7           gtable_0.3.3           codetools_0.2-19      
[61] DBI_1.1.3              R6_2.5.1               knitr_1.42            
[64] dplyr_1.1.1            fastmap_1.1.1          bit_4.0.5             
[67] utf8_1.2.3             clue_0.3-64            shape_1.4.6           
[70] parallel_4.2.3         Rcpp_1.0.10            vctrs_0.6.1           
[73] geneplotter_1.76.0     png_0.1-8              tidyselect_1.2.0      
[76] xfun_0.38             
```

:::::::::::::::::::::::::::::::::::::::: keypoints

- Exploratory analysis is essential for quality control and to detect potential problems with a data set. 
- Different classes of exploratory analysis methods expect differently preprocessed data. The most commonly used methods expect counts to be normalized and log-transformed (or similar- more sensitive/sophisticated), to be closer to homoskedastic. Other methods work directly on the raw counts.  


::::::::::::::::::::::::::::::::::::::::::::::::::


