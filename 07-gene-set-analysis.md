---
source: Rmd
title: Gene set analysis
teaching: XX
exercises: XX
---



::::::::::::::::::::::::::::::::::::::: objectives

- Explain how to find differentially expressed pathways with gene set analysis in R.
- Understand how differentially expressed genes can enrich a gene set.
- Explain how to perform a gene set analysis in R, using clusterProfiler.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do we find differentially expressed pathways?

::::::::::::::::::::::::::::::::::::::::::::::::::



:::::::::::::::::::::::::::::::::::::::::  callout

### Contribute!

This episode is intended to introduce the concept of how to carry out a
functional analysis of a subset of differentially expressed (DE) genes,
by means of assessing how significantly DE genes enrich gene sets of
our interest.


::::::::::::::::::::::::::::::::::::::::::::::::::

First, we are going to explore the basic concept of enriching a gene set with
differentially expressed (DE) genes. Recall the differential expression analysis.


```r
library(SummarizedExperiment)
library(DESeq2)
```


```r
se <- readRDS("data/GSE96870_se.rds")
```


```r
dds <- DESeq2::DESeqDataSet(se[, se$tissue == "Cerebellum"],
                            design = ~ sex + time)
```

```{.warning}
Warning in DESeq2::DESeqDataSet(se[, se$tissue == "Cerebellum"], design = ~sex
+ : some variables in design formula are characters, converting to factors
```


```r
dds <- DESeq2::DESeq(dds)
```

Fetch results for the contrast between male and female mice.


```r
resSex <- DESeq2::results(dds, contrast = c("sex", "Male", "Female"))
summary(resSex)
```

```{.output}

out of 32652 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 53, 0.16%
LFC < 0 (down)     : 71, 0.22%
outliers [1]       : 10, 0.031%
low counts [2]     : 13717, 42%
(mean count < 6)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

Select DE genes between males and females with FDR \< 5%.


```r
sexDE <- as.data.frame(subset(resSex, padj < 0.05))
dim(sexDE)
```

```{.output}
[1] 54  6
```

```r
sexDE <- sexDE[order(abs(sexDE$log2FoldChange), decreasing=TRUE), ]
head(sexDE)
```

```{.output}
               baseMean log2FoldChange     lfcSE      stat        pvalue
Eif2s3y       1410.8750       12.62514 0.5652155  22.33685 1.620659e-110
Kdm5d          692.1672       12.55386 0.5936267  21.14773  2.895664e-99
Uty            667.4375       12.01728 0.5935911  20.24505  3.927797e-91
Ddx3y         2072.9436       11.87241 0.3974927  29.86825 5.087220e-196
Xist         22603.0359      -11.60429 0.3362822 -34.50761 6.168523e-261
LOC105243748    52.9669        9.08325 0.5976242  15.19893  3.594320e-52
                      padj
Eif2s3y      1.022366e-106
Kdm5d         1.370011e-95
Uty           1.486671e-87
Ddx3y        4.813782e-192
Xist         1.167393e-256
LOC105243748  1.133708e-48
```

```r
sexDEgenes <- rownames(sexDE)
head(sexDEgenes)
```

```{.output}
[1] "Eif2s3y"      "Kdm5d"        "Uty"          "Ddx3y"        "Xist"        
[6] "LOC105243748"
```

```r
length(sexDEgenes)
```

```{.output}
[1] 54
```

## Enrichment of a curated gene set

:::::::::::::::::::::::::::::::::::::::::  callout

### Contribute!

Here we illustrate how to assess the enrichment of one gene set
we curate ourselves with our subset of DE genes with sex-specific
expression. Here we form such a gene set with genes from sex
chromosomes. Could you think of another more accurate gene set formed
by genes with sex-specific expression?


::::::::::::::::::::::::::::::::::::::::::::::::::

Build a gene set formed by genes located in the sex chromosomes X and Y.


```r
xygenes <- rownames(se)[decode(seqnames(rowRanges(se)) %in% c("X", "Y"))]
length(xygenes)
```

```{.output}
[1] 2123
```

Build a contingency table and conduct a one-tailed Fisher's exact test that
verifies the association between genes being DE between male and female mice
and being located in a sex chromosome.


```r
N <- nrow(se)
n <- length(sexDEgenes)
m <- length(xygenes)
k <- length(intersect(xygenes, sexDEgenes)) 
dnames <- list(GS=c("inside", "outside"), DE=c("yes", "no"))
t <- matrix(c(k, n-k, m-k, N+k-n-m),
                        nrow=2, ncol=2, dimnames=dnames)
t
```

```{.output}
         DE
GS        yes    no
  inside   18  2105
  outside  36 39627
```

```r
fisher.test(t, alternative="greater")
```

```{.output}

	Fisher's Exact Test for Count Data

data:  t
p-value = 7.944e-11
alternative hypothesis: true odds ratio is greater than 1
95 percent confidence interval:
 5.541517      Inf
sample estimates:
odds ratio 
  9.411737 
```

## Gene ontology analysis with clusterProfiler

:::::::::::::::::::::::::::::::::::::::::  callout

### Contribute!

Here we illustrate how to assess the enrichment on the entire
collection of Gene Ontology (GO) gene sets using the package
clusterProfiler. Could we illustrate any missing important feature
of this package for this analysis objective? Could we briefly
mention other packages that may be useful for this task?


::::::::::::::::::::::::::::::::::::::::::::::::::

Second, let's perform a gene set analysis for an entire collection of gene sets
using the Bioconductor package
[clusterProfiler](https://bioconductor.org/packages/clusterProfiler). For this
purpose, we will fetch the results for the contrast between two time points.


```r
resTime <- DESeq2::results(dds, contrast = c("time", "Day8", "Day0"))
summary(resTime)
```

```{.output}

out of 32652 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 4472, 14%
LFC < 0 (down)     : 4276, 13%
outliers [1]       : 10, 0.031%
low counts [2]     : 8732, 27%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

Select DE genes between `Day0` and `Day8` with FDR \< 5% and minimum 1.5-fold
change.


```r
timeDE <- as.data.frame(subset(resTime, padj < 0.05 & abs(log2FoldChange) > log2(1.5)))
dim(timeDE)
```

```{.output}
[1] 2110    6
```

```r
timeDE <- timeDE[order(abs(timeDE$log2FoldChange), decreasing=TRUE), ]
head(timeDE)
```

```{.output}
              baseMean log2FoldChange     lfcSE      stat       pvalue
LOC105245444  2.441873       4.768938 0.9013067  5.291138 1.215573e-07
LOC105246405  9.728219       4.601505 0.6101832  7.541186 4.657174e-14
4933427D06Rik 1.480365       4.556126 1.0318402  4.415535 1.007607e-05
A930006I01Rik 2.312732      -4.353155 0.9176026 -4.744053 2.094837e-06
LOC105245223  3.272536       4.337202 0.8611255  5.036666 4.737099e-07
A530053G22Rik 1.554735       4.243903 1.0248977  4.140806 3.460875e-05
                      padj
LOC105245444  1.800765e-06
LOC105246405  2.507951e-12
4933427D06Rik 9.169093e-05
A930006I01Rik 2.252139e-05
LOC105245223  6.047199e-06
A530053G22Rik 2.720142e-04
```

```r
timeDEgenes <- rownames(timeDE)
head(timeDEgenes)
```

```{.output}
[1] "LOC105245444"  "LOC105246405"  "4933427D06Rik" "A930006I01Rik"
[5] "LOC105245223"  "A530053G22Rik"
```

```r
length(timeDEgenes)
```

```{.output}
[1] 2110
```

Call the `enrichGO()` function from
[clusterProfiler](https://bioconductor.org/packages/clusterProfiler)
as follows.


```r
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

resTimeGO <- enrichGO(gene = timeDEgenes,
                      keyType = "SYMBOL",
                      universe = rownames(se),
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01)
dim(resTimeGO)
```

```{.output}
[1] 17  9
```

```r
head(resTimeGO)
```

```{.output}
                   ID                                           Description
GO:0035456 GO:0035456                           response to interferon-beta
GO:0071674 GO:0071674                            mononuclear cell migration
GO:0035458 GO:0035458                  cellular response to interferon-beta
GO:0050900 GO:0050900                                   leukocyte migration
GO:0030595 GO:0030595                                  leukocyte chemotaxis
GO:0002523 GO:0002523 leukocyte migration involved in inflammatory response
           GeneRatio   BgRatio       pvalue     p.adjust       qvalue
GO:0035456   17/1260  54/21417 6.434058e-09 3.254347e-05 3.082930e-05
GO:0071674   32/1260 179/21417 1.596474e-08 4.037482e-05 3.824814e-05
GO:0035458   15/1260  47/21417 4.064983e-08 6.749342e-05 6.393832e-05
GO:0050900   50/1260 375/21417 5.337558e-08 6.749342e-05 6.393832e-05
GO:0030595   34/1260 222/21417 2.880367e-07 2.913780e-04 2.760302e-04
GO:0002523   10/1260  25/21417 6.944905e-07 5.854555e-04 5.546176e-04
                                                                                                                                                                                                                                                                                             geneID
GO:0035456                                                                                                                                                                            Tgtp1/Tgtp2/F830016B08Rik/Iigp1/Ifitm6/Igtp/Gm4951/Bst2/Oas1c/Irgm1/Gbp6/Ifi47/Aim2/Ifitm7/Irgm2/Ifit1/Ifi204
GO:0071674                                                                                                    Tnfsf18/Aire/Ccl17/Ccr7/Nlrp12/Ccl2/Retnlg/Apod/Il12a/Ccl5/Fut7/Ccl7/Spn/Itgb3/Grem1/Ptk2b/Lgals3/Adam8/Dusp1/Ch25h/Nbl1/Alox5/Padi2/Plg/Calr/Ager/Slamf9/Ccl6/Mdk/Itga4/Hsd3b7/Trpm4
GO:0035458                                                                                                                                                                                        Tgtp1/Tgtp2/F830016B08Rik/Iigp1/Igtp/Gm4951/Oas1c/Irgm1/Gbp6/Ifi47/Aim2/Ifitm7/Irgm2/Ifit1/Ifi204
GO:0050900 Tnfsf18/Aire/Ccl17/Ccr7/Nlrp12/Bst1/Ccl2/Retnlg/Ppbp/Cxcl5/Apod/Il12a/Ccl5/Umodl1/Fut7/Ccl7/Ccl28/Spn/Sell/Itgb3/Grem1/Cxcl1/Ptk2b/Lgals3/Adam8/Pf4/Dusp1/Ch25h/S100a8/Nbl1/Alox5/Padi2/Plg/Edn3/Il33/Ptn/Ada/Enpp1/Calr/Ager/Slamf9/Ccl6/Prex1/Aoc3/Itgam/Mdk/Itga4/Hsd3b7/P2ry12/Trpm4
GO:0030595                                                                                        Tnfsf18/Ccl17/Ccr7/Bst1/Ccl2/Retnlg/Ppbp/Cxcl5/Il12a/Ccl5/Ccl7/Sell/Grem1/Cxcl1/Ptk2b/Lgals3/Adam8/Pf4/Dusp1/Ch25h/S100a8/Nbl1/Alox5/Padi2/Edn3/Ptn/Calr/Slamf9/Ccl6/Prex1/Itgam/Mdk/Hsd3b7/Trpm4
GO:0002523                                                                                                                                                                                                                                     Ccl2/Ppbp/Fut7/Adam8/S100a8/Alox5/Ptn/Aoc3/Itgam/Mdk
           Count
GO:0035456    17
GO:0071674    32
GO:0035458    15
GO:0050900    50
GO:0030595    34
GO:0002523    10
```

Let's build a more readable table of results.


```r
library(kableExtra)

resTimeGOtab <- as.data.frame(resTimeGO)
resTimeGOtab$ID <- NULL
resTimeGOtab$geneID <- sapply(strsplit(resTimeGO$geneID, "/"), paste, collapse=", ")
ktab <- kable(resTimeGOtab, row.names=TRUE, caption="GO results for DE genes between time points.")
kable_styling(ktab, bootstrap_options=c("stripped", "hover", "responsive"), fixed_thead=TRUE)
```

<table class="table table-hover table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>GO results for DE genes between time points.</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Description </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> GeneRatio </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> BgRatio </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> pvalue </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adjust </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> qvalue </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> geneID </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Count </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> GO:0035456 </td>
   <td style="text-align:left;"> response to interferon-beta </td>
   <td style="text-align:left;"> 17/1260 </td>
   <td style="text-align:left;"> 54/21417 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 0.0000325 </td>
   <td style="text-align:right;"> 0.0000308 </td>
   <td style="text-align:left;"> Tgtp1, Tgtp2, F830016B08Rik, Iigp1, Ifitm6, Igtp, Gm4951, Bst2, Oas1c, Irgm1, Gbp6, Ifi47, Aim2, Ifitm7, Irgm2, Ifit1, Ifi204 </td>
   <td style="text-align:right;"> 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0071674 </td>
   <td style="text-align:left;"> mononuclear cell migration </td>
   <td style="text-align:left;"> 32/1260 </td>
   <td style="text-align:left;"> 179/21417 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 0.0000404 </td>
   <td style="text-align:right;"> 0.0000382 </td>
   <td style="text-align:left;"> Tnfsf18, Aire, Ccl17, Ccr7, Nlrp12, Ccl2, Retnlg, Apod, Il12a, Ccl5, Fut7, Ccl7, Spn, Itgb3, Grem1, Ptk2b, Lgals3, Adam8, Dusp1, Ch25h, Nbl1, Alox5, Padi2, Plg, Calr, Ager, Slamf9, Ccl6, Mdk, Itga4, Hsd3b7, Trpm4 </td>
   <td style="text-align:right;"> 32 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0035458 </td>
   <td style="text-align:left;"> cellular response to interferon-beta </td>
   <td style="text-align:left;"> 15/1260 </td>
   <td style="text-align:left;"> 47/21417 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 0.0000675 </td>
   <td style="text-align:right;"> 0.0000639 </td>
   <td style="text-align:left;"> Tgtp1, Tgtp2, F830016B08Rik, Iigp1, Igtp, Gm4951, Oas1c, Irgm1, Gbp6, Ifi47, Aim2, Ifitm7, Irgm2, Ifit1, Ifi204 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0050900 </td>
   <td style="text-align:left;"> leukocyte migration </td>
   <td style="text-align:left;"> 50/1260 </td>
   <td style="text-align:left;"> 375/21417 </td>
   <td style="text-align:right;"> 1.00e-07 </td>
   <td style="text-align:right;"> 0.0000675 </td>
   <td style="text-align:right;"> 0.0000639 </td>
   <td style="text-align:left;"> Tnfsf18, Aire, Ccl17, Ccr7, Nlrp12, Bst1, Ccl2, Retnlg, Ppbp, Cxcl5, Apod, Il12a, Ccl5, Umodl1, Fut7, Ccl7, Ccl28, Spn, Sell, Itgb3, Grem1, Cxcl1, Ptk2b, Lgals3, Adam8, Pf4, Dusp1, Ch25h, S100a8, Nbl1, Alox5, Padi2, Plg, Edn3, Il33, Ptn, Ada, Enpp1, Calr, Ager, Slamf9, Ccl6, Prex1, Aoc3, Itgam, Mdk, Itga4, Hsd3b7, P2ry12, Trpm4 </td>
   <td style="text-align:right;"> 50 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0030595 </td>
   <td style="text-align:left;"> leukocyte chemotaxis </td>
   <td style="text-align:left;"> 34/1260 </td>
   <td style="text-align:left;"> 222/21417 </td>
   <td style="text-align:right;"> 3.00e-07 </td>
   <td style="text-align:right;"> 0.0002914 </td>
   <td style="text-align:right;"> 0.0002760 </td>
   <td style="text-align:left;"> Tnfsf18, Ccl17, Ccr7, Bst1, Ccl2, Retnlg, Ppbp, Cxcl5, Il12a, Ccl5, Ccl7, Sell, Grem1, Cxcl1, Ptk2b, Lgals3, Adam8, Pf4, Dusp1, Ch25h, S100a8, Nbl1, Alox5, Padi2, Edn3, Ptn, Calr, Slamf9, Ccl6, Prex1, Itgam, Mdk, Hsd3b7, Trpm4 </td>
   <td style="text-align:right;"> 34 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0002523 </td>
   <td style="text-align:left;"> leukocyte migration involved in inflammatory response </td>
   <td style="text-align:left;"> 10/1260 </td>
   <td style="text-align:left;"> 25/21417 </td>
   <td style="text-align:right;"> 7.00e-07 </td>
   <td style="text-align:right;"> 0.0005855 </td>
   <td style="text-align:right;"> 0.0005546 </td>
   <td style="text-align:left;"> Ccl2, Ppbp, Fut7, Adam8, S100a8, Alox5, Ptn, Aoc3, Itgam, Mdk </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0050953 </td>
   <td style="text-align:left;"> sensory perception of light stimulus </td>
   <td style="text-align:left;"> 26/1260 </td>
   <td style="text-align:left;"> 162/21417 </td>
   <td style="text-align:right;"> 2.90e-06 </td>
   <td style="text-align:right;"> 0.0020909 </td>
   <td style="text-align:right;"> 0.0019808 </td>
   <td style="text-align:left;"> Aipl1, Vsx2, Nxnl2, Lrit3, Cryba2, Bfsp2, Lrat, Gabrr2, Lum, Rlbp1, Pde6g, Gpr179, Col1a1, Cplx3, Best1, Ush1g, Rs1, Rdh5, Guca1b, Th, Ppef2, Rbp4, Olfm2, Rom1, Vsx1, Rpe65 </td>
   <td style="text-align:right;"> 26 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0071675 </td>
   <td style="text-align:left;"> regulation of mononuclear cell migration </td>
   <td style="text-align:left;"> 21/1260 </td>
   <td style="text-align:left;"> 117/21417 </td>
   <td style="text-align:right;"> 4.20e-06 </td>
   <td style="text-align:right;"> 0.0026851 </td>
   <td style="text-align:right;"> 0.0025436 </td>
   <td style="text-align:left;"> Tnfsf18, Aire, Ccr7, Ccl2, Apod, Il12a, Ccl5, Ccl7, Spn, Itgb3, Grem1, Ptk2b, Lgals3, Adam8, Dusp1, Nbl1, Padi2, Calr, Ager, Mdk, Itga4 </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0007601 </td>
   <td style="text-align:left;"> visual perception </td>
   <td style="text-align:left;"> 25/1260 </td>
   <td style="text-align:left;"> 158/21417 </td>
   <td style="text-align:right;"> 5.80e-06 </td>
   <td style="text-align:right;"> 0.0032508 </td>
   <td style="text-align:right;"> 0.0030795 </td>
   <td style="text-align:left;"> Aipl1, Vsx2, Nxnl2, Lrit3, Cryba2, Bfsp2, Lrat, Gabrr2, Lum, Rlbp1, Pde6g, Gpr179, Col1a1, Cplx3, Best1, Rs1, Rdh5, Guca1b, Th, Ppef2, Rbp4, Olfm2, Rom1, Vsx1, Rpe65 </td>
   <td style="text-align:right;"> 25 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0034341 </td>
   <td style="text-align:left;"> response to interferon-gamma </td>
   <td style="text-align:left;"> 22/1260 </td>
   <td style="text-align:left;"> 135/21417 </td>
   <td style="text-align:right;"> 1.28e-05 </td>
   <td style="text-align:right;"> 0.0056936 </td>
   <td style="text-align:right;"> 0.0053937 </td>
   <td style="text-align:left;"> Ccl17, Gbp4, Ccl2, Tgtp1, H2-Q7, Il12rb1, Ifitm6, Ccl5, Ccl7, Igtp, Nos2, Nlrc5, Bst2, Irgm1, Gbp6, Capg, Ifitm7, Gbp9, Gbp5, Irgm2, Ccl6, Aqp4 </td>
   <td style="text-align:right;"> 22 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0060326 </td>
   <td style="text-align:left;"> cell chemotaxis </td>
   <td style="text-align:left;"> 38/1260 </td>
   <td style="text-align:left;"> 308/21417 </td>
   <td style="text-align:right;"> 1.30e-05 </td>
   <td style="text-align:right;"> 0.0056936 </td>
   <td style="text-align:right;"> 0.0053937 </td>
   <td style="text-align:left;"> Tnfsf18, Ccl17, Ccr7, Bst1, Ccl2, Retnlg, Ppbp, Cxcl5, Nr4a1, Il12a, Ccl5, Ccl7, Ccl28, Sell, Grem1, Cxcl1, Ptk2b, Lgals3, Adam8, Pf4, Dusp1, Ch25h, S100a8, Nbl1, Alox5, Padi2, Edn3, Ptn, Plxnb3, Calr, Lpar1, Slamf9, Ccl6, Prex1, Itgam, Mdk, Hsd3b7, Trpm4 </td>
   <td style="text-align:right;"> 38 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0002685 </td>
   <td style="text-align:left;"> regulation of leukocyte migration </td>
   <td style="text-align:left;"> 31/1260 </td>
   <td style="text-align:left;"> 230/21417 </td>
   <td style="text-align:right;"> 1.41e-05 </td>
   <td style="text-align:right;"> 0.0056936 </td>
   <td style="text-align:right;"> 0.0053937 </td>
   <td style="text-align:left;"> Tnfsf18, Aire, Ccr7, Bst1, Ccl2, Apod, Il12a, Ccl5, Fut7, Ccl7, Ccl28, Spn, Sell, Itgb3, Grem1, Ptk2b, Lgals3, Adam8, Dusp1, Nbl1, Padi2, Edn3, Il33, Ptn, Ada, Calr, Ager, Aoc3, Mdk, Itga4, P2ry12 </td>
   <td style="text-align:right;"> 31 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0036336 </td>
   <td style="text-align:left;"> dendritic cell migration </td>
   <td style="text-align:left;"> 9/1260 </td>
   <td style="text-align:left;"> 27/21417 </td>
   <td style="text-align:right;"> 1.46e-05 </td>
   <td style="text-align:right;"> 0.0056936 </td>
   <td style="text-align:right;"> 0.0053937 </td>
   <td style="text-align:left;"> Tnfsf18, Ccr7, Nlrp12, Retnlg, Il12a, Alox5, Calr, Slamf9, Trpm4 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:1990266 </td>
   <td style="text-align:left;"> neutrophil migration </td>
   <td style="text-align:left;"> 21/1260 </td>
   <td style="text-align:left;"> 127/21417 </td>
   <td style="text-align:right;"> 1.59e-05 </td>
   <td style="text-align:right;"> 0.0057426 </td>
   <td style="text-align:right;"> 0.0054401 </td>
   <td style="text-align:left;"> Ccl17, Ccr7, Bst1, Ccl2, Ppbp, Cxcl5, Ccl5, Umodl1, Fut7, Ccl7, Sell, Cxcl1, Lgals3, Adam8, Pf4, S100a8, Edn3, Ccl6, Prex1, Itgam, Mdk </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0016264 </td>
   <td style="text-align:left;"> gap junction assembly </td>
   <td style="text-align:left;"> 7/1260 </td>
   <td style="text-align:left;"> 16/21417 </td>
   <td style="text-align:right;"> 1.71e-05 </td>
   <td style="text-align:right;"> 0.0057811 </td>
   <td style="text-align:right;"> 0.0054766 </td>
   <td style="text-align:left;"> Gjd3, Aplnr, Gja5, Agt, Hopx, Cav1, Ace </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0030593 </td>
   <td style="text-align:left;"> neutrophil chemotaxis </td>
   <td style="text-align:left;"> 18/1260 </td>
   <td style="text-align:left;"> 100/21417 </td>
   <td style="text-align:right;"> 1.94e-05 </td>
   <td style="text-align:right;"> 0.0061338 </td>
   <td style="text-align:right;"> 0.0058107 </td>
   <td style="text-align:left;"> Ccl17, Ccr7, Bst1, Ccl2, Ppbp, Cxcl5, Ccl5, Ccl7, Sell, Cxcl1, Lgals3, Pf4, S100a8, Edn3, Ccl6, Prex1, Itgam, Mdk </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:1903028 </td>
   <td style="text-align:left;"> positive regulation of opsonization </td>
   <td style="text-align:left;"> 6/1260 </td>
   <td style="text-align:left;"> 12/21417 </td>
   <td style="text-align:right;"> 2.78e-05 </td>
   <td style="text-align:right;"> 0.0082825 </td>
   <td style="text-align:right;"> 0.0078462 </td>
   <td style="text-align:left;"> Pla2g5, Masp2, C4b, Colec11, C3, Cfp </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
</tbody>
</table>

:::::::::::::::::::::::::::::::::::::::: keypoints

- Key point 1

::::::::::::::::::::::::::::::::::::::::::::::::::


