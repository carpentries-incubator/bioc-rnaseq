---
title: Setup
---

Check you have the most recent versions of R and RStudio. See the instructions under "If you already have R and RStudio installed" in [Introduction to R](https://carpentries-incubator.github.io/bioc-intro/#r-and-rstudio). You will also need to install the packages below that will be used during the workshop. 

```r
install.packages("BiocManager")
BiocManager::install(c("tidyverse", "SummarizedExperiment",
                       "ExploreModelMatrix", "org.Mm.eg.db",
                       "DESeq2", "vsn", "ComplexHeatmap",
                       "RColorBrewer", "hexbin", "cowplot",
                       "clusterProfiler", "enrichplot", "kableExtra"))
```








