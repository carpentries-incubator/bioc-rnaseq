suppressPackageStartupMessages({
    library(GenomicRanges)
    library(SummarizedExperiment)
})

if (!file.exists("data/GSE96870_counts_cerebellum.csv")) {
    dir.create("data", showWarnings = FALSE)
    download.file(
        url = "https://github.com/Bioconductor/bioconductor-teaching/blob/master/data/GSE96870/GSE96870_counts_cerebellum.csv?raw=true", 
        destfile = "data/GSE96870_counts_cerebellum.csv"
    )
}

if (!file.exists("data/GSE96870_coldata_cerebellum.csv")) {
    dir.create("data", showWarnings = FALSE)
    download.file(
        url = "https://github.com/Bioconductor/bioconductor-teaching/blob/master/data/GSE96870/GSE96870_coldata_cerebellum.csv?raw=true", 
        destfile = "data/GSE96870_coldata_cerebellum.csv"
    )
}

if (!file.exists("data/GSE96870_coldata_all.csv")) {
    dir.create("data", showWarnings = FALSE)
    download.file(
        url = "https://github.com/Bioconductor/bioconductor-teaching/blob/master/data/GSE96870/GSE96870_coldata_all.csv?raw=true", 
        destfile = "data/GSE96870_coldata_all.csv"
    )
}

if (!file.exists("data/GSE96870_rowranges.tsv")) {
    dir.create("data", showWarnings = FALSE)
    download.file(
        url = "https://github.com/Bioconductor/bioconductor-teaching/blob/master/data/GSE96870/GSE96870_rowranges.tsv?raw=true", 
        destfile = "data/GSE96870_rowranges.tsv"
    )
}


if (!file.exists("data/GSE96870_se2.rds")) {
    counts_cerebellum <- read.csv("data/GSE96870_counts_cerebellum.csv",
                                  row.names = 1)
    coldata_cerebellum <- read.csv("data/GSE96870_coldata_cerebellum.csv",
                                   row.names = 1)
    rowranges <- read.delim("data/GSE96870_rowranges.tsv", sep = "\t", 
                            colClasses=c(ENTREZID="character"),
                            header = TRUE, quote = "", row.names = 5)

    stopifnot(rownames(rowranges) == rownames(counts_cerebellum),
              rownames(coldata_cerebellum) == colnames(counts_cerebellum))
    
    se <- SummarizedExperiment(
        assays = list(counts = as.matrix(counts_cerebellum)),
        rowRanges = as(rowranges, "GRanges"),
        colData = coldata_cerebellum
    )
    se$Label <- paste(se$sex, se$time, se$mouse, sep = "_")
    colnames(se) <- se$Label
    se$Group <- paste(se$sex, se$time, sep = "_")
    se$Group <- factor(se$Group, levels = c("Female_Day0","Male_Day0", 
                                            "Female_Day4","Male_Day4",
                                            "Female_Day8","Male_Day8"))
    se <- se[, order(se$Group)]
    se$Label <- factor(se$Label, levels = se$Label)
    
    saveRDS(se, file = "data/GSE96870_se.rds")
}

