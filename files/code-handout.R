#' ## Introduction to RNA-seq
#' 
#' ## Challenge: Discuss the following points with your neighbor
#' 
#' 1. Which of the mentioned RNA-Seq quantification tools have you heard about? Do you know other pros and cons of the methods?
#' 2. Have you done your own RNA-Seq experiment? If so, what quantification tool did you use and why did you choose it?
#' 3. Do you have access to specific tools / local bioinformatics expert / computational resources for quantification? If you don't, how might you gain access?
#' 
#' ## RStudio Project and Experimental Data
#' 
#' ## Create the directories for subsequent episodes
#' 
#' Create a directory on your computer to serve as the working directory for the rest
#' of this episode and lesson (the workshop example uses a directory called `bio_rnaseq`).
#' Then, within this chosen directory, create the four fundamental directories previously discussed
#' (`data`, `scripts`, `documents`, and `output`).
#' 
#' ## Download the remaining data set files
#' 
#' There are three more data set files we need to download for the remainder of this lesson.
#' 
#' | URL | Filename                        | 
#' | --- | ------------------------------- |
#' | [https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870\_coldata\_cerebellum.csv](https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870_coldata_cerebellum.csv)    | GSE96870\_coldata\_cerebellum.csv | 
#' | [https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870\_coldata\_all.csv](https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870_coldata_all.csv)    | GSE96870\_coldata\_all.csv        | 
#' | [https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870\_rowranges.tsv](https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870_rowranges.tsv)    | GSE96870\_rowranges.tsv          | 
#' 
#' Use the `download.file` function to download the files into the `data` folder in your
#' working directory.
#' 
#' ## Importing and annotating quantified data into R
#' 
## ---- purl=TRUE---------------------------------------------------------------
suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.Mm.eg.db)
    library(hgu95av2.db)
    library(SummarizedExperiment)
})

#' 
## ---- purl=TRUE---------------------------------------------------------------
counts <- read.csv("data/GSE96870_counts_cerebellum.csv", 
                   row.names = 1)
dim(counts)
# View(counts)

#' 
## ---- purl=TRUE---------------------------------------------------------------
coldata <- read.csv("data/GSE96870_coldata_cerebellum.csv",
                    row.names = 1)
dim(coldata)
# View(coldata)

#' 
## ---- purl=TRUE---------------------------------------------------------------
rowranges <- read.delim("data/GSE96870_rowranges.tsv", 
                        sep = "\t", 
                        colClasses = c(ENTREZID = "character"),
                        header = TRUE, 
                        quote = "", 
                        row.names = 5)
dim(rowranges)
# View(rowranges)

#' 
## ---- purl=TRUE---------------------------------------------------------------
table(rowranges$gbkey)

#' 
#' ## Challenge: Discuss the following points with your neighbor
#' 
#' 1. How are the 3 objects `counts`, `coldata` and `rowranges` related to each other in terms of their rows and columns?
#' 2. If you only wanted to analyse the mRNA genes, what would you have to do keep just those (generally speaking, not exact codes)?
#' 3. If you decided the first two samples were outliers, what would you have to do to remove those (generally speaking, not exact codes)?
#' 
## ---- purl=TRUE---------------------------------------------------------------
all.equal(colnames(counts), rownames(coldata)) # samples
all.equal(rownames(counts), rownames(rowranges)) # genes

# If the first is not TRUE, you can match up the samples/columns in
# counts with the samples/rows in coldata like this (which is fine
# to run even if the first was TRUE):

tempindex <- match(colnames(counts), rownames(coldata))
coldata <- coldata[tempindex, ]

# Check again:
all.equal(colnames(counts), rownames(coldata)) 


#' 
#' If the features (i.e., genes) in the assay (e.g., `counts`) and the gene
#' annotation table (e.g., `rowranges`) are different, how can we fix them?
#' Write the codes.
#' 
## ---- purl=TRUE---------------------------------------------------------------
# One final check:
stopifnot(rownames(rowranges) == rownames(counts), # features
          rownames(coldata) == colnames(counts)) # samples

se <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    rowRanges = as(rowranges, "GRanges"),
    colData = coldata
)

#' 
## ---- purl=TRUE---------------------------------------------------------------
# wrong number of samples:

bad1 <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    rowRanges = as(rowranges, "GRanges"),
    colData = coldata[1:3,]
)

#' 
## ---- purl=TRUE---------------------------------------------------------------
# same number of genes but in different order:

bad2 <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts)),
  rowRanges = as(rowranges[c(2:nrow(rowranges), 1),], "GRanges"),
  colData = coldata
)


#' 
## ---- purl=TRUE---------------------------------------------------------------
# Access the counts
head(assay(se))
dim(assay(se))

# The above works now because we only have one assay, "counts"
# But if there were more than one assay, we would have to specify
# which one like so:

head(assay(se, "counts"))

# Access the sample annotations
colData(se)
dim(colData(se))

# Access the gene annotations
head(rowData(se))
dim(rowData(se))

# Make better sample IDs that show sex, time and mouse ID:

se$Label <- paste(se$sex, se$time, se$mouse, sep = "_")
se$Label
colnames(se) <- se$Label

# Our samples are not in order based on sex and time
se$Group <- paste(se$sex, se$time, sep = "_")
se$Group

# change this to factor data with the levels in order 
# that we want, then rearrange the se object:

se$Group <- factor(se$Group, levels = c("Female_Day0","Male_Day0", 
                                        "Female_Day4","Male_Day4",
                                        "Female_Day8","Male_Day8"))
se <- se[, order(se$Group)]
colData(se)

# Finally, also factor the Label column to keep in order in plots:

se$Label <- factor(se$Label, levels = se$Label)



#' 
## ---- purl=TRUE---------------------------------------------------------------
saveRDS(se, "data/GSE96870_se.rds")
rm(se) # remove the object!
se <- readRDS("data/GSE96870_se.rds")

#' 
#' ## Challenge: How to subset to mRNA genes
#' 
#' Before, we conceptually discussed subsetting to only the mRNA genes. Now that we have our `SummarizedExperiment` object, it becomes much easier to write the codes to subset `se` to a new object called `se_mRNA` that contains only the genes/rows where the `rowData(se)$gbkey` is equal to mRNA. Write the codes and then check you correctly got the 16,859 mRNA genes:
#' 
## ---- purl=TRUE---------------------------------------------------------------
mapIds(org.Mm.eg.db, keys = "497097", column = "SYMBOL", keytype = "ENTREZID")

#' 
## ---- purl=TRUE---------------------------------------------------------------
keys <- head(keys(hgu95av2.db, "ENTREZID"))
last <- function(x){x[[length(x)]]}

mapIds(hgu95av2.db, keys = keys, column = "ALIAS", keytype = "ENTREZID")

# When there is 1:many mapping, the default behavior was 
# to output the first match. This can be changed to a function,
# which we defined above to give us the last match:

mapIds(hgu95av2.db, keys = keys, column = "ALIAS", keytype = "ENTREZID", multiVals = last)

# Or we can get back all the many mappings:

mapIds(hgu95av2.db, keys = keys, column = "ALIAS", keytype = "ENTREZID", multiVals = "list")

#' 
## ---- purl=TRUE---------------------------------------------------------------
sessionInfo()

#' 
#' ## Exploratory analysis and quality control
#' 
#' ## Challenge: What kind of genes survived this filtering?
#' 
#' Last episode we discussed subsetting down to only mRNA genes. Here we subsetted based on a minimal expression level. How many of each type of gene survived the filtering?
#' 
#' ## Differential expression analysis
#' 
#' What will be the default **contrast**, **reference level** and **"last level"** for comparisons when running `results(dds)` for the example used in this lesson?
#' 
#' Explore the DE genes between males and females independent of time.
#' 
#' ::::::::::::::::::::::::::::::::::::: challenge
#' 
#' By default `independentFiltering` is set to `TRUE`. What happens without filtering lowly expressed genes? Use the `summary()` function to compare the results. Most of the lowly expressed genes are not significantly differential expressed (blue in the above MA plots). What could cause the difference in the results then?
#' 
#' Check the heatmap and top DE genes. Do you find something expected/unexpected?
#' 
#' ## Gene set enrichment analysis
#' 
#' Can you convert between different gene set representations? E.g. convert a
#' list to a two-column data frame?
#' 
#' ## Next steps
#' 
#' 
#' ## Extra exploration of design matrices
#' 
#' ### Challenge
#' 
#' Based on this visualization, would you say that the data set is balanced, or are there combinations of predictor variables that are severely over- or underrepresented?
#' 
#' ### Challenge
#' 
#' With this design, what is the interpretation of the `sexMale` coefficient?
#' 
#' ### Challenge
#' 
#' Set up the design formula to compare male and female spinal cord samples from Day0 as above, but instruct R to not include an intercept in the model. How does this change the interpretation of the coefficients? What contrast would have to be specified to compare the mean expression of a gene between male and female mice?
#' 
#' ### Challenge
#' 
#' Set up the design formula to compare the three time points (Day0, Day4, Day8) in the male spinal cord samples, and visualize it using `ExploreModelMatrix`.
#' 
