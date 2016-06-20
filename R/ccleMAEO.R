###################### ccleMAEO.R ######################
#
# A script for creating a MultiAssayExperiment Object
# (MAEO) from Cancer Cell Line Encyclopedia (CCLE) data
#
# Authors: Marcel Ramos, Lucas Schiffer
#
########################################################

# load necessary packages
library(readr)
library(S4Vectors)
library(BiocInterfaces) # https://github.com/waldronlab/BiocInterfaces
library(Biobase)
library(SummarizedExperiment)
library(GenomicRanges)
library(MultiAssayExperiment) # https://github.com/vjcitn/MultiAssayExperiment

# read in data sets
if (!dir.exists("rawdata")) {dir.create("rawdata")}
if (!file.exists("rawdata/CCLE_copynumber_byGene_2012-09-29.txt")) {stop("The file CCLE_copynumber_byGene_2012-09-29.txt must exist in the rawdata directory. \n It can be downloaded from http://www.broadinstitute.org/ccle/data/browseData")}
if (!file.exists("rawdata/CCLE_Expression_Entrez_2012-09-29.gct")) {stop("The file CCLE_Expression_Entrez_2012-09-29.gct must exist in the rawdata directory. \n It can be downloaded from http://www.broadinstitute.org/ccle/data/browseData")}
if (!file.exists("rawdata/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf")) {stop("The file CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf must exist in the rawdata directory. \n It can be downloaded from http://www.broadinstitute.org/ccle/data/browseData")}
if (!file.exists("rawdata/CCLE_NP24.2009_Drug_data_2012.02.20.csv")) {stop("The file CCLE_NP24.2009_Drug_data_2012.02.20.csv must exist in the rawdata directory. \n It can be downloaded from http://www.broadinstitute.org/ccle/data/browseData")}

DNAcopyNumber <- read_delim("rawdata/CCLE_copynumber_byGene_2012-09-29.txt", delim = "\t")
mRNAexpression <- read_delim("rawdata/CCLE_Expression_Entrez_2012-09-29.gct", delim = "\t", skip = 2)
mutations <- read_delim("rawdata/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf", delim = "\t", na = "<NA>")
pData <- read_csv("rawdata/CCLE_NP24.2009_Drug_data_2012.02.20.csv")
pData <- DataFrame(pData)
splitData <- S4Vectors::split(pData, pData$CCLE.Cell.Line.Name)
source("R/drugDataFrame.R")
pData <- drugDataFrame(splitData, c("Doses..uM.", "Activity.Data..median.", "Activity.SD"))
pData$TissueOrigin <- gsub("^[^_]+_", "", rownames(pData), perl = TRUE)

# add rownames to mRNAexpression
mRNAexpression <- DataFrame(mRNAexpression)
rownames(mRNAexpression) <- mRNAexpression$Name
mRNAexpression <- mRNAexpression[, -which(names(mRNAexpression) == "Name")]
annoteFeatures <- mRNAexpression[, "Description"]
annoteFeatures <- data.frame(annoteFeatures)

mRNAexpression <- mRNAexpression[, -which(names(mRNAexpression) == "Description")]
mRNAexpression <- as.matrix(mRNAexpression)

mRNAEx <- mRNAexpression[, colnames(mRNAexpression) %in% rownames(pData)]
rownames(annoteFeatures) <- rownames(mRNAEx)

mRNAEset <- ExpressionSet(assayData = mRNAEx, featureData = AnnotatedDataFrame(annoteFeatures))

# use IDs as the rownames of the pData # see line 29
# no function needed to translate cellLine names
# rownames(pData)

# create a RangeSummarizedExperiment from DNAcopyNumber
newRSE <- makeRangedSummarizedExperimentFromDataFrame(DataFrame(DNAcopyNumber), seqnames.field = "NumChr",
                                                      start.field = "txStart", end.field = "txEnd")

# create a GRangesList from mutations
newMut <- makeGRangesListFromDataFrame(as.data.frame(mutations, stringsAsFactors = FALSE),
                                        partitioning.field = "Tumor_Sample_Barcode")
newMut <- RangedRaggedAssay(newMut)

dataList <- list(CNA = newRSE, Mutations = newMut, mRNA = mRNAEset)

## Convenient function for removing the first part of the CCLE ID:
getCellID <- function(fullCCLE) {
    stopifnot(is.character(fullCCLE))
    cellID <- gsub("(^[^_]+)_\\w+", "\\1", fullCCLE, perl = TRUE)
    return(cellID)
}

CCMap <- generateMap(dataList, pData)
newList <- PrepMultiAssay(dataList, pData, CCMap)

# call constructor, passing Elist, pData, and sampleMap arguments
ccleMAEO <- MultiAssayExperiment(newList$Elist, newList$pData, newList$sampleMap)
