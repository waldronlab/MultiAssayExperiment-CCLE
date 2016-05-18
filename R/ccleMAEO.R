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
library(TCGAmisc)
library(Biobase)
library(SummarizedExperiment)
library(GenomicRanges)
library(MultiAssayExperiment)

# read in data sets
if (!dir.exists("rawdata")) {dir.create("rawdata")}
DNAcopyNumber <- read_delim("rawdata/CCLE_copynumber_byGene_2012-09-29.txt", delim = "\t")
mRNAexpression <-read_delim("rawdata/CCLE_Expression_Entrez_2012-09-29.gct", delim = "\t", skip = 2)
mutations <- read_delim("rawdata/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf", delim = "\t", na = "<NA>")
pData <- read_csv("rawdata/CCLE_NP24.2009_Drug_data_2012.02.20.csv")
pData <- DataFrame(pData)
splitData <- S4Vectors::split(pData, pData$CCLE.Cell.Line.Name)
source("R/drugDataFrame.R")
pData <- drugDataFrame(splitData, c("Doses..uM.", "Activity.Data..median.", "Activity.SD"))
pData$TissueOrigin <- gsub("^[^_]+_", "", rownames(pData), perl = TRUE)

# add rownames to mRNAexpression
rownames(mRNAexpression) <- mRNAexpression$Name
mRNAexpression <- mRNAexpression[, -which(names(mRNAexpression) == "Name")]
annoteFeatures <- mRNAexpression[, "Description"]
annoteFeatures <- data.frame(annoteFeatures)

mRNAexpression <- mRNAexpression[, -which(names(mRNAexpression) == "Description")]
mRNAexpression <- as.matrix(mRNAexpression)

mRNAEx <- mRNAexpression[, colnames(mRNAexpression) %in% rownames(pData)]
rownames(annoteFeatures) <- rownames(mRNAEx)

mRNAEset <- ExpressionSet(assayData = mRNAEx, featureData = AnnotatedDataFrame(annoteFeatures))

# save objects as rsd files
if (!dir.exists("rdsdata")) {dir.create("rdsdata")}
saveRDS(DNAcopyNumber, file = "rdsdata/DNAcopyNumber.rds")
saveRDS(mRNAexpression, file = "rdsdata/mRNAexpression.rds")
saveRDS(mutations, file = "rdsdata/mutations.rds")
saveRDS(pData, file = "rdsdata/pData.rds")

# load objects from rsd files
DNAcopyNumber <- readRDS("rdsdata/DNAcopyNumber.rds")
mRNAexpression <- readRDS("rdsdata/mRNAexpression.rds")
mutations <- readRDS("rdsdata/mutations.rds")
pData <- readRDS("rdsdata/pData.rds")

# use IDs as the rownames of the pData # see line 29
# no function needed to translate cellLine names
rownames(pData)


# create a ExpressionSet or RangeSummarizedExperiment from
# DNAcopyNumber, available from SummarizedExperiment package
source("R/makeRSE.R")
newRSE <- makeRSE(DNAcopyNumber)

# create a GRangesList from mutations
?GRangesList()

# split the tummor sample barcode and make list thereof
# duplicate makeGRangesList method from TCGAmics package
?makeGRangesList()

# generate CCLE MAEO by passing these arguments to the MultiAssayExperiment method
# specimine ID (short spec), assay (long spec), assayname (e.g. DNAcopyNumber)
?MultiAssayExperiment()