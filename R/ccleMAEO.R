###################### ccleMAEO.R ######################
#
# A script for creating a MultiAssayExperiment Object
# (MAEO) from Cancer Cell Line Encyclopedia (CCLE) data
#
# Authors: Marcel Ramos, Lucas Schiffer
#
########################################################

# Load libraries ----------------------------------------------------------
library(readr)
library(S4Vectors)
library(TCGAutils) # https://github.com/waldronlab/TCGAutils
library(Biobase)
library(SummarizedExperiment)
library(GenomicRanges)
library(MultiAssayExperiment) # https://github.com/vjcitn/MultiAssayExperiment
library(RaggedExperiment)
library(downloader)

# Check files and directories ---------------------------------------------
dataURL <- "https://data.broadinstitute.org/ccle_legacy_data"
folders <- c("dna_copy_number", "hybrid_capture_sequencing",
    "mRNA_expression", "pharmacological_profiling")

filesOfInterest <- c("CCLE_copynumber_byGene_2013-12-03.txt",
    "CCLE_Expression_Entrez_2012-09-29.gct",
    "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf",
    "CCLE_NP24.2009_Drug_data_2015.02.24.csv")

if (!dir.exists("rawdata")) { dir.create("rawdata") }

for (i in seq_along(filesOfInterest)) {
    if (!file.exists(file.path("rawdata", filesOfInterest[i]))) {
        message(filesOfInterest[i], " not found in the 'rawdata' folder")
        download(file.path(dataURL, folders[i], filesOfInterest[i]),
                 destfile = file.path("rawdata", filesOfInterest[i]))
    }
}

# DNA Copy Number ---------------------------------------------------------
DNAcopyNumber <- read_delim("rawdata/CCLE_copynumber_byGene_2012-09-29.txt", delim = "\t")

# create a RangeSummarizedExperiment from DNAcopyNumber
DNAcopyNumber <- DataFrame(DNAcopyNumber)
rowRanges <- unlist(makeGRangesListFromDataFrame(DNAcopyNumber,
    names.field = "SYMBOL", seqnames.field = "CHR", start.field = "CHRLOC",
    end.field = "CHRLOCEND"), use.names = FALSE)
EGID <- DNAcopyNumber[, "EGID"]
# Extract SITE info from column names
# lapply(strsplit(colnames(DNAcopyNumber), "_"), `[`, 2L)
newRSE <- makeSummarizedExperimentFromDataFrame(DNAcopyNumber[,
    -which(names(DNAcopyNumber) %in% c("EGID", "SYMBOL"))],
    seqnames.field = "CHR", start.field = "CHRLOC", end.field = "CHRLOCEND")
rowRanges(newRSE) <- rowRanges
rowData(newRSE) <- EGID

# mRNA Expression Entrez --------------------------------------------------
mRNAexpression <- read_delim("rawdata/CCLE_Expression_Entrez_2012-09-29.gct", delim = "\t", skip = 2)

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


# Mutations ---------------------------------------------------------------
types <- c("c", "i", "c", "i", "c", "i", "i", "c", "c", "c", "c", "c",
  "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c",
  "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c",
  "c", "c", "c", "c", "c", "c", "c", "c", "c", "i", "i", "c", "c" )
types2 <- paste(types, collapse = "")
mutations <- read_delim("rawdata/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf",
                        delim = "\t", na = "<NA>", col_types = types2)
# create a GRangesList from mutations
newMut <- makeGRangesListFromDataFrame(as.data.frame(mutations, stringsAsFactors = FALSE),
                                      names.field = "Hugo_Symbol",
                                      split.field = "Tumor_Sample_Barcode",
                                      start.field = "Start_position",
                                      end.field = "End_position",
                                      keep.extra.columns = TRUE)
newMut <- RaggedExperiment(newMut)
genome(newMut) <- TCGAutils:::.getHGBuild("37")


# primary DataFrame -------------------------------------------------------
pData <- read_csv("rawdata/CCLE_NP24.2009_Drug_data_2012.02.20.csv")
pData <- DataFrame(pData)
splitData <- S4Vectors::split(pData, pData$CCLE.Cell.Line.Name)
source("R/drugDataFrame.R")
pData <- drugDataFrame(splitData, c("Doses..uM.", "Activity.Data..median.", "Activity.SD"))
pData$TissueOrigin <- gsub("^[^_]+_", "", rownames(pData), perl = TRUE)

# use IDs as the rownames of the pData # see line 29
# no function needed to translate cellLine names
# rownames(pData)


# Experiment List ---------------------------------------------------------
dataList <- list(CNA = newRSE, Mutations = newMut, mRNA = mRNAEset)

## Convenient function for removing the first part of the CCLE ID:
getCellID <- function(fullCCLE) {
    stopifnot(is.character(fullCCLE))
    cellID <- gsub("(^[^_]+)_\\w+", "\\1", fullCCLE, perl = TRUE)
    return(cellID)
}

# Map Creation ------------------------------------------------------------
CCMap <- generateMap(dataList, pData)


# Prepare the MultiAssayExperiment ----------------------------------------
newList <- PrepMultiAssay(dataList, pData, CCMap)


# Create MultiAssayExperiment ---------------------------------------------
ccleMAEO <- MultiAssayExperiment(newList$Elist, newList$pData, newList$sampleMap)
