###################### ccleMAEO.R ######################
#
# A script for creating a MultiAssayExperiment Object
# (MAEO) from Cancer Cell Line Encyclopedia (CCLE) data
#
# Authors: Marcel Ramos
# Contributor: Lucas Schiffer
#
########################################################

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages({
    library(readr)
    library(S4Vectors)
    library(TCGAutils) # https://github.com/waldronlab/TCGAutils
    library(Biobase)
    library(SummarizedExperiment)
    library(GenomicRanges)
    library(MultiAssayExperiment)
    library(RaggedExperiment)
    library(downloader)
})

# Check files and directories ---------------------------------------------
dataURL <- "https://data.broadinstitute.org/ccle_legacy_data"
folders <- c("dna_copy_number", "mRNA_expression",
    "hybrid_capture_sequencing", "pharmacological_profiling")

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
DNAcopyNumber <- read_delim("rawdata/CCLE_copynumber_byGene_2013-12-03.txt", delim = "\t")

# create a RangeSummarizedExperiment from DNAcopyNumber
DNAcopyNumber <- DataFrame(DNAcopyNumber)
rowRanges <- unlist(makeGRangesListFromDataFrame(DNAcopyNumber,
    names.field = "SYMBOL", seqnames.field = "CHR", start.field = "CHRLOC",
    end.field = "CHRLOCEND"), use.names = FALSE)

# Save EntrezGene IDs for rowData
EGID <- DNAcopyNumber[, "EGID"]

otherNames <- c("EGID", "SYMBOL")
# Take only relevant names
DNAcopyNumber <- DNAcopyNumber[, !colnames(DNAcopyNumber) %in% otherNames]

# Find necessary range info variables
rangeVars <- -which(names(DNAcopyNumber) %in% c("CHR", "CHRLOC", "CHRLOCEND"))

# Extract SITE info from column names
SITE <- gsub("(^[^_]+)_(\\w+$)", "\\2",
    colnames(DNAcopyNumber[, rangeVars]), perl = TRUE)

# Create SummarizedExperiment
newRSE <- makeSummarizedExperimentFromDataFrame(DNAcopyNumber,
    seqnames.field = "CHR", start.field = "CHRLOC", end.field = "CHRLOCEND")

# Add rowRanges > RangedSummarizedExperiment
rowRanges(newRSE) <- rowRanges

# Annotate
rowData(newRSE) <- DataFrame(EGID)
colData(newRSE) <- DataFrame(SITE, row.names = colnames(DNAcopyNumber[, rangeVars]))
assayNames(newRSE) <- "copyNumber"

# mRNA Expression Entrez --------------------------------------------------
mRNAexpression <- read_delim("rawdata/CCLE_Expression_Entrez_2012-09-29.gct",
    delim = "\t", skip = 2)

# add rownames to mRNAexpression
mRNAexpression <- DataFrame(mRNAexpression)
rownames(mRNAexpression) <- mRNAexpression[["Name"]]
mRNAexpression <- mRNAexpression[, -which(names(mRNAexpression) == "Name")]

rowDat <- DataFrame(description = mRNAexpression[, "Description"])
rownames(rowDat) <- rownames(mRNAexpression)

mRNAexpression <- mRNAexpression[, -which(names(mRNAexpression) == "Description")]
mRNAexpression <- as.matrix(mRNAexpression)

mRNASE <- SummarizedExperiment(mRNAexpression, rowData = rowDat)

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
genome(newMut) <- TCGAutils:::.getHGBuild("37")
newMut <- RaggedExperiment(newMut)

# primary DataFrame -------------------------------------------------------
colData <- read_csv("rawdata/CCLE_NP24.2009_Drug_data_2015.02.24.csv")
colData <- DataFrame(colData)
splitData <- S4Vectors::split(colData, colData[["CCLE.Cell.Line.Name"]])
drugVars <- c("Compound", "EC50..uM.", "IC50..uM.", "Amax", "ActArea")
doseVars <- c("Compound", "Doses..uM.", "Activity.Data..median.", "Activity.SD")


getDrugData <- function(dataList, variables) {
    allCompounds <- sort(Reduce(union, lapply(dataList, function(x)
        x[["Compound"]])
    ))
    allCompounds <- DataFrame(Compound = allCompounds)
    splitData <- lapply(dataList, function(x) {
        drugDF <- merge(allCompounds, x[, variables],
                        by = "Compound", all = TRUE)
        rownames(drugDF) <- drugDF[["Compound"]]
        drugDF[, -which(colnames(drugDF) == "Compound")]
    })

    varVect <- Reduce(unique, lapply(splitData, seq_along))
    names(varVect) <- Reduce(intersect, lapply(splitData, names))

    lapply(varVect, function(y) {
        do.call(cbind, lapply(seq_along(splitData), function(i, x, y) {
            col <- x[[i]][, y, drop = FALSE]
            names(col) <- names(splitData[i])
            col
        }, x = splitData, y = y))
    })
}

## Produce drug array
DrugList <- getDrugData(splitData, drugVars)

nrows <- unique(vapply(DrugList, nrow, numeric(1L)))
ncols <- unique(vapply(DrugList, ncol, numeric(1L)))
llen <- length(DrugList)

allrnames <- vapply(DrugList, rownames, character(nrows))
dfrnames <- as.data.frame(allrnames)
if (all(vapply(dfrnames[, seq_along(dfrnames)[-1]],
               function(x) identical(dfrnames[[1]], x), logical(1L))))
    rnames <- allrnames[, 1L, drop = TRUE]

allcnames <- vapply(DrugList, colnames, character(ncols))
dfcnames <- as.data.frame(allcnames)
if (all(vapply(dfcnames[, seq_along(dfcnames)[-1]],
               function(x) identical(dfcnames[[1]], x), logical(1L))))
    cnames <- allcnames[, 1L, drop = TRUE]

dnames <- names(DrugList)

drugArray <- array(unlist(lapply(DrugList, data.matrix)),
                   dim = c(nrows, ncols, llen),
                   dimnames = list(rnames, cnames, dnames))

DoseList <- getDrugData(splitData, doseVars)

.replaceShorts <- function(splitList) {
    shorts <- lengths(splitList) < 8L &
        lengths(splitList) > 0L &
        !is.na(splitList)
    nas <- is.na(splitList)
    if (any(shorts)) {
    allVals <- Reduce(union, splitList[!shorts])
    elemIdx <- which(shorts)
    for (idx in elemIdx) {
        notInShorts <- !allVals %in% splitList[[idx]]
        valVec <- rep(NA, length(allVals))
        valVec[notInShorts] <- NA
        valVec[!notInShorts] <- splitList[[idx]]
        splitList[[idx]] <- valVec
        }
    }
    if (any(nas)) {
        elemIdx <- which(nas)
        lenSplit <- unique(lengths(splitList[!nas]))
        stopifnot(S4Vectors::isSingleInteger(lenSplit))
        valVect <- rep(NA, lenSplit)
        for (idx in elemIdx) {
            splitList[[idx]] <- valVect
        }
    }
    splitList
}

# TODO Fix BUG: Some dose concentrations have 9 values instead of 8
# for (i in seq_along(DoseList)) {
    i <- 1L
    lapply(DoseList[[i]], function(x) {
        splitList <- lapply(strsplit(x, ","), as.numeric)
        splitList <- .replaceShorts(splitList = splitList)
        midDF <- as.data.frame(splitList, col.names = seq_len(length(splitList)))
        plyr::aaply(t(midDF), 1L, as.numeric)
    })
# }


# source("R/drugDataFrame.R")
# colData <- drugDataFrame(splitData, doseVars)

## FIX ME see drugDataFrame()
colData[["TissueOrigin"]] <- gsub("^[^_]+_", "", rownames(colData), perl = TRUE)

# use IDs as the rownames of the colData # see line 29
# no function needed to translate cellLine names
# rownames(colData)


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
ccleMAEO <- do.call(MultiAssayExperiment, newList)
