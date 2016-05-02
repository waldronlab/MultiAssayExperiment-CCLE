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
library(MultiAssayExperiment)


# orininal working directory
working0 <- getwd()

# data directories
data <- "~/Source/MultiAssayExperiment-CCLE/data"
#data <- "rawdata"

# set working directory as data directory
setwd(data)

# read in data sets
DNAcopyNumber <- read_delim("CCLE_copynumber_byGene_2012-09-29.txt", delim = "\t")
mRNAexpression <-read_delim("CCLE_Expression_Entrez_2012-09-29.gct", delim = "\t", skip = 2)
mutations <- read_delim("CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf", delim = "\t", na = "<NA>")
pData <- read_csv("CCLE_NP24.2009_Drug_data_2012.02.20.csv")

# add rownames to mRNAexpression
rownames(mRNAexpression) <- mRNAexpression$Name
mRNAexpression <- mRNAexpression[, -which(names(mRNAexpression) == "Name")]

# save objects as rsd files
saveRDS(DNAcopyNumber, file = paste0(data, "/DNAcopyNumber.rds"))
saveRDS(mRNAexpression, file = paste0(data, "/mRNAexpression.rds"))
saveRDS(mutations, file = paste0(data, "/mutations.rds"))
saveRDS(pData, file = paste0(data, "/pData.rds"))

# restore orininal working directory
setwd(working0)

# load objects from rsd files
readRDS(paste0(data, "/DNAcopyNumber.rds"))
readRDS(paste0(data, "/mRNAexpression.rds"))
readRDS(paste0(data, "/mutations.rds"))
readRDS(paste0(data, "/pData.rds"))
