# load readr
library(readr)
#library(CePa)
#library(GSRI) # for readGct function

# orininal working directory
working0 <- getwd()

# data directory
data <- "~/Source/MultiAssayExperiment-CCLE/data"

# set working directory as data directory
setwd(data)

# read in data sets
pdata <- read_csv("CCLE_NP24.2009_Drug_data_2012.02.20.csv")
myMAF <- read_delim("CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf", delim = "\t", na = "<NA>")
newEx <- read_delim("CCLE_Expression_Entrez_2012-09-29.gct", delim = "\t")
#NewMatrix <- readGct("CCLE_Expression_Entrez_2012-09-29.gct")
#testM <- read.gct("CCLE_Expression_Entrez_2012-09-29.gct")
fastMat <- read_delim("CCLE_Expression_Entrez_2012-09-29.gct", delim = "\t", skip = 2)
copyNum <- read_delim("CCLE_copynumber_byGene_2012-09-29.txt", delim = "\t")

# add rownames to fastMat
rownames(fastMat) <- fastMat$Name
fastMat <- fastMat[, -which(names(fastMat) == "Name")]

# restore orininal working directory
setwd(working0)