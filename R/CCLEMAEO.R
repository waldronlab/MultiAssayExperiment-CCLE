library(readr)
library(GSRI) # for readGct function
pdata <- read_csv("rawdata/CCLE_NP24.2009_Drug_data_2012.02.20.csv")
myMAF <- read_delim("rawdata/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf", delim = "\t", na = "<NA>")
newEx <- read_delim("rawdata/CCLE_Expression_Entrez_2012-09-29.gct", delim = "\t")
# NewMatrix <- readGct("rawdata/CCLE_Expression_Entrez_2012-09-29.gct")
fastMat <- read_delim("rawdata/CCLE_Expression_Entrez_2012-09-29.gct", delim = "\t", skip = 2)
rownames(fastMat) <- fastMat$Name
fastMat <- fastMat[, -which(names(fastMat) == "Name")]

copyNum <- read_delim("rawdata/CCLE_copynumber_byGene_2012-09-29.txt", delim = "\t")
