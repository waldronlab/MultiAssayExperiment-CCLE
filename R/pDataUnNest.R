## Example for nesting pData
library(readr)
library(S4Vectors)

pData <- read_csv("rawdata/CCLE_NP24.2009_Drug_data_2012.02.20.csv")
names(pData) <- gsub(" ", "", names(pData))

## Un-Nest approach

pData <- DataFrame(pData)
CellLineName <- pData$CCLECellLineName
pData <- pData[, -which(c("CCLECellLineName", "PrimaryCellLineName") %in% names(pData))]

pData <- S4Vectors::split(pData, CellLineName)
class(pData)

pData
