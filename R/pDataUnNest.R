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
names(pData[[1]])

combos <- Reduce(union, lapply(pData, function(CellLine) {
    apply(CellLine, 1, function(dataset) {
        paste(dataset["Compound"], dataset["Target"], sep = ":")
    })
}))

## Test all unique despite case
length(unique(toupper(combos))) == length(combos)

## Not working
getDrugMatrix <- function(dataList, drugCol) {
extractedVectors <- lapply(dataList,
                           function(DF, columnRequest) {
                               DF[, columnRequest]
                               }, columnRequest = drugCol)
fullVector <- Reduce(union, extractedVectors)
FullList <- lapply(extractedVectors,
                   function(CellLine, AllIn)
                       { AllIn[match(CellLine, AllIn)]}, AllIn = fullVector)
drugList <- do.call(cbind, FullList)
drugMatrix <- as.matrix(drugList)
return(drugMatrix)
}
