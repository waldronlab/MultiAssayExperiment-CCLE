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

# Create a list of concatenated compound and target name combinations
concatCol <- lapply(pData, function(CellLine) {
    apply(CellLine, 1, function(dataset) {
        paste(dataset["Compound"], dataset["Target"], sep = ":")
    })
})

numIndex <- seq_along(pData)
names(numIndex) <- names(pData)

# Add unique compound:target vector as pData column
pData2 <- lapply(numIndex, function(i, CellLine, comboTarget) {
    cbind(CellLine[[i]], DataFrame(CompTarget = comboTarget[[i]]))
}, CellLine = pData, comboTarget = concatCol)

# Get all possible drug-target combinations
combos <- Reduce(union, concatCol)

## Test that all are unique regardless of case
length(unique(toupper(combos))) == length(combos)

## Create a matrix by extracting one column from pData for all CellLines and
## drug-target combinations
getDrugMatrix <- function(dataList, targetColumn, features) {
FullList <- lapply(dataList, function(DF) {
    DF[, targetColumn][match(features, DF[, "CompTarget"])]
})
drugMatrix <- do.call(cbind, FullList)
rownames(drugMatrix) <- features
colnames(drugMatrix) <- names(dataList)
return(drugMatrix)
}

names(pData2[[1]])
getDrugMatrix(pData2, "EC50.uM.", combos)
