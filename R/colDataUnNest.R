## Example for nesting colData
library(readr)
library(S4Vectors)

colData <- read_csv("rawdata/CCLE_NP24.2009_Drug_data_2015.02.24.csv")
names(colData) <- gsub(" ", "", names(colData))

## Un-Nest approach

colData <- DataFrame(colData)
CellLineName <- colData$CCLECellLineName
colData <- colData[, -which(c("CCLECellLineName", "PrimaryCellLineName") %in% names(colData))]

colData <- S4Vectors::split(colData, CellLineName)
class(colData)

colData
names(colData[[1]])

# Create a list of concatenated compound and target name combinations
concatCol <- lapply(colData, function(CellLine) {
    apply(CellLine, 1, function(dataset) {
        paste(dataset["Compound"], dataset["Target"], sep = ":")
    })
})

numIndex <- seq_along(colData)
names(numIndex) <- names(colData)

# Add unique compound:target vector as colData column
colData2 <- lapply(numIndex, function(i, CellLine, comboTarget) {
    cbind(CellLine[[i]], DataFrame(CompTarget = comboTarget[[i]]))
}, CellLine = colData, comboTarget = concatCol)

# Get all possible drug-target combinations
combos <- Reduce(union, concatCol)

## Test that all are unique regardless of case
length(unique(toupper(combos))) == length(combos)

## Create a matrix by extracting one column from colData for all CellLines and
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

names(colData2[[1]])
getDrugMatrix(colData2, "EC50.uM.", combos)
combos
