drugDataFrame <- function(cellSplits, doseNames) {
    otherVars <- !(names(cellSplits[[1]]) %in% doseNames)
    otherCols <- lapply(cellSplits, function(cellLine) {
        otherSub <- cellLine[otherVars]
        otherList <- lapply(otherSub, function(subCol) {
            if (length(unique(subCol)) == 1L) {
                return(unique(subCol))
            } else {
                names(subCol) <- cellLine$Compound
                subCol <- S4Vectors::SimpleList(subCol)
                return(subCol)
            }
        })
        return(otherList)
    })

    newDF <- do.call(rbind, lapply(otherCols, function(x) {
        S4Vectors::DataFrame(x)
        })
    )

    newList <- lapply(splitData, function(cellLine, doseColumns) {
        vars <- vector("list", length(doseColumns))
        names(vars) <- doseColumns
        for (var in doseColumns) {
            groups <- rep(seq_len(nrow(cellLine)), cellLine$Num.Data)
            doseVector <- rapply(strsplit(cellLine[, var], ","), as.numeric)
            compList <- split(doseVector, groups)
            names(compList) <- cellLine$Compound
            vars[[var]] <- S4Vectors::SimpleList(compList)
        }
        return(vars)
    }, doseColumns = doseNames)

    doseVars <- lapply(newList, function(x) { DataFrame(do.call(cbind, x)) })
    doseVars <- do.call(rbind, doseVars)
    drugDat <- cbind(newDF, doseVars)
    rownames(drugDat) <- drugDat[, "CCLE.Cell.Line.Name"]
    drugDat <- drugDat[, -which(names(drugDat) == "CCLE.Cell.Line.Name")]
    return(drugDat)
}

