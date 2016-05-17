doseVariables <- function(cellSplits, doseNames) {
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
return(doseVars)
}
