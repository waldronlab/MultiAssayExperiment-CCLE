drugDataFrame <- function(cellSplits, doseNames) {
    compressedListCandidates <- lapply(cellSplits, function(CL) {
        apply(CL, 2, function(x) all(grepl(",", x)))
    })
    drugData <- Map(function(x, y) {
        compListFrame <- x[y]
        compListFrame[] <- apply(compListFrame, 2L, function(col) {
            if (is.integer(col))
                IRanges::IntegerList(strsplit(col, ","))
            else if (is.numeric(col))
                IRanges::NumericList(strsplit(col, ","))
            else (is.character(col))
                IRanges::CharacterList(strsplit(col, ","))
        })
        compListFrame
    }, x = cellSplits, y = compressedListCandidates)
    allNames <- Reduce(intersect, lapply(cellSplits, names))
    otherVars <- allNames[!(allNames %in% doseNames)]
    ## FIX HERE
    otherCols <- lapply(cellSplits, function(cellLine) {
        otherSub <- cellLine[, otherVars]
        otherList <- lapply(otherSub, function(subCol) {
            if (length(unique(subCol)) == 1L) {
                return(unique(subCol))
            } else {
                names(subCol) <- cellLine[["Compound"]]
                if (is.integer(subCol))
                    subCol <- IRanges::IntegerList(subCol)
                else if (is.character(subCol))
                    subCol <- IRanges::CharacterList(subCol)
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

## FIX ME HERE
    doseVars <- lapply(newList, function(x) { DataFrame(do.call(cbind, x)) })
    doseVars <- do.call(rbind, doseVars)
    drugDat <- cbind(newDF, doseVars)
    rownames(drugDat) <- drugDat[, "CCLE.Cell.Line.Name"]
    drugDat <- drugDat[, -which(names(drugDat) == "CCLE.Cell.Line.Name")]
    return(drugDat)
}

