makeRSE <- function(dataset) {
    grangeNames <- grep("name|chr|start|end$", names(dataset),
                        value = TRUE, ignore.case = TRUE)
    if (length(grangeNames) != 4L) {warning("not all rowRange columns found")}
    DFranges <- DataFrame(dataset[, grangeNames])
    chrname <- grep("chr", names(DFranges), value = TRUE, ignore.case = TRUE)
    genName <- grep("name", names(DFranges), value = TRUE, ignore.case = TRUE)
    if (all(grepl("^[0-9]+$", sample(DFranges[, chrname], size = 5)))) {
        DFranges[, chrname] <- paste0("chr", DFranges[, chrname])
    }
    if (!chrname %in% c("seqnames", "chr", "chrom")) {
        names(DFranges)[which(names(DFranges) == chrname)] <- "chrom"
    }
    RowRanges <- as(DFranges[,!grepl("name", grangeNames, ignore.case = TRUE)],
                    "GRanges")
    if(length(unique(DFranges[, genName])) == dim(DFranges)[1]) {
        names(RowRanges) <- DFranges[, genName]
        rNames <- DFranges[, genName]
    }
    dm <- dataset[, !(names(dataset) %in% grangeNames)]
    if(exists("rNames")) {
        rownames(dm) <- rNames
    }
    newSE <- SummarizedExperiment::SummarizedExperiment(
        assays = SimpleList(counts = dm), rowRanges = RowRanges)
    return(newSE)
}
