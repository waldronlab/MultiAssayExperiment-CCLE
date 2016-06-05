# save objects as rsd files
if (!dir.exists("rdsdata")) {dir.create("rdsdata")}
saveRDS(DNAcopyNumber, file = "rdsdata/DNAcopyNumber.rds")
saveRDS(mRNAexpression, file = "rdsdata/mRNAexpression.rds")
saveRDS(mutations, file = "rdsdata/mutations.rds")
saveRDS(pData, file = "rdsdata/pData.rds")
