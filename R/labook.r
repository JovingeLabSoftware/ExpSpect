library(rhdf5)

info2hdf5 <- function(file="inst_info.gctx") {
  unlink(file)
  h5createFile(file)
  h5createGroup(file, "0")
  h5createGroup(file, "0/DATA")
  h5createGroup(file, "0/DATA/0")
  h5createGroup(file, "0/META")
  h5createGroup(file, "0/META/COL")
  h5createGroup(file, "0/META/ROW")
  # we are transposing to match the structure of the gene expression matrix
  h5createDataset(file, "0/META/COL/id", c(nrow(info)), size=46, chunk=1000, level=9, storage.mode = "character")
  h5write(info[,1], file, "0/META/COL/id")
  h5createDataset(file, "0/META/ROW/id", c(ncol(info)), size=37, storage.mode = "character")
  h5write(colnames(info), file, "0/META/ROW/id")
  h5createDataset(file, "0/DATA/0/matrix", c(ncol(info), nrow(info)), size=100, 
                  storage.mode = "character", chunk=c(1, 10000), level=9)

  # make sure you transpose info first
  h5write(info, file, "0/DATA/0/matrix")
}