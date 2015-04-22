library(rhdf5)

info2hdf5 <- function(infile="/mnt/lincs/inst.info", outfile="/mnt/lincs/inst_info.gctx") {

  cat("Loading datafile\n")
  info <- read.delim(infile, sep="\t", header=TRUE, as.is=TRUE)

  unlink(outfile)
  h5createFile(outfile)
  h5createGroup(outfile, "0")
  h5createGroup(outfile, "0/DATA")
  h5createGroup(outfile, "0/DATA/0")
  h5createGroup(outfile, "0/META")
  h5createGroup(outfile, "0/META/COL")
  h5createGroup(outfile, "0/META/ROW")

  # we are transposing to match the structure of the gene expression matrix
  h5createDataset(outfile, "0/META/COL/id", c(nrow(info)), size=46, chunk=1000, level=7, storage.mode = "character")
  h5write(info[,1], outfile, "0/META/COL/id")
  h5createDataset(outfile, "0/META/ROW/id", c(ncol(info)), size=37, storage.mode = "character")
  h5write(colnames(info), outfile, "0/META/ROW/id")

  cat("Transposing and saving into hdf5 format\n")
  h5createDataset(outfile, "0/DATA/0/matrix", dims = c(ncol(info), nrow(info)), size=50, chunk=c(1, 20000), level=7, storage.mode = "character")
  h5WriteMatrix(outfile, info2, "0/DATA/0/matrix")
}

# writes in column slabs to circumvent segfault with huge data set (>20M datapoints) 
# dataset "name" must already have been created (i.e. by h5createDataset)
h5WriteMatrix <- function(hdffile, data, name, colchunk=250000) {
  if(colchunk > ncol(data)) {
    colchunk = ncol(data)
  }  
  for (i in 1:ceiling(ncol(data)/colchunk)) {
    start <- (i-1)*colchunk+1
    end <- min(ncol(data), i*colchunk)
    print(paste(i, start, end, sep=":"))
    h5write(data[,c(start:end)], hdffile, name, index=list(NULL,start:end))  
    H5close()
  }
}
