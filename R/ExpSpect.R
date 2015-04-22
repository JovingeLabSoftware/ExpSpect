library(R6)
library("rhdf5")
library("hgu133plus2.db")

LINCS <- R6Class("LINCS",
   public = list(
     dataFile = NA,
     infoFile = NA,
     nrow = NA,
     ncol = NA,
     
     initialize = function(dataFile, infoFile) {
       if (!missing(dataFile)) {
         self$setDataFile(dataFile)
       } else {
         self$setDataFile('/mnt/lincs/zspc_n1328098x22268.gctx')         
       }

       if (!missing(infoFile)) {
         self$infoFile <- infoFile         
       } else {
         self$infoFile <- '/mnt/lincs/inst_info.gctx'
       }
     },

     setInfoFile = function(val) {
       self$infoFile <- val
     },

     setDataFile = function(val) {
       self$dataFile <- val
       loc <- H5Fopen(self$dataFile)
       ds <- H5Dopen(loc, "0/DATA/0/matrix")
       dim <- H5Sget_simple_extent_dims(H5Dget_space(ds))$size
       self$nrow = dim[1]
       self$ncol = dim[2]
     }

   ),
   private = list(
     getMetaRowNames = function() {
       h5read(self$infoFile, "/0/META/ROW")[[1]]
     }
   )
)


LINCS$set("public", "colsByIndex", function(cols) {
  result = numeric()
  ncol = as.numeric(unlist(strsplit(unlist(h5ls(self$dataFile)[5])[4], " x ")))[2]
  for(i in 1:length(cols)) {
    name <- h5read(self$dataFile, "/0/META/COL", start=cols[i], count=1)
    result <- c(result, name)
  }  
  result
})

LINCS$set("public", "getProbes", function() {
  unlist(lapply(h5read(self$dataFile, "/0/META/ROW"), function(x) { gsub(" ", "", x) } ) )
})

LINCS$set("public", "getData", function(rows, cols) {
  first = TRUE
  count = 1
  for(c in cols) {  
    cat(paste("Processing column ", count, " of ", length(cols), "\n", sep=""))
    count <- count+1
    column <- h5read(self$dataFile, "/0/DATA/0", start=c(1,c), count=c(self$nrow,1))
    if(first) { # there is a better way to do this I just don't know what it is
      result <- as.data.frame(column[[1]][rows]) 
      first = FALSE
    } else {  
      result <- cbind(result, as.data.frame(column[[1]][rows]))                          
    }
  }
  colnames(result) <- self$colsByIndex(cols)
  rownames(result) <- self$getProbes()[rows]
  result
})


LINCS$set("public", "getMetaData", function(rows, cols) {
  first = TRUE;
  count = 1;
  nrow = length(private$getMetaRowNames())
  print(nrow)
  for(c in cols) {  
    cat(paste("Processing column ", count, " of ", length(cols), "\n", sep=""))
    count <- count+1
    column <- h5read(self$infoFile, "/0/DATA/0/matrix", start=c(1,c), count=c(nrow,1))
    if(first) { # there is a better way to do this I just don't know what it is
      result <- as.data.frame(column[[1]][rows]) 
      first = FALSE
    } else {  
      result <- cbind(result, as.data.frame(column[[1]][rows]))                          
    }
    H5close()
  }
  colnames(result) <- self$colsByIndex(cols)
  #rownames(result) <- private$getMetaRowNames()
  result
})


LINCS$set("public", "rowsByProbe", function(probes) {                      
  ix <- which(self$getProbes() %in% probes)
  ix
})


LINCS$set("public", "colsById", function(ids) {
  result = numeric()
  # let's take this in manageable chunks...
  cat("Looking up columns in > 1 million columns of data...this will take a few seconds.\n")
  for(i in 1:ceiling(self$ncol/100000)) {
    start <- ((i-1)*100000)+1
    count = min(100000, self$ncol-start)
    names <- unlist(lapply(h5read(self$dataFile, "/0/META/COL", start=start, count=count), function(x) { gsub(" ", "", x) } ) )
    result <- c(result, (start-1) + which(names %in% ids))
  }  
  H5close()
  result
})


LINCS$set("public", "l1000ids", function() {
  library(hgu133plus2.db)

  entrezids <- as.list(hgu133plus2ENTREZID)
  l1000_md <- read.delim("data/Landmark_Genes_n978.csv", sep=",", header=TRUE, as.is=TRUE)

  entrezids_l1000 <- entrezids[which(entrezids %in% l1000_md$Entrez.Gene.ID)]
  ix <- match(entrezids_l1000, l1000_md$Entrez.Gene.ID)
  probesets_l1000 <- cbind(as.data.frame(names(entrezids_l1000), stringsAsFactors = FALSE), as.data.frame(unlist(entrezids_l1000), stringsAsFactors = FALSE), l1000_md[ix,])[,-5]

  # add row numbers from data matrix in gctx file for convenience
  probesets_l1000 <- cbind(match(probesets_l1000[,1], self$getProbes()), probesets_l1000);
  ixna <- which(is.na(probesets_l1000[,1]))
  probesets_l1000 <- probesets_l1000[-ixna,]  
  probesetx_l1000 <- probesets_l1000[-grep("_x_", probesets_l1000$probeset),]
  rownames(probesets_l1000) <- 1:nrow(probesets_l1000)
  colnames(probesets_l1000) <- c("row", "probeset", "entrez_id", "L1000_id", "symbol", "title")
  probesetx_l1000 <- probesets_l1000[-grep("_x_", probesets_l1000$probeset),]
  rownames(probesets_l1000) <- 1:nrow(probesets_l1000)  
  probesets_l1000

})

### summarize
# ix <- match(rownames(data), l1000$probeset)
# summarized <- apply(data, 2, function(x) { tapply(x, l1000$entrez_id[ix], mean) } )



LINCSUTIL <- R6Class("LINCSUTIL",
   public = list(
     dataFile = NA,
     infoFile = NA,
     nrow = NA,
     ncol = NA,
     initialize = function() {
     }
   )
)

LINCSUTIL$set("public", "info2hdf5", function(infile="/mnt/lincs/inst.info", 
                                              outfile="/mnt/lincs/inst_info.gctx") {
  
  cat("Loading datafile (this will take a minute)\n")
  info <- read.delim(infile, sep="\t", header=TRUE, as.is=TRUE)
  
  h5createFile(outfile)
  h5createGroup(outfile, "0")
  h5createGroup(outfile, "0/DATA")
  h5createGroup(outfile, "0/DATA/0")
  h5createGroup(outfile, "0/META")
  h5createGroup(outfile, "0/META/COL")
  h5createGroup(outfile, "0/META/ROW")
  
  # we are transposing to match the structure of the gene expression matrix
  h5createDataset(outfile, "0/META/COL/id", 
                  c(nrow(info)), 
                  size=46, 
                  chunk=1000, 
                  level=7, storage.mode = "character")
  h5write(info[,1], outfile, "0/META/COL/id")
  
  h5createDataset(outfile, "0/META/ROW/id", 
                  c(ncol(info)), 
                  size=37, 
                  storage.mode = "character")
  h5write(colnames(info), outfile, "0/META/ROW/id")
  
  cat("Transposing and saving into hdf5 format (this will also take a minute)\n")
  h5createDataset(outfile, "0/DATA/0/matrix", 
                  dims = c(ncol(info), nrow(info)), 
                  size=50, 
                  chunk=c(1, 20000), 
                  level=7, storage.mode = "character")
  
  self$h5WriteMatrix(outfile, t(info), "0/DATA/0/matrix")
})

# writes in column slabs to circumvent segfault with huge data set (>20M datapoints) 
# dataset "name" must already have been created (i.e. by h5createDataset)
LINCSUTIL$set("public", "h5WriteMatrix", function(hdffile, data, name, colchunk=250000) {
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
})



