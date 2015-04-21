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
                     }
                     if (!missing(infoFile)) self$infoFile <- infoFile
                   },
                   setInfoFile = function(val) {
                     self$infoFile <- val
                   },
                   setDataFile = function(val) {
                     self$dataFile <- val
                     self$nrow = as.numeric(unlist(strsplit(unlist(h5ls(self$dataFile)[5])[4], " x ")))[1]
                     self$ncol = as.numeric(unlist(strsplit(unlist(h5ls(self$dataFile)[5])[4], " x ")))[2]
                   },
                   getProbes = function() {
                     unlist(lapply(h5read(self$dataFile, "/0/META/ROW"), function(x) { gsub(" ", "", x) } ) )
                   },
                   getData = function(rows, cols) {
                     first = TRUE;
                     count = 1;
                     for(c in cols) {  
                       cat(paste("Processing column ", count, " of ", length(cols), "\n", sep=""))
                       count <- count+1
                       column <- h5read(self$dataFile, "/0/DATA/0", start=c(1,c), count=c(self$nrow,1))
                       if(first) { # there is a better way to do this I just don't know what it is
                         result <- as.data.frame(column[[1]][rows]) 
                         first = FALSE;
                       } else {  
                         result <- cbind(result, as.data.frame(column[[1]][rows]))                          
                       }
                     }
                     colnames(result) <- self$colsByIndex(cols)
                     rownames(result) <- self$getProbes()[rows]
                     result
                   },
                   rowsByProbe = function(probes) {                      
                     ix <- which(self$getProbes() %in% probes)
                     ix
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


go <- function() {
  t1 <- LINCS$new(dataFile = "/mnt/jovinge.lab/zspc_n1328098x22268.gctx");
  #  t1$colsById("ASG001_MCF7_24H_X1_B7_DUO52HI53LO:P12")
  t1$l1000ids()
}







