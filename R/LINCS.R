#' LINCS class
#' 
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @description A class to facilitate working with the level 4 data from the LINCS project.  This 
#' class does not provide any of the data (which must be obtained under individual agreement with 
#' LINCS).  However, once obtained, working with their large binary data files will is facilitated 
#' by this class.
#'
#' @format \code{\link{LINCS}} class generator
#'
#' @usage \code{lincs = LINCS$new(dataFile, infoFile)}
#' \code{ids <- lincs$getProbes()[1:3]}
#'
#' @param dataFile path to the level for LINCS data file in hdf5 format (.gctx)
#' @param infoFile path to the instance info data file for the LINCS data file \emph{in hdf5 format (.gctx)}. \strong{Note:} 
#' LINCS provides this data in a tab delimited format (inst.info).  However that file is unwieldy to work with 
#' in that format as it takes a long time to load and takes a considerable chunk of memory.  It is much 
#' more efficient to convert it to an hdf5 file with a structure parallel to that of the datafile.  This class
#' provides a utility function \code{info2hdf5} that will accomplish this for you.  See methods section below.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{getProbes()}}{Returns a vector of 
#'       \href{http://www.bioconductor.org/packages/release/data/annotation/html/hgu133plus2.db.html}{HG
#'        U133 Plus 2.0} probesets for which data has been imputed from the 978
#'        landmark genes measured directly on the luminex bead arrays by the LINCS }
#' }
#'
#' @keywords data
#' 
LINCS <- R6Class("LINCS",
                 public = list(
                   dataFile = NA,
                   infoFile = NA,
                   nrow = NA,
                   ncol = NA,
                   metadataRows = NA,
                   
                   initialize = function(dataFile, infoFile) {
                     if (!missing(dataFile)) {
                       self$setDataFile(dataFile)
                     } else {
                       self$setDataFile('/mnt/lincs/zspc_n1328098x22268.gctx')         
                     }
                     
                     if (!missing(infoFile)) {
                       self$infoFile <- infoFile         
                     } else {
                       self$setInfoFile('/mnt/lincs/inst_info.gctx')         
                     }
                   },
                   
                   setInfoFile = function(val) {
                     self$infoFile <- val
                     loc <- H5Fopen(self$infoFile)
                     ds <- H5Dopen(loc, "0/DATA/0/matrix")
                     dim <- H5Sget_simple_extent_dims(H5Dget_space(ds))$size
                     self$metadataRows = dim[1]
                     H5close()
                   },
                   
                   setDataFile = function(val) {
                     self$dataFile <- val
                     loc <- H5Fopen(self$dataFile)
                     ds <- H5Dopen(loc, "0/DATA/0/matrix")
                     dim <- H5Sget_simple_extent_dims(H5Dget_space(ds))$size
                     self$nrow = dim[1]
                     self$ncol = dim[2]
                     private$colset = 1:self$ncol
                     private$rowset = 1:self$nrow
                     H5close()
                     private$.rownames = gsub(" ", "", h5read(self$dataFile, "/0/META/ROW/id"))[private$rowset]                     
                   }
                   
                   
                 ),
                 private = list(
                   getMetaRowNames = function() {
                     rn <- h5read(self$infoFile, "/0/META/ROW")[[1]]
                     H5close()
                     rn
                   },
                   rowset = NA,
                   colset = NA,
                   .data = NULL,
                   .metadata = NULL,    
                   .rownames = NULL,
                   geneids = 'probeset',
                   dataIsStale = TRUE,
                   metadataIsStale = TRUE
                 )
)

# getProbes
#
# returns the (imputed) affymetrix HG U133 Plus 2.0 probeset ids from the 
# gctx datafile.
#
LINCS$set("public", "getGeneIds", function() {
  return(private$.rownames)
  ps
})

LINCS$set("private", "getProbes", function() {
  probes <- gsub(" ", "", h5read(self$dataFile, "/0/META/ROW/id"))
  H5close()
  return(probes)
})

LINCS$set("public", "setId", function(id=c('probeset', 'entrez')) {
    private$geneids = id;
    private$.rownames = gsub(" ", "", h5read(self$dataFile, "/0/META/ROW/id"))[private$rowset]
    H5close()
    if(private$geneids == 'entrez') {
      entrezids <- as.list(hgu133plus2ENTREZID)
      private$.rownames = unique(sort(unlist(entrezids[private$.rownames])))
    } 
    private$dataIsStale = TRUE;    
})


LINCS$set("public", "setCols", function(cols=NA) {
  if(length(cols)==1 && is.na(cols)) {
    loc <- H5Fopen(self$dataFile)
    ds <- H5Dopen(loc, "0/DATA/0/matrix")
    dim <- H5Sget_simple_extent_dims(H5Dget_space(ds))$size
    cols = 1:dim[2]
    H5close()
  }
  private$colset = cols
  self$ncol = length(cols)
  private$dataIsStale = TRUE;
  private$metadataIsStale = TRUE;
})

LINCS$set("public", "setRows", function(rows=NA) {
  if(length(rows)==1 && is.na(rows)) {
    loc <- H5Fopen(self$dataFile)
    ds <- H5Dopen(loc, "0/DATA/0/matrix")
    dim <- H5Sget_simple_extent_dims(H5Dget_space(ds))$size
    rows = 1:dim[1]
    H5close()
  }
  private$rowset = rows
  self$nrow = length(rows)
  self$setId(private$geneids) # refresh rownames
  private$dataIsStale = TRUE;
  private$metadataIsStale = TRUE;
})

LINCS$set("public", "data", function(rows) {
  if(private$dataIsStale) {
    private$loadData()
    private$dataIsStale = FALSE;
  }
  if(private$metadataIsStale) {
    private$loadMetadata()
    private$metadataIsStale = FALSE;
  }
  return(private$.data)
})


LINCS$set("public", "metadata", function(rows) {  
  if(private$metadataIsStale) {
    private$loadMetadata()
    private$metadataIsStale = FALSE;
  }
  
  return(private$.metadata)
})


LINCS$set("private", "loadData", function(verbose=TRUE) {
  rows = private$rowset
  cols = private$colset
  print("Loading data from file into memory, this may take a few minutes")
  private$.data <- h5read(self$dataFile, "0/DATA/0/matrix", index=list(rows,cols))
  H5close()
  colnames(private$.data) <- self$colnamesByIndex(cols)
  if(private$geneids == 'entrez') {
    print("Summarizing data to entrez gene id, using mean() as summary method")
    ids = as.list(hgu133plus2ENTREZID)
    rownames <- unlist(ids[private$getProbes()[rows]])
    private$.data <- apply(private$.data, 2, function(x) { tapply(x, rownames, mean)})
    self$nrow = nrow(private$.data)
  } else {
    rownames(private$.data) <- private$getProbes()[rows]    
  }
})

LINCS$set("private", "loadMetadata", function() {
  cols = private$colset
  rows = self$metadataRows
  private$.metadata <- h5read(self$infoFile, "0/DATA/0/matrix", index=list(1:rows, cols))
  H5close()
  colnames(private$.metadata) <- self$colnamesByIndex(cols)
  rownames(private$.metadata) <- private$getMetaRowNames()
})


LINCS$set("public", "rowsByProbe", function(probes) {                      
  ix <- which(private$getProbes() %in% probes)
  ix
})



LINCS$set("public", "colnamesByIndex", function(cols) {
  result <- h5read(self$dataFile, "/0/META/COL/id", index=list(cols))
  result <- gsub(" ", "", result)
  result
})



LINCS$set("public", "colnamesById", function(ids) {
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
  result <- gsub(" ", "", result)
  result
})

LINCS$set("public", "l1000ProbeSets", function() {
  return(self$l1000ids()$probeset)
})

LINCS$set("public", "l1000Rows", function() {
  return(self$rowsByProbe(self$l1000ids()$probeset))
})

LINCS$set("public", "l1000ids", function() {
  #library(hgu133plus2.db)
  
  entrezids <- as.list(hgu133plus2ENTREZID)
  data(l1000Genes)
  entrezids_l1000 <- entrezids[which(entrezids %in% l1000Genes$EntrezGeneID)]
  ix <- match(entrezids_l1000, l1000Genes$EntrezGeneID)
  probesets_l1000 <- cbind(as.data.frame(names(entrezids_l1000), stringsAsFactors = FALSE), as.data.frame(unlist(entrezids_l1000), stringsAsFactors = FALSE), l1000Genes[ix,])[,-5]
  
  # add row numbers from data matrix in gctx file for convenience
  probesets_l1000 <- cbind(match(probesets_l1000[,1], private$getProbes()), probesets_l1000);
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


# info2hdf5
#
# Utility to convert the inst.info tab delimited instance metadata file from LINCS to 
# hdf5 format that directly mirrors the format of the expression data file, so by 
# default the resulting file is give the 'gctx' suffice rather than 'hdf'.
#
# The resulting file is 1/10th the size and much faster to work with due to random 
# access provided by hdf5.
#
# This only needs to be done once.
#
LINCS$set("public", "info2hdf5", function(infile="/mnt/lincs/inst.info", 
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
  
  private$h5WriteMatrix(outfile, t(info), "0/DATA/0/matrix")
})


# Workaround of the bug in rhdf5 (or possible libhdf5) that produces a segfaul when 
# writing large data sets (> ~20M datapoints) to hdf5 file.  
#
# name: name of the dataset in the hdffile, which already have been created
#       (i.e. by h5createDataset)
#
LINCS$set("private", "h5WriteMatrix", function(hdffile, data, name, colchunk=250000) {
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


# cellInfo: method to retreve cell line metadata via LINCS web API
#
# By default returs the normal cell lines, but your own query can be provided. 
# See echo http://api.lincscloud.org/a2/docs/cellinfo for details,
#
LINCS$set("public", "cellInfo", function(query=NA) {
  if(is.na(query)) {
    query <- '"cell_type":{"$regex":"[^-666|cancer]"}'
  }
  info <- getURL(paste('http://api.lincscloud.org/a2/cellinfo?q={', query, '}&user_key=lincsdemo', sep=""))
  return(private$JSON2df(info))
})

# usage: lincs$filterInstances('cell_id', 'MCF7')
# set reset = TRUE to search all columns.  Otherwise only currently selected columns are queried.
LINCS$set("public", "filterInstances", function(field, values, reset=FALSE) {
  if(reset) {
    self$setCols();     
    print(paste("Querying entire dataset (", self$ncol, " columns) for values of interest.", sep=""))
  } else {
    print(paste("Querying currently selected dataset (", self$ncol, " columns) for values of interest.", sep=""))    
  }
  private$loadMetadata()
  ix <- which(private$.metadata[field,] %in% values)
  self$setCols(private$colsset[ix])
  print(paste(length(ix), " columns identified and selected in data object."))
})


LINCS$set("public", "summarizeGenes", function(datamatrix, FUN=mean) {  
  eg <- unlist(as.list(hgu133plus2ENTREZID)[rownames(datamatrix)])
  apply(datamatrix, 2, function(x){tapply(x, eg, FUN)})
})

# fdaDrugs: method to retreve perturbagens that are predominently clinical drugs
#
# This is a very rough approximation, likely with low sensitivity but high specificity
LINCS$set("public", "fdaDrugs", function(query=NA) {
  first = TRUE
  for(i in 1:3) {
    data <- getURL(paste('http://api.lincscloud.org/a2/pertinfo?q={"pert_type":"trt_cp"}',
                         '&f={"pert_iname":1,"pert_id":1}',
                         '&l=10000&sk=', (i-1)*10000, '&user_key=lincsdemo', sep=""))
    #drugs <- JSON2df(data)
    #View(drugs)
    if(first) {
      drugs <- private$JSON2df(data)      
      first <- FALSE
    } else {
      drugs <- rbind(drugs, private$JSON2df(data))
    }
  }
  data(fdadrugs)
  ix <- which(drugs[,3] %in% fdadrugs)
  drugs[ix,]
})


# JSON2df: private utility method
#
# convert array of simple JSON objects to dataframe, ensuring that each column is accounted for 
# in each element of the JSON array.
#
# json: the data sructure, which must be an array of anonmyous JSON objects, like this:
#
# [
#   {
#     one: "some data",
#     two: "some more data"
#   },
#   {
#     one: "another bit of data",
#     two: "yet more",
#     three: "and this element has a third member!"
#   }
# ]
#              

LINCS$set("private", "JSON2df", function(json) {
  data <- fromJSON(json)
  res <- character()
  ids <- unique(unlist(lapply(data, names)))
  for(i in 1:length(data)) {
    row <- data[[i]]
    res <- rbind(res, row[ids])
  }
  colnames(res) <- gsub("_", "", ids)
  res[res %in% c('-666', 'NULL')] <- NA
  as.data.frame(res)
})

# grab the data
# raw_data <- getURL('http://api.lincscloud.org/a2/cellinfo?q={"cell_type":{"$regex":"[^-666|cancer]"}}&f={"cell_id":1,"cell_lineage":1,"cell_type":1,"gender":1}&user_key=lincsdemo')
# raw_data <- getURL('http://api.lincscloud.org/a2/cellinfo?q={"cell_type":{"$regex":"[^-666|cancer]"}}&user_key=lincsdemo')
# raw_data
# 
# unlist(fromJSON(raw_data))
# final_data <- do.call(rbind, fromJSON(raw_data))
# 
# final_data
