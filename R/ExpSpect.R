#' ExpSpect class
#' 
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @description A class to perform expression based 'spectrometry'.  
#'
#' @format \code{\link{ExpSpect}} class generator
#'
#' @usage \code{exps = ExpSpect$new()}
#'
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
ExpSpect <- R6Class("ExpSpect",
  public = list(
    calcScores = function(exp, treated, untreated, lincs) {
      treated <- apply(exps(exp)[,treated], 1, mean)
      untreated <- apply(exps(exp)[,untreated], 1, mean)
      if(sum(rownames(treated) %in% lincs$getGeneIds) > 0.5 * nrow(treated)) {
        warning("Less than 50% of gene ids in expression set are present in lincs object.\nPossible gene id mismatch?\n")
      }
      exp_r <- treated/untreated
      ix <- match(lincs$getGeneIds(), rownames(treated))
      exp_r <- order(exp_r[ix], decreasing=TRUE)
      rankmatrix <- apply(lincs$data(), 2, order, decreasing=TRUE)
      cor <- apply(rankmatrix, 2, function(x) { return(cor.test(x, exp_r)$statistic) })
    }
  )
)

