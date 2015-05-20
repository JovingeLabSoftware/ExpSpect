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
      treated <- apply(exprs(exp)[,treated], 1, mean)
      untreated <- apply(exprs(exp)[,untreated], 1, mean)
      if(sum(names(treated) %in% lincs$getGeneIds()) < (0.5 * length(lincs$getGeneIds()))) {
        warning("Less than 50% of gene ids in expression set are present in lincs object.\nPossible gene id mismatch?\n")
      }
      exp_r <- treated/untreated
      ix <- match(lincs$getGeneIds(), rownames(treated))
      exp_r <- order(exp_r[ix], decreasing=TRUE)
      rankmatrix <- apply(lincs$data(), 2, order, decreasing=TRUE)
      cor <- apply(rankmatrix, 2, function(x) { return(cor.test(x, exp_r)$estimate) })
      names(cor) <- lincs$metadata()['pert_id',]
      cor
    }
  )
)

ExpSpect$set("public", "excos", function(treatment, condition, ngenes=c(100,100,100,100)) {
#  XSet=(UpInDisease∪DownInDisease)∩(ChangedByCompound)
  rank_tr <- order(treatment, decreasing=TRUE)
  rank_cnd <- order(condition, decreasing=TRUE)
  gene_tr <- names(treatment)[c(rank_tr[1:ngenes[1]], rank_tr[(length(treatment)-ngenes[2]):length(treatment)])]
  gene_cnd <- names(condition)[c(rank_cnd[1:ngenes[3]], rank_cnd[(length(condition)-ngenes[4]):length(condition)])]
  xset <- gene_tr[which(gene_tr %in% gene_cnd)]
  if(length(xset) == 0) {
    stop("No common genes found between gene sets")
  }
})

