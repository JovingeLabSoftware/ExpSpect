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
#'  \item{\code{calcScores(data, treated, untreated, sigs, method=c('ks', 'excos'), normalize=TRUE, ...)}}{Given a gene expression matrix \code{data}, computes fold change based on 
#'       columns given in vectors \code{treated} and \code{untreated}, and then computes enrichment scores against each column of 
#'       matrix \code{sigs}.  The \code{sigs} matrix is simply a matrix of gene expression data, one column per pertubagen which 
#'       can be generated in any number of ways--assembled from GEO, downloaded from CMAP, extracted from LINCS, etc.  See 
#'       Vignette for an example.  Data will be renormalized after fold change calculation prior to processing by \code{method} if \code{normalize=TRUE}.
#'       Extra-parameters (\code{...}) will be sent to selected function (see below for details)}
#'  \item{\code{cmap(data, y=NULL, up, down)}}{Given a gene expression matrix \code{data}, computes the connectivity 
#'       score based on either a vector of expression data (\code{y}), or explicit sets of \code{up} and \code{down} 
#'       regulated genes pertinent to a biological state or perturbation of interest.  Returns the list of connectivity 
#'       scores for each column of \code{data}. If providing a vector of expression values, up and down genes are s
#'       selected based on the provided \code{threshold}, after optional normalization (\code{normalize = TRUE}) step.}
#' }
#' 
#'
#' @keywords data
#' 
ExpSpect <- R6Class("ExpSpect",
  public = list(
    calcScores = function(data, treated, untreated, method=c('cmap', 'excos'), normalize=TRUE, ...) {
      
      if(normalize) {
        data <- (data - mean(data)) / sd(data)
        sigs <- (data - mean(data)) / sd(data)
      }
      treated <- apply(treated, 1, mean)
      untreated <- apply(untreated, 1, mean)
      exp_r <- treated/untreated
      if(normalize) {
        exp_r <- (exp_r - mean(exp_r)) / sd(exp_r)
      }
       cl <- call(paste('self[["', method, '"]](x=data, y=exp_r, normalize=normalize, ...)', sep=""))
       eval(parse(text=cl))
     }
  )
)

  

ExpSpect$set("public", "excos", function(x, y, threshold=1.96, normalize=TRUE) {
  # coerce to numeric just in case
  dm <- matrix(as.numeric(x), ncol=ncol(x))
  colnames(dm) <- colnames(x)
  rownames(dm) <- rownames(x)
  data <- dm
  y.nm <- names(y)
  y <- as.numeric(y)
  names(y) <- y.nm
  
  if(normalize) {
    y <- (y-mean(y))/sd(y)
    data <- (data-mean(data))/sd(data)
  }  
  # identify signature genes
  up <- names(y)[which(y > threshold)]
  up <- up[which(up %in% rownames(data))]
  down <- names(y)[which(y < -threshold)]
  down <- down[which(down %in% rownames(data))]
  # ensure that we only use gene ids common to both data and y 
  xset <- c(up, down)
  if(length(xset) == 0) {
    stop("No common gene ids between data matrix and UP and DOWN regulated genes in y")
  }
  
  # calculate the 'extreme cosine'
  exc <- t(data[xset,]) %*% y[xset]
  return(list(excos=as.vector(exc), up=up, down=down))

})

  
ExpSpect$set("public", "cmap", function(x, y, threshold=1.96, normalize=TRUE, rescale=TRUE) {
  if(is.list(y)) {
    up <- y$up;
    down <- y$down;
  } else {
    # sometimes imported data comes accross as character instead of numeric, will ensure numeric...
    names <- names(y)
    y <- as.numeric(y)
    if(normalize) {
      y <- (y-mean(y))/sd(y)
    }
    up <- names[which(y > threshold)]
    down <- names[which(y < -threshold)]
    if(length(up) == 0 || length(down) == 0) {
      stop("No genes met threshold criteria.")
    }    
  }
  
  s_i <- apply(x, 2, function(xx) {
    # again, coercing to numeric just in case
    ix <- sort(as.numeric(xx), index.return=TRUE)$ix
    xx <- names(xx)[ix]
    t_u <- length(up)
    t_d <- length(down)
    
    j_u <- 1:t_u
    j_d <- 1:t_d
    
    n <- length(xx)
    V_u <- sort(which(xx %in% up))
    a_u  = max(j_u/t_u - V_u / n)
    b_u  = max(V_u / n - (j_u - 1) / t_u)
    if (a_u > b_u) {
      ks_u = a_u
    } else {
      ks_u = -b_u
    }
    
    V_d <- sort(which(xx %in% down))
    a_d  = max(j_d/t_d - V_d / n)
    b_d  = max(V_d / n - (j_d - 1) / t_d)
    if (a_d > b_d) {
      ks_d = a_d
    } else {
      ks_d = -b_d
    }
    (ks_u - ks_d) * ((sign(ks_u) * sign(ks_d)) < 0)
  })
  
  if(rescale) {
    # scale to -1..1 like the original CMAP paper...
    p <- max(s_i)
    q <- min(s_i)
    
    s_i[s_i > 0] <- s_i[s_i > 0] / p
    s_i[s_i < 0] <- -(s_i[s_i < 0] / q)    
  }
  
  return(s_i)
})


