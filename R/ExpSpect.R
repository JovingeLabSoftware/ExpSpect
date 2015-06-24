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
    calcScores = function(data, treated, untreated, method=c('cmap', 'excos', 'excount'), ...) {
      treated <- apply(treated, 1, mean)
      untreated <- apply(untreated, 1, mean)
      treated <- (treated - mean(treated))/sd(treated)
      untreated <- (untreated - mean(untreated))/sd(untreated)

      exp_r <- treated-untreated
      exp_r <- (exp_r - mean(exp_r)) / sd(exp_r)
      cl <- call(paste('self[["', method, '"]](x=data, y=exp_r, ...)', sep=""))
      eval(parse(text=cl))
     },
     permScores = function(x, y, method=c('cmap', 'excos', 'excount'), sample=0.2, i=100) {
       scores <- numeric();
       if(sample > 0.5) {
         warning("sample must be <= 0.5.  Using 0.5.")
       }
       smp <- floor(sample * ncol(y))
       for(i in 1:i) {
        print(i)
        s <- sample(1:ncol(y), ncol(y))
        treated <- s[1:smp]
        untreated <- s[(smp+1):(2*smp)]
        data <- data[sample(1:nrow(data), nrow(data)),]
        scores <- rbind(scores, self$calcScores(x, data[,treated], data[,untreated], method)$scores)
      } 
      scores
     }
    
  )
)

ExpSpect$set("public", "excount", function(x, y, n = 100
                                           ) {
  # coerce to numeric just in case
  y <- y[which(names(y) %in% rownames(x))]
  dm <- matrix(as.numeric(x), ncol=ncol(x))
  colnames(dm) <- colnames(x)
  rownames(dm) <- rownames(x)
  data <- dm
  y.nm <- names(y)
  y <- as.numeric(y)
  names(y) <- y.nm
  
  #if(normalize) {
  #  y <- (y-mean(y))/sd(y)
  #  data <- (data-mean(data))/sd(data)
  #}  
  # identify signature genes
  up <- names(y)[which(rank(y) > (length(y)-n))]
  up <- up[which(up %in% rownames(data))]
  down <- names(y)[which(rank(y) < n)]
  down <- down[which(down %in% rownames(data))]
  
  # calculate the 'extreme count'
  f <- function(x) {
    x <- rank(x)
    (sum(names(x[which(x > (length(x) - n))]) %in% up) + 
    sum(names(x[which(x < n)]) %in% down))  / (2*n)
  }
  exc <- apply(x, 2, f)
  return(list(scores=as.vector(exc), up=up, down=down))  
})

  

ExpSpect$set("public", "excos", function(x, y, threshold=1.96) {
  # coerce to numeric just in case
  dm <- matrix(as.numeric(x), ncol=ncol(x))
  colnames(dm) <- colnames(x)
  rownames(dm) <- rownames(x)
  data <- dm
  y.nm <- names(y)
  y <- as.numeric(y)
  names(y) <- y.nm
  
  #if(normalize) {
  #  y <- (y-mean(y))/sd(y)
  #  data <- (data-mean(data))/sd(data)
  #}  
  # identify signature genes
  up <- names(y)[which(y > threshold)]
  up <- up[which(up %in% rownames(data))]
  down <- names(y)[which(y < -threshold)]
  down <- down[which(down %in% rownames(data))]
  # ensure that we only use gene ids common to both data and y 
  xset <- c(up, down)
  if(length(xset) < 2) {
    stop("No common gene ids between data matrix and UP and DOWN regulated genes in y")
  }
  
  # calculate the 'extreme cosine'
  exc <- t(data[xset,]) %*% y[xset]
  return(list(scores=as.vector(exc), up=up, down=down))

})

ExpSpect$set("public", "threshold", function(x, limit=0.01, maxit=1000000) {
  conv <- FALSE
  m <- mean(x)
  c <- 0
  while(!conv) {
    c <- c+1
    m.p <- m
    m1 <- mean(x[which(x<m)])
    m2 <- mean(x[which(x>m)])
    m <- (m1 + m2)/2
    if(abs(m-m.p) < limit){
      conv=TRUE
    }
    if(c > maxit) {
      warning("threshold alogorithm did not converge...returning best guess")
      conv=TRUE
    }
  }
  return(list(m=m, m1=m1, m2=m2))
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

  up <- up[which(up %in% rownames(x))]
  down <- down[which(down %in% rownames(x))]
  if(length(up) == 0 || length(down) == 0) {
    stop("No genes in up and/or down list found in rownames of data matrix.  Gene id mismatch?")
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
    
    V_d <- which(xx %in% down)
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
  
  return(list(scores=s_i, up=up, down=down))  
})

ExpSpect$set("public", "fetchATC", function(drugs) {
  codes <- read.table("http://www.ebi.ac.uk/Rebholz-srv/atc/public/ontologies/atc.owl", sep="\n", as.is=TRUE, header=FALSE)
  res <- rep("", 2)
  for(l in codes[,1]) {
    if(length(grep("owl:Class.*/atc/.*>", l))) {
      c <- gsub("<owl:Class.*/atc/(.*)>", "\\1", l)
      c <- gsub("\\s", "", c)
    }
    if(length(grep("rdfs:label.*string>.*<", l))) {
      d <- gsub('.*string\\"*>(.*?)<\\/.*', "\\1", l)
      res <- rbind(res, c(c,tolower(d)))
    }
  }
  entry <- rep("", 7)
  atc <- entry
  codes <- character()
  for(i in 1:nrow(res)) {
    ontology[nchar(res[i,1])] <- res[i,2]
    if(nchar(res[i,1]) == 7) {
      atc <- rbind(atc, ontology)
      codes <- c(codes, res[i,1])
    }
  }
  atc <- atc[-1,]
  atc <- atc[,-c(2,6)]
  rownames(atc) <- codes
  colnames(atc) <- c("main", "sub1", "sub2", "sub3", "substance") 
  
  atc
})


