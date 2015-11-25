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
    scores = NA,
    labels = NA,
    
    calcScores = function(x, 
                          treated, 
                          control, 
                          method = c('cmap', 'excos', 'excount'), 
                          returnScores = FALSE,
                          perm = FALSE,
                          parallel = FALSE,
                          ...) {
      scores <- private$.score(x, treated, control, method, perm, parallel, ...)

      self$scores <- scores
      self$labels <- colnames(x)
      if(returnScores) {
        return(scores)
      }      
     }
  ),
  private = list(
    .score = function(x, 
                      treated,
                      control, 
                      method, 
                      perm,
                      parallel,
                      ...) {
      cl <- call(paste('self[["', method, '"]](x=x, treated=treated, control=control, parallel=parallel, perm=perm, ...)', sep=""))  
      eval(parse(text=cl))
    },
    pb = NA,
    pb.init = function(max) {
      private$pb = txtProgressBar(max=max, style=3, width=30)
    },
    pb.step = function() {
      setTxtProgressBar(private$pb, getTxtProgressBar(private$pb) + 1)
    },
    pb.close = function() {
      close(private$pb)
    }
  )
)

#
# FUN must be a function that takes the following arguments:
#   x:        a matrix of reference perturbation data such as extracted from LINCS
#   treated:  a matrix of gene expression data for "treated" samples 
#   control:  a matrix of gene expression data for "control" samples 
#   ...:      optional additional arguments that will be passed through to FUN
#
ExpSpect$set("private", ".perm", function(x, treated, control, FUN, iterations=100, parallel=FALSE, sample=0.5, ...) {
  exc <- FUN(x, treated, control, ...)
  data <- cbind(treated, control)
  scores <- numeric();
  if(sample > 0.5) {
    warning("sample must be <= 0.5.  Using 0.5.")
    sample <- 0.5
  }
  smp <- floor(sample * ncol(data))
  if(parallel) {
    cat("Permuting in parallel (progress updates not available)...\n")
    require(parallel)
    f <- function(i) {
      s <- sample(1:ncol(data), ncol(data))
      FUN(x, 
          data[, s[1:smp]], 
          data[, s[(smp+1):(2*smp)]], 
          ...)
    }  
    scores <- mclapply(1:iterations, f, mc.cores=detectCores(), mc.preschedule=FALSE)
    print("...complete.")
    scores <- matrix(unlist(scores), nrow=iterations, byrow = TRUE)
  } else {
    cat("Permuting...")
    private$pb.init(iterations)
    for(i in 1:iterations) {
      private$pb.step();
      s <- sample(1:ncol(data), ncol(data))
      scores <- rbind(scores, FUN(x, 
                                  data[, s[1:smp]], 
                                  data[, s[(smp+1):(2*smp)]], 
                                  ...))
    } 
    private$pb.close()
  }
  exc <- apply(sweep(scores, 2, exc, '>'), 2, sum) / iterations
  exc <- -log2(exc)
  ix <- which(is.infinite(exc))
  if(length(ix)) {
    exc[ix] <- ceiling(max(exc[-ix]))    
  }
  exc
})


ExpSpect$set("public", "excount", function(x, 
                                           treated, 
                                           control, 
                                           n=100,
                                           perm = FALSE, 
                                           parallel=FALSE, 
                                           iterations=100, 
                                           sample=0.5) {
  
  treated <- treated[which(rownames(treated) %in% rownames(x)),]
  control <- control[which(rownames(control) %in% rownames(x)),]
  x <- apply(x, 2, rank)    
    
  if(perm) {
    exc <- private$.perm(x, 
                         treated, 
                         control, 
                         private$.excount, 
                         sample=sample, 
                         parallel=parallel,
                         iterations=iterations, 
                         n=n)
  } else {
    exc <- private$.excount(x, treated, control, n)
  }
  exc  
})

ExpSpect$set("private", ".excount", function(x, treated, control, n) {
  tr <- apply(treated, 1, mean)
  ct <- apply(control, 1, mean)
  y <- rank(tr / ct)
  
  up <- names(y)[which(y > (length(y)-n))]
  down <- names(y)[which(rank(y) < n)]  
  f <- function(x) {
    (sum(names(x[which(x > (length(x) - n))]) %in% up) + 
       sum(names(x[which(x < n)]) %in% down))  / (2*n)
  }
  exc <- apply(x, 2, f)    
  exc
})



ExpSpect$set("private", ".excos", function(x, treated, control, threshold) {
  tr <- apply(treated, 1, mean)
  ct <- apply(control, 1, mean)
  y <- log2(tr / ct)
  y <- (y-median(y))/sd(y)
  up <- names(y)[which(y > threshold)]
  up <- up[which(up %in% rownames(data))]
  down <- names(y)[which(y < -threshold)]
  down <- down[which(down %in% rownames(data))]
  if(length(up) < 2 || length(down) < 2) {
    warning("too few genes met threshold...returning NA")
    return(rep(NA, ncol(x)))
  }
  xset <- c(up, down)
  exc <- t(x[xset,]) %*% y[xset]
  exc[,1]
})

  
ExpSpect$set("public", "excos", function(x, 
                                           treated, 
                                           control, 
                                           threshold=1.96,
                                           perm = FALSE, 
                                           iterations=100, 
                                           parallel=FALSE,
                                           sample=0.5) {  

  
  treated <- treated[which(rownames(treated) %in% rownames(x)),]
  control <- control[which(rownames(control) %in% rownames(x)),]
  
  if(perm) {
    exc <- private$.perm(x, 
                         treated, 
                         control, 
                         private$.excos, 
                         sample=sample, 
                         parallel=parallel,
                         iterations=iterations,
                         threshold=threshold)
  } else {
    exc <- private$.excos(x, treated, control)
  }
  exc
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

ExpSpect$set("public", "cmap", function(x, 
                                         treated, 
                                         control, 
                                         threshold=1.96,
                                         perm = FALSE, 
                                         iterations=100, 
                                         sample=0.5,
                                         parallel=FALSE, 
                                         rescale=TRUE) {  
  
  
  treated <- treated[which(rownames(treated) %in% rownames(x)),]
  control <- control[which(rownames(control) %in% rownames(x)),]
  
  if(perm) {
    exc <- private$.perm(x, 
                         treated, 
                         control, 
                         private$.cmap, 
                         sample=sample, 
                         iterations=iterations,
                         threshold=threshold,
                         parallel=parallel,
                         rescale=rescale)
  } else {
    exc <- private$.cmap(x, treated, control, threshold, rescale=rescale)
  }
  exc
})

ExpSpect$set("private", ".cmap", function(x, treated, control, threshold=1.96, normalize=TRUE, rescale=TRUE) {
  names <- rownames(treated)
  
  tr <- apply(treated, 1, mean)
  ct <- apply(control, 1, mean)

  y <- tr/ct
  y <- (y-mean(y))/sd(y)
  
  
  up <- names[which(y > threshold)]

  down <- names[which(y < -threshold)]
  up <- up[which(up %in% rownames(x))]
  down <- down[which(down %in% rownames(x))]
  if(length(up) == 0 || length(down) == 0) {
    return(rep(0, ncol(x)))
    warn("No up and/or down regulated genes identified.  Returning 0 as cmap scores.")
  }    
  
  s_i <- apply(x, 2, function(xx) {
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
  return(s_i)  
})

ExpSpect$set("public", "fetchATC", function(drugs) {
  codes <- read.table("http://www.ebi.ac.uk/Rebholz-srv/atc/public/ontologies/atc.owl", sep="\n", as.is=TRUE, header=FALSE)
  data(atc)
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

ExpSpect$set("public", "plot", function(squelch = 7, highlight = NA, labels=10, aggregate=FALSE, FUN=max, palette=rainbow, ...) {
  drugs <- tolower(gsub(" .*", "", names(self$scores)))
  
  if(aggregate) {
    data <- tapply(self$scores, self$labels, FUN)
    ix <- sort(data, index.return=TRUE)$ix
    data <- data[ix]
  } else {
    data <- self$scores
    names(data) <- self$labels
  }

  col <- rep(NA, length(data))
  if(length(highlight) == 1) {
    col[grep(highlight, names(data))] <- 'red'    
  } else  {
    col[which(names(data) %in% highlight)] <- 'red'        
  }
  
  par(mgp=c(3,.5,0), mar=c(8,3,3,0))
  barplot(data**squelch, col='orange', border='orange',  space = 0, width=1, names=NA, xlab="drug", ylab="ExpSpect Score", axes=FALSE, ...) 
  par(new=TRUE)
  barplot(data**squelch, col=col, border=col,  names=NA, xlab="drug", width=1, space=0, ylab="ExpSpect Score", axes=FALSE, ...) 
  
  ix <- which(data >= sort(data, decreasing = TRUE)[labels])
  ax <- ix
  spacing <- ix[-length(ix)] - ix[2:length(ix)] > -10
  ax[which(spacing)] <- ax[which(spacing)] - (0.01 * length(data))
  text(x = ax, par("usr")[3], labels = names(data)[ix], adj=c(0,0), srt = 290, cex=0.6, xpd = TRUE)
  axis(2, cex.axis=0.7)
})

# mmatch  -  find all ocurrences of the first argument in the second 
# returns a dataframe with columns comprised of elements of x and indices of table
#
mmatch <- function(x, table) {
  f <- function(a, b) { 
    if(a %in% b) { 
      data.frame(element=a, index=which(b == a), stringsAsFactors=FALSE) 
    } else { 
      data.frame(element=a, index=NA)
    }
  }
  do.call(rbind, lapply(x, f, table))
}


