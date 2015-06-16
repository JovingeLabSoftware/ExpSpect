Test <- R6Class("Test",
                    public = list(
                      testEval = function(...) {
                        cl = call(paste0("self[['", 'evalMe', "']]", "('ok', ...)"))
                        eval(parse(text=cl))
                      },
                      evalMe = function(msg) {
                        print(msg)
                      }                        
                    )             
)

excos <- function(data, y, threshold=1.96, normalize=TRUE) {
  dm <- matrix(as.numeric(data), ncol=ncol(data))
  colnames(dm) <- colnames(data)
  rownames(dm) <- rownames(data)
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
  down <- names(y)[which(y < -threshold)]
  # ensure that we only use gene ids common to both data and y 
  xset <- rownames(data)[which(rownames(data) %in% c(up, down))]
  if(length(xset) == 0) {
    stop("No common gene ids between data matrix and UP and DOWN regulated genes in y")
  }
  
  # calculate the 'extreme cosine'
  exc <- t(data[xset,]) %*% y[xset]
  print(str(exc))
  return(as.vector(exc))
}


cmap <- function(x, y, threshold=1.96, normalize, rescale) {
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
}