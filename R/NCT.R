NCT <- function(data1, data2, it = 100, 
                          default=c("association", "concentration", "EBICglasso", "IsingFit", "custom"), 
                          paired=FALSE, weighted=TRUE, progressbar=TRUE, test.edges=FALSE, edges,
                          custom_func, AND=TRUE, global.strength=c("raw", "absolute_value")){
  
  if(match.arg(global.strength) == "aboslute_value"){
    gl.str <- function(x){
      return(abs(sum(abs(x[upper.tri(x)]))))
    } 
  } else if (match.arg(global.strength) == "raw") {
    gl.str <- function(x){
      return(sum(x[upper.tri(x)]))
    } 
  }
  
  ## From boot
  x1 <- data.frame(data1)
  x2 <- data.frame(data2)
  nobs1 <- nrow(x1)
  nobs2 <- nrow(x2)
  dataall <- rbind(x1,x2)
  data.list <- list(x1,x2)
  b <- 1:(nobs1+nobs2)
  nvars <- ncol(x1)
  nedges <- nvars*(nvars-1)/2
  
  ## From original
  if (progressbar == TRUE) 
    pb <- txtProgressBar(max = it, style = 3)
  x1 <- data1
  x2 <- data2
  nobs1 <- nrow(x1)
  nobs2 <- nrow(x2)
  dataall <- rbind(x1, x2)
  data.list <- list(x1, x2)
  b <- 1:(nobs1 + nobs2)
  nvars <- ncol(x1)
  nedges <- nvars * (nvars - 1)/2
  glstrinv.perm <- glstrinv.real <- nwinv.real <- nwinv.perm <- c()
  diffedges.perm <- matrix(0, it, nedges)
  diffedges.permtemp <- matrix(0, nvars, nvars)
  einv.perm.all <- array(NA, dim = c(nvars, nvars, it))
  edges.pval.HBall <- matrix(NA, nvars, nvars)
  edges.pvalmattemp <- matrix(0, nvars, nvars)
  

  
  ##### Determine which function to use ####
  if(match.arg(default)=="association"){
    fun <- cor
  } else if(match.arg(default)=="concentration"){
    fun <- function(x){
      ppcor::pcor(x)$estimate
    }
  } else if(match.arg(default)=="EBICglasso"){
    fun <- function(x){
      return(qgraph::EBICglasso(cor(x), n=dim(x)[1]))
    }
  } else if(match.arg(default)=="IsingFit"){
    fun <- function(x){
      return(IsingFit::IsingFit(x, AND=AND, progressbar=FALSE, plot=F)$weiadj)
    }
  } else if(match.arg(default)=="custom"){
    fun <- custom_func
  }
  
  
  ##### Run the permutation #####

    nw1 <- fun(x1)
    nw2 <- fun(x2)
    if (weighted == FALSE) {
      nw1 = (nw1 != 0) * 1
      nw2 = (nw2 != 0) * 1
    }
    glstrinv.real <- abs(sum(abs(nw1[upper.tri(nw1)])) - 
                           sum(abs(nw2[upper.tri(nw2)])))
    glstrinv.sep <- c(sum(abs(nw1[upper.tri(nw1)])), sum(abs(nw2[upper.tri(nw2)])))
    diffedges.real <- abs(nw1 - nw2)[upper.tri(abs(nw1 - 
                                                     nw2))]
    diffedges.realmat <- matrix(diffedges.real, it, nedges, 
                                byrow = TRUE)
    diffedges.realoutput <- abs(nw1 - nw2)
    nwinv.real <- max(diffedges.real)
    for (i in 1:it) {
      if (paired == FALSE) {
        s <- sample(1:(nobs1 + nobs2), nobs1, replace = FALSE)
        x1perm <- dataall[s, ]
        x2perm <- dataall[b[-s], ]
        r1perm <- fun(x1perm)
        r2perm <- fun(x2perm)
        if (weighted == FALSE) {
          r1perm = (r1perm != 0) * 1
          r2perm = (r2perm != 0) * 1
        }
      }
      if (paired == TRUE) {
        s <- sample(c(1, 2), nobs1, replace = TRUE)
        x1perm <- x1[s == 1, ]
        x1perm <- rbind(x1perm, x2[s == 2, ])
        x2perm <- x2[s == 1, ]
        x2perm <- rbind(x2perm, x1[s == 2, ])
        r1perm <- fun(x1perm)
        r2perm <- fun(x2perm)
        if (weighted == FALSE) {
          r1perm = (r1perm != 0) * 1
          r2perm = (r2perm != 0) * 1
        }
      }
      glstrinv.perm[i] <- abs(sum(abs(r1perm[upper.tri(r1perm)])) - 
                                sum(abs(r2perm[upper.tri(r2perm)])))
      diffedges.perm[i, ] <- abs(r1perm - r2perm)[upper.tri(abs(r1perm - 
                                                                  r2perm))]
      diffedges.permtemp[upper.tri(diffedges.permtemp, 
                                   diag = FALSE)] <- diffedges.perm[i, ]
      diffedges.permtemp <- diffedges.permtemp + t(diffedges.permtemp)
      einv.perm.all[, , i] <- diffedges.permtemp
      nwinv.perm[i] <- max(diffedges.perm[i, ])
      if (progressbar == TRUE) 
        setTxtProgressBar(pb, i)
    }
    if (test.edges == TRUE) {
      edges.pvaltemp <- colSums(diffedges.perm >= diffedges.realmat)/it
      edges.pvalmattemp[upper.tri(edges.pvalmattemp, diag = FALSE)] <- edges.pvaltemp
      edges.pvalmattemp <- edges.pvalmattemp + t(edges.pvalmattemp)
      if (is.character(edges)) {
        ept.HBall <- p.adjust(edges.pvaltemp, method = "holm")
        edges.pval.HBall[upper.tri(edges.pval.HBall, 
                                   diag = FALSE)] <- ept.HBall
        rownames(edges.pval.HBall) <- colnames(edges.pval.HBall) <- colnames(data1)
        einv.pvals <- melt(edges.pval.HBall, na.rm = TRUE, 
                           value.name = "p-value")
        einv.perm <- einv.perm.all
        einv.real <- diffedges.realoutput
      }
      if (is.list(edges)) {
        einv.perm <- matrix(NA, it, length(edges))
        colnames(einv.perm) <- edges
        uncorrpvals <- einv.real <- c()
        for (j in 1:length(edges)) {
          uncorrpvals[j] <- edges.pvalmattemp[edges[[j]][1], 
                                              edges[[j]][2]]
          einv.real[j] <- diffedges.realoutput[edges[[j]][1], 
                                               edges[[j]][2]]
          for (l in 1:it) {
            einv.perm[l, j] <- einv.perm.all[, , l][edges[[j]][1], 
                                                    edges[[j]][2]]
          }
        }
        HBcorrpvals <- p.adjust(uncorrpvals, method = "holm")
        einv.pvals <- HBcorrpvals
      }
      edges.tested <- colnames(einv.perm)
      res <- list(glstrinv.real = glstrinv.real, glstrinv.sep = glstrinv.sep, 
                  glstrinv.pval = sum(glstrinv.perm >= glstrinv.real)/it, 
                  glstrinv.perm = glstrinv.perm, nwinv.real = nwinv.real, 
                  nwinv.pval = sum(nwinv.perm >= nwinv.real)/it, 
                  nwinv.perm = nwinv.perm, edges.tested = edges.tested, 
                  einv.real = einv.real, einv.pvals = einv.pvals, 
                  einv.perm = einv.perm, method="permute", global.strength=match.arg(global.strength),
                  default=match.arg(default))
    }
    if (progressbar == TRUE) 
      close(pb)
    if (test.edges == FALSE) {
      res <- list(glstrinv.real = glstrinv.real, glstrinv.sep = glstrinv.sep, 
                  glstrinv.pval = sum(glstrinv.perm >= glstrinv.real)/it, 
                  glstrinv.perm = glstrinv.perm, nwinv.real = nwinv.real, 
                  nwinv.pval = sum(nwinv.perm >= nwinv.real)/it, 
                  nwinv.perm = nwinv.perm, method="permute", global.strength=match.arg(global.strength),
                  default=match.arg(default))
    }
  
  class(res) <- "NCT"
  return(res)
  message("This is the development version of NCT. Association networks and sum connectivity are used by default")
}

