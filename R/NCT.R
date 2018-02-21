NCT <- function(data1, data2, it = 100, 
                          default=c("association", "concentration", "EBICglasso", "IsingFit", "custom"), 
                          paired=FALSE, weighted=TRUE, progressbar=TRUE, test.edges=FALSE, edges,
                          custom_func, AND=TRUE, global=c("expectedInf", "strength"),test.centralities=FALSE,
                          centralities=c("Strength", "ExpectedInfluence", "Betweenness", "Closeness"), nodes,
                          adjust=c("holm","fdr","BY", "none")){
  
  global <- match.arg(global)
  gl.str <- switch(global,
                   strength = function(x){
                     return(abs(sum(abs(x[upper.tri(x)]))))
                   } ,
                   expectedInf = function(x){
                     return(sum(x[upper.tri(x)]))
                   } )
  #if(match.arg(global.strength) == "aboslute_value"){
   # gl.str <- function(x){
    #  return(abs(sum(abs(x[upper.tri(x)]))))
  #  } 
  #} else if (match.arg(global.strength) == "raw") {
   # gl.str <- function(x){
    #  return(sum(x[upper.tri(x)]))
    #} 
  #}
  
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
  
  ## Centrality
  if(missing(nodes)){
    nodes <- "all"
  }
  if(is.list(nodes)){
    nodes<- unlist(nodes)
  }
  if(nodes=="all"){
    nnodes <- nvars
  } else {
    nnodes <- length(nodes)
  } 
  diffcen.perm <- matrix(0, it, nnodes*length(centralities))
  
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
  cen.pvalmattemp <- matrix(0, nvars, 4)
  

  
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
    
    glstrinv.real <- gl.str(nw1) - gl.str(nw2)
    glstrinv.sep <- c(gl.str(nw1), gl.str(nw2))
    diffedges.real <- abs(nw1 - nw2)[upper.tri(abs(nw1 - 
                                                     nw2))]
    diffedges.realmat <- matrix(diffedges.real, it, nedges, 
                                byrow = TRUE)
    diffedges.realoutput <- abs(nw1 - nw2)
    nwinv.real <- max(diffedges.real)
    if(test.centralities==TRUE){
      cen1 <- centrality_auto(nw1)$node.centrality
      cen2 <- centrality_auto(nw2)$node.centrality
      cen1$ExpectedInfluence <- networktools::expectedInf(nw1)$step1
      cen2$ExpectedInfluence <- networktools::expectedInf(nw2)$step1
      diffcen.real <- abs(as.matrix(cen1) - as.matrix(cen2))
    }
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
      glstrinv.perm[i] <- gl.str(r1perm) - gl.str(r2perm) 
      diffedges.perm[i, ] <- abs(r1perm - r2perm)[upper.tri(abs(r1perm - 
                                                                  r2perm))]
      diffedges.permtemp[upper.tri(diffedges.permtemp, 
                                   diag = FALSE)] <- diffedges.perm[i, ]
      diffedges.permtemp <- diffedges.permtemp + t(diffedges.permtemp)
      einv.perm.all[, , i] <- diffedges.permtemp
      nwinv.perm[i] <- max(diffedges.perm[i, ])
      if(test.centralities==TRUE){
        cen1permtemp <- centrality_auto(r1perm)$node.centrality
        cen2permtemp <- centrality_auto(r2perm)$node.centrality
        cen1permtemp$ExpectedInfluence <- networktools::expectedInf(r1perm)$step1
        cen2permtemp$ExpectedInfluence <- networktools::expectedInf(r2perm)$step1
        diffcen.permtemp <- abs(as.matrix(cen1permtemp) - as.matrix(cen2permtemp))
        if(nodes=="all"){
          diffcen.perm[i,] <- melt(diffcen.permtemp[,centralities])$value
        } else {
          diffcen.perm[i,] <- melt(diffcen.permtemp[nodes,centralities])$value
        } 
      }
      if (progressbar == TRUE) 
        setTxtProgressBar(pb, i)
    }
    if (test.edges == TRUE) {
      edges.pvaltemp <- colSums(diffedges.perm >= diffedges.realmat)/it
      edges.pvalmattemp[upper.tri(edges.pvalmattemp, diag = FALSE)] <- edges.pvaltemp
      edges.pvalmattemp <- edges.pvalmattemp + t(edges.pvalmattemp)
      if (is.character(edges)) {
        if(match.arg(adjust)=="none"){
          ept.HBall <- edges.pvaltemp
        } else {ept.HBall <- p.adjust(edges.pvaltemp, method = martch.arg(adjust))}
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
        if(match.arg(adjust)=="none"){
          HBcorrpvals <- uncorrpvals
        } else {HBcorrpvals <- p.adjust(uncorrpvals, method = match.arg(adjust))}
        einv.pvals <- HBcorrpvals
      }
      edges.tested <- colnames(einv.perm)

    }
    if(test.centralities==TRUE){
      if(nodes=="all"){
        diffcen.real.vec <- melt(diffcen.real[,centralities])$value
        nnodes <- nvars
      } else {
        diffcen.real.vec <- melt(diffcen.real[nodes,centralities])$value
        nnodes <- length(nodes)
      } 
      diffcen.pvaltemp <- colSums(diffcen.perm >= diffcen.real.vec)/it
      if(match.arg(adjust)=="none"){
        diffcen.HBall <- diffcen.pvaltemp
      } else {diffcen.HBall <- p.adjust(diffcen.pvaltemp, method = match.arg(adjust))}
      diffcen.pval <- matrix(diffcen.HBall, nnodes, length(centralities))
      colnames(diffcen.pval) <- centralities
      if(nodes=="all"){
        rownames(diffcen.pval) <- rownames(diffcen.real)
      } else {rownames(diffcen.pval) <- nodes}
    }
    if (progressbar == TRUE) 
      close(pb)
      res <- list(glstrinv.real = glstrinv.real, glstrinv.sep = glstrinv.sep, 
                  glstrinv.pval = sum(abs(glstrinv.perm) >= abs(glstrinv.real))/it, 
                  glstrinv.perm = glstrinv.perm, nwinv.real = nwinv.real, 
                  nwinv.pval = sum(nwinv.perm >= nwinv.real)/it, 
                  nwinv.perm = nwinv.perm, method="permute", global=match.arg(global),
                  default=match.arg(default))
    if(test.edges == TRUE){
      res[["edges.tested"]] <- edges.tested
      res[["einv.real"]] <- einv.real
      res[["einv.pvals"]] <- einv.pvals
      res[["einv.perm"]] <- einv.perm
    }
    if(test.centralities==TRUE){
      res[["diffcen.real"]] <- diffcen.real
      res[["diffcen.perm"]] <- diffcen.perm
      res[["diffcen.pval"]] <- diffcen.pval
    }
  class(res) <- "NCT"
  return(res)
  message("This is the development version of NCT. Association networks and global expected influence are used by default")
}

