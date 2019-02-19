NCT <- function(data1, 
                data2, 
                it = 500, 
                estimation_method=c("EBICglasso", "association", "concentration", "IsingFit", "custom"), 
                paired=FALSE, 
                weighted=TRUE, 
                progressbar=TRUE, 
                test_edges=FALSE, 
                edges,
                nodes,
                custom_estimation_method, 
                centrality_global=c("ExpectedInfluence", "Strength"),
                test_centrality_node=FALSE,
                centrality_node=c("ExpectedInfluence", "Strength", "Betweenness", "Closeness"), 
                p_adjust_methods=c("none", "holm","fdr","BY")){
  
  if (progressbar == TRUE) 
    pb <- txtProgressBar(max = it, style = 3)
 
  ## Set up data structures
  x1 <- data.frame(data1)
  x2 <- data.frame(data2)
  nobs1 <- nrow(x1)
  nobs2 <- nrow(x2)
  dataall <- rbind(x1,x2)
  data.list <- list(x1,x2)
  b <- 1:(nobs1+nobs2)
  nvars <- ncol(x1)
  nedges <- nvars*(nvars-1)/2
  glstrinv.perm <- glstrinv.real <- nwinv.real <- nwinv.perm <- c()
  diffedges.perm <- matrix(0, it, nedges)
  diffedges.permtemp <- matrix(0, nvars, nvars)
  einv.perm.all <- array(0, dim = c(nvars, nvars, it))
  edges.pval.HBall <- matrix(0, nvars, nvars)
  edges.pvalmattemp <- matrix(0, nvars, nvars)
  cen.pvalmattemp <- matrix(0, nvars, 4)
  
  if(!is.null(colnames(data1))){
    nodenames <- colnames(data1)
  } else{
    nodenames <- 1:(dim(data1)[2])
  }
  
  ## Centrality
  if(missing(nodes)){
    nodes <- "all"
  }
  if(nodes=="all"){
    nnodes <- nvars
  } else {
    nodes<- unlist(nodes)
    nnodes <- length(nodes)
  } 
  centrality_global <- match.arg(centrality_global)
  centrality_node <- match.arg(centrality_node)
  gl.str <- switch(centrality_global,
                   Strength = function(x){
                     return(abs(sum(abs(x[upper.tri(x)]))))
                   } ,
                   ExpectedInfluence = function(x){
                     return(sum(x[upper.tri(x)]))
                   } )
  diffcen.perm <- matrix(NA, it, nnodes*length(centrality_node))

  ## Determine which function to use 
  if(match.arg(estimation_method)=="association"){
    fun <- cor
  } else if(match.arg(estimation_method)=="concentration"){
    fun <- function(x){
      ppcor::pcor(x)$estimate
    }
  } else if(match.arg(estimation_method)=="EBICglasso"){
    fun <- function(x){
      return(qgraph::EBICglasso(cor(x), n=dim(x)[1]))
    }
  } else if(match.arg(estimation_method)=="IsingFit"){
    fun <- function(x){
      return(IsingFit::IsingFit(x, progressbar=FALSE, plot=F)$weiadj)
    }
  } else if(match.arg(estimation_method)=="custom"){
    fun <- custom_estimation_method
  }
  
  ## Run the permutation
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
    if(test_centrality_node==TRUE){
      cen1 <- centrality_auto(nw1)$node.centrality
      cen2 <- centrality_auto(nw2)$node.centrality
      cen1$ExpectedInfluence <- networktools::expectedInf(nw1)$step1
      cen2$ExpectedInfluence <- networktools::expectedInf(nw2)$step1
      diffcen.real <- as.matrix(cen1) - as.matrix(cen2)
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
      diffedges.permtemp <- matrix(0, nvars, nvars)
      diffedges.permtemp[upper.tri(diffedges.permtemp, diag = FALSE)] <- diffedges.perm[i, ]
      diffedges.permtemp <- diffedges.permtemp + t(diffedges.permtemp)
      einv.perm.all[, , i] <- diffedges.permtemp
      nwinv.perm[i] <- max(diffedges.perm[i, ])
      if(test_centrality_node==TRUE){
        cen1permtemp <- centrality_auto(r1perm)$node.centrality
        cen2permtemp <- centrality_auto(r2perm)$node.centrality
        cen1permtemp$ExpectedInfluence <- networktools::expectedInf(r1perm)$step1
        cen2permtemp$ExpectedInfluence <- networktools::expectedInf(r2perm)$step1
        diffcen.permtemp <- as.matrix(cen1permtemp) - as.matrix(cen2permtemp)
        if(nodes=="all"){
          diffcen.perm[i,] <- melt(diffcen.permtemp[,centrality_node])$value
        } else {
          diffcen.perm[i,] <- melt(diffcen.permtemp[nodes,centrality_node])$value
        } 
      }
      if (progressbar == TRUE) 
        setTxtProgressBar(pb, i)
    }
    
    if (test_edges == TRUE) {
      if(missing(edges)){
        edges <- "all"
      }
      edges.pvaltemp <- colSums(diffedges.perm >= diffedges.realmat)/it
      edges.pvalmattemp[upper.tri(edges.pvalmattemp, diag = FALSE)] <- edges.pvaltemp
      edges.pvalmattemp <- edges.pvalmattemp + t(edges.pvalmattemp)
      rownames(edges.pvalmattemp) <- colnames(edges.pvalmattemp) <- nodenames
      rownames(diffedges.realoutput) <- colnames(diffedges.realoutput) <- nodenames
      rownames(einv.perm.all) <- colnames(einv.perm.all) <- nodenames
      ## is.character - purpose is to check if edges == "all"
      if (is.character(edges)) {
        if(p_adjust_methods[1]=="none"){
          ept.HBall <- edges.pvaltemp
        } else {ept.HBall <- p.adjust(edges.pvaltemp, method = p_adjust_methods[1])}
        edges.pval.HBall[upper.tri(edges.pval.HBall, 
                                   diag = FALSE)] <- ept.HBall
        edges.pval.HBall <- edges.pval.HBall + t(edges.pval.HBall)
        colnames(edges.pval.HBall) <- rownames(edges.pval.HBall) <- nodenames
        einv.pvals <- melt(edges.pval.HBall, na.rm = TRUE, 
                           value.name = "p-value")
        einv.perm <- einv.perm.all
        einv.real <- diffedges.realoutput
        edges.tested <- "all"
      }
      if (is.list(edges)) {
        einv.perm <- matrix(NA, it, length(edges))
        colnames(einv.perm) <- edges
        uncorrpvals <- einv.real <- vector()
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
        if(p_adjust_methods[1]=="none"){
          HBcorrpvals <- uncorrpvals
        } else {HBcorrpvals <- p.adjust(uncorrpvals, method = p_adjust_methods[1])}
        einv.pvals <- HBcorrpvals
        edges.tested <- colnames(einv.perm)
      }

    }
    if(test_centrality_node==TRUE){
      if(nodes=="all"){
        diffcen.real.vec <- melt(diffcen.real[,centrality_node])$value
        nnodes <- nvars
      } else {
        diffcen.real.vec <- melt(diffcen.real[nodes,centrality_node])$value
        nnodes <- length(nodes)
      } 
      diffcen.realmat <- matrix(diffcen.real.vec, it, nnodes, 
                                  byrow = TRUE)
      diffcen.pvaltemp <- colSums(abs(diffcen.perm) >= abs(diffcen.realmat))/it
      if(p_adjust_methods[1]=="none"){
        diffcen.HBall <- diffcen.pvaltemp
      } else {diffcen.HBall <- p.adjust(diffcen.pvaltemp, method = p_adjust_methods[1])}
      diffcen.pval <- matrix(diffcen.HBall, nnodes, length(centrality_node))
      colnames(diffcen.pval) <- centrality_node
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
                  nwinv.perm = nwinv.perm, method="permute", centrality_global=match.arg(centrality_global),
                  estimation_method=match.arg(estimation_method))
    if(test_edges == TRUE){
      res[["edges.tested"]] <- edges.tested
      res[["einv.real"]] <- einv.real
      res[["einv.pvals"]] <- einv.pvals
      res[["einv.perm"]] <- einv.perm
    }
    if(test_centrality_node==TRUE){
      res[["diffcen.real"]] <- matrix(diffcen.real.vec, nrow=nnodes,ncol=1, 
                                      dimnames=list(rownames(diffcen.pval), match.arg(centrality_node)))
      res[["diffcen.perm"]] <- diffcen.perm
      res[["diffcen.pval"]] <- diffcen.pval
    }
  class(res) <- "NCT"
  return(res)
  message("This is the development version of NCT. Association networks and global expected influence are used by default")
}

