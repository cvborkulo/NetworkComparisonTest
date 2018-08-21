NCT_bootstrap <- function(data1, 
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
                          bootcut=c("none", "cutEqual"),
                          include_t=FALSE){
  centrality_node <- match.arg(centrality_node)
  if (progressbar==TRUE) pb <- txtProgressBar(max=it, style = 3)
  x1 <- data.frame(data1)
  x2 <- data.frame(data2)
  nobs1 <- nrow(x1)
  nobs2 <- nrow(x2)
  dataall <- rbind(x1,x2)
  data.list <- list(x1,x2)
  b <- 1:(nobs1+nobs2)
  ncolall <- 1:ncol(x1) ## This stores a vector for the length of the original number of columns in the data (before adding "group", below)
  nvars <- ncol(x1)
  nedges <- nvars*(nvars-1)/2
  if(match.arg(centrality_global) == "Strength"){
    gl.str <- function(x){
      return(abs(sum(abs(x[upper.tri(x)]))))
    } 
    warning("Global strength centrality uses absolute value, which can be problematic when bootstrapping")
  } else if (match.arg(centrality_global) == "ExpectedInfluence") {
    gl.str <- function(x){
      return(sum(x[upper.tri(x)]))
    } 
  }

  ## The next 7 lines generate an empty edgelist  
  ## We will store edge invariances in edgelist format, for ease with CIs
  if(missing(edges)){edges <- "all"}
  edgelistVec <- vector()
  edgelist <- matrix(NA,nvars,nvars) 
  edgelist[upper.tri(edgelist, diag=FALSE)] <- rep(1, nedges)
  rownames(edgelist) <- colnames(edgelist) <- colnames(data1)
  edgelist <- reshape2::melt(edgelist, na.rm=TRUE, value.name= "Edge Strength")
  for(i in 1:nedges) {
    edgelistVec[i] <- paste(edgelist[i, 1], edgelist[i, 2], sep = " -- ")
  }
  if(is.list(edges)){
    reduced_edgelistVec <- vector()
      for(i in 1:length(edges)){
        reduced_edgelistVec[i] <- paste(edges[[i]][1], edges[[i]][2], sep = " -- ")
      }
  } else{
    reduced_edgelistVec <- paste(edges[1], edges[2], sep = " -- ")
  }
  
  if(missing(nodes)){
    nodes <- "all"
  }
  if(is.list(nodes)){
    nodes<- unlist(nodes)
  }
  if(nodes=="all"){
    nodes <- colnames(data1)
    nnodes <- nvars
  } else {
    nnodes <- length(nodes)
  } 
  
  glstrinv <- list()
  einv <- list()
  diffcen <- list()
  
  glstrinv$t <- numeric(0) ## Empty vector to insert bootstrapped values
  einv$t <- data.frame(matrix(NA, it, length(edgelistVec))) ## Empty data frame. Each column is an empty vector representing an edge
  colnames(einv$t) <- edgelistVec ## Name the columns by edge in my empty data frame
  einv$boot <- data.frame(matrix(NA, nedges, 7)) ## Create an empty data frame for confidence intervals of edges
  einv$boot[,1] <- edgelistVec
  colnames(einv$boot) <- c("Edge","Network1", "Network2", "RealInv", "Estimate", "2.5% CI", "97.5% CI")
  
  diffcen$t <- data.frame(matrix(NA,it,nnodes)) ## Empty data frame. Each column is an empty vector representing a node
  colnames(diffcen$t) <- nodes
  diffcen$boot <- data.frame(matrix(NA, nnodes,7)) ## Empty data frame for confidence intervals of centrality differences
  diffcen$boot[,1] <- nodes
  colnames(diffcen$boot) <- c("Node","Network1", "Network2", "RealInv", "Estimate", "2.5% CI", "97.5% CI")
  
  
  ##### Determine which function to use ####
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
      return(IsingFit::IsingFit(x, AND=AND, progressbar=FALSE, plot=F)$weiadj)
    }
  } else if(match.arg(estimation_method)=="custom"){
    fun <- custom_estimation_method
  }
  
  ##### Calculate real values #####
  nw1.real <- fun(x1) 
  nw2.real <- fun(x2)
  glstrinv.real <- gl.str(nw1.real) - gl.str(nw2.real)
  glstrinv.sep <- c(gl.str(nw1.real), gl.str(nw2.real))
  diffedges.real <- abs(nw1.real - nw2.real)[upper.tri(abs(nw1.real - nw2.real))]
  diffedges.realmat <- matrix(diffedges.real, it, nedges, 
                              byrow = TRUE)
  diffedges.realoutput <- abs(nw1.real - nw2.real)
  einv$boot[,2] <- (nw1.real)[upper.tri((nw1.real), diag=FALSE)]
  einv$boot[,3] <- (nw2.real)[upper.tri((nw2.real), diag=FALSE)]
  einv$boot[,4] <- (nw1.real-nw2.real)[upper.tri((nw1.real-nw2.real), diag=FALSE)]
  einv$einv.real <- (nw1.real-nw2.real)
  
  if(test_centrality_node==TRUE){
    diffcen$boot[,2] <- centrality_auto(nw1.real)$node.centrality[nodes,centrality_node]
    diffcen$boot[,3] <- centrality_auto(nw2.real)$node.centrality[nodes,centrality_node]
    diffcen.real <- as.matrix(diffcen$boot[,2]) - as.matrix(diffcen$boot[,3])
    rownames(diffcen.real) <- nodes
    diffcen$boot[,4] <- diffcen.real
  }
  
  
  ##### Run the bootstrap #####
  for(i in 1:it) {
    if(match.arg(bootcut)=="cutEqual") {
      min_obs1 <- min_obs2 <- min(nrow(x1), nrow(x2))
    } else {
      min_obs1 <- nrow(x1)
      min_obs2 <- nrow(x2)
    }
    if(paired){
      samp1 <- samp2 <- sample(1:nrow(x1), size=min(c(min_obs1, min_obs2)), replace = TRUE)
    } else {
      samp1 <- sample(1:nrow(x1), size=min_obs1, replace = TRUE)
      samp2 <- sample(1:nrow(x2), size=min_obs2, replace = TRUE)
    }
    div1 <- x1[samp1,]
    div2 <- x2[samp2,]
    nw1 <- fun(div1) ## Generate a network for div1
    nw2 <- fun(div2) ## Generate a network for div2
    
    if(weighted==FALSE){      ## If unweighted, reassign values as 0 or 1
      nw1=(nw1!=0)*1
      nw2=(nw2!=0)*1 
    }
    
    glstrinv$t[i] <- gl.str(nw1)-gl.str(nw2) ## Calculate the global strength invariance for this bootstrapped sample
    einv$t[i,] <- (nw1-nw2)[upper.tri((nw1-nw2), diag=FALSE)]  ## Compute a vector of all edge invariances for this sample
    
    ## Compute centrality invariance (only works for 1 type of centrality at a time)
    if(test_centrality_node==TRUE){
      cen1 <- centrality_auto(nw1)$node.centrality[nodes,centrality_node]
      cen2 <- centrality_auto(nw2)$node.centrality[nodes,centrality_node]
      diffcen$t[i,] <- as.matrix(cen1) - as.matrix(cen2)
    }

    if (progressbar==TRUE) setTxtProgressBar(pb, i)
  }
  
  
  ##### Calculate estimated true values #####
  
  glstrinv.est <- quantile(glstrinv$t, probs = 0.5)
  for(i in 1:nedges) {
    einv$boot[i, 5] <- quantile(einv$t[,i], probs = 0.5)
  }
  if(test_centrality_node==TRUE){
    for(i in 1:nnodes){
      diffcen$boot[i,5] <- quantile(diffcen$t[,i], probs = 0.5)
    }
  }
  
  ##### Calculate confidence intervals #####
  
  glstrinv$ci <- quantile(glstrinv$t, probs = c(0.025, 0.975))
  for(i in 1:nedges) {
    einv$boot[i, 6:7] <- quantile(einv$t[,i], probs = c(0.025, 0.975))[1:2]
  }
  if(test_centrality_node==TRUE){
    for(i in 1:nnodes){
      diffcen$boot[i,6:7] <- quantile(diffcen$t[,i], probs = c(0.025, 0.975))[1:2]
    }
  }
  
  ##### Calculate network invariance #####
  
  nwinv.sig <- einv$boot$`2.5% CI` * einv$boot$`97.5% CI` #Trick to see if CI "crosses 0"
  sigfun <- function(x){if(x<0){FALSE}else{TRUE}}
  nwinv.sig2 <- sapply(nwinv.sig, FUN=sigfun)
  if(TRUE %in% nwinv.sig2){
    nwinv.num <- which(einv$boot$RealInv==max(abs(einv$boot$RealInv[nwinv.sig2])) | einv$boot$RealInv== -max(abs(einv$boot$RealInv[nwinv.sig2])))[1]
    nwinv.real <- einv$boot$RealInv[nwinv.num]
    nwinv.sep <- c(einv$boot$Network1[nwinv.num], einv$boot$Network2[nwinv.num])
    nwinv.est <- einv$boot$Estimate[nwinv.num]
    nwinv.ci <- c(einv$boot$`2.5% CI`[nwinv.num], einv$boot$`97.5% CI`[nwinv.num])
  }
  else{
    nwinv.num <- nwinv.real <- nwinv.sep <- nwinv.ci<- NA
    nwinv.est <- "Network invariance not significant"
  }

  ## Return instructions
  res <- list()
  res$glstrinv.real <- glstrinv.real
  res$glstrinv.sep <- glstrinv.sep
  res$glstrinv.est <- glstrinv.est
  res$glstrinv.ci <- glstrinv$ci
  if(include_t){
    res$glstrinv.t <- glstrinv$t
  }
  res$nwinv.real <- nwinv.real
  res$nwinv.sep <- nwinv.sep
  res$nwinv.est <- nwinv.est
  res$nwinv.ci <- nwinv.ci
  res$edges.tested <- edges
  if(test_edges){
    if(is.list(edges)){
      res$einv.mat <- einv$boot[einv$boot$Edge %in% reduced_edgelistVec,]
      if(include_t){
        res$einv.t <- einv$t[,colnames(einv$t) %in% reduced_edgelistVec]
      }
    } else {
      res$einv.mat <- einv$boot
      if(include_t){
        res$einv.t <- einv$t
      }
    }
    res$einv.real <- einv$einv.real
  }
  if(test_centrality_node){
    res$diffcen.real <- diffcen.real
    res$diffcen.mat <- diffcen$boot
    if(include_t){
      res$diffcen.t <- diffcen$t
    }
  }
  res$method <- "bootstrap"
  
  class(res) <- "NCT"
  return(res)
}

