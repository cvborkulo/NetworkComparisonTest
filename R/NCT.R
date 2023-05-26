#' NetworkComparisonTest: Statistical Comparison of Two Networks Based on Several Invariance Measures
#'
#' @description This permutation based hypothesis test, suited for several types of data supported by the estimateNetwork function of the bootnet package (Epskamp & Fried, 2018), assesses the difference between two networks based on several invariance measures (network structure invariance, global strength invariance, edge invariance, several centrality measures, etc.). Network structures are estimated with l1-regularization. The Network Comparison Test is suited for comparison of independent (e.g., two different groups) and dependent samples (e.g., one group that is measured twice).
#'
#' @name NCT
#' 
#' @param data1 One of two datasets. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables. Must be cross-sectional data. Can also be the result of \code{estimateNetwork} from the bootnet package.
#' @param data2 The other of two datasets. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables. Must be cross-sectional data.  Can also be the result of \code{estimateNetwork} from the bootnet package.
#' @param gamma A single value between 0 and 1. When not entered, gamma is set to 0.25 for binary data and 0.50 for gaussian data. Networks are estimated with this value for hyperparameter gamma in the extended BIC.
#' @param it The number of iterations (permutations).
#' @param binary.data Logical. Can be TRUE or FALSE to indicate whether the data is binary or not. If binary.data is FALSE, the data is regarded gaussian. This argument is ignored when using estimateNetwork() output as input for NCT.
#' @param paired Logical. Can be TRUE of FALSE to indicate whether the samples are dependent or not. If paired is TRUE, relabeling is performed within each pair of observations. If paired is FALSE, relabeling is not restricted to pairs of observations. Note that, currently, dependent data is assumed to entail one group measured twice.
#' @param weighted Logical. Can be TRUE of FALSE to indicate whether the networks to be compared should be weighted of not. If not, the estimated networks are dichotomized. Defaults to TRUE.
#' @param AND Logical. Can be TRUE of FALSE to indicate whether the AND-rule or the OR-rule should be used to define the edges in the network. Defaults to TRUE. Only necessary for binary data.
#' @param abs Logical. Should global strength consider the absolute value of edge weights, or the raw value (i.e., global expected influence)?
#' @param test.edges Logical. Can be TRUE of FALSE to indicate whether or not differences in individual edges should be tested.
#' @param edges Character or list. When 'all', differences between all individual edges are tested. When provided a list with one or more pairs of indices referring to variables, the provided edges are tested.
#' @param progressbar Logical. Should the pbar be plotted in order to see the progress of the estimation procedure? Defaults to TRUE.
#' @param make.positive.definite If \code{make.positive.definite = TRUE}, the covariance matrices used for the glasso are projected to the nearest positive definite matrices, if they are not yet positive definite. This is useful for small n, for which it is very likely that at least one of the bootstrap comparisons involves a covariance matrix that is not positive definite.
#' @param p.adjust.methods Character. Can be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", or "none". To control (or not) for testing of multiple edges. Defaults to "none".
#' @param test.centrality Logical. Should centrality metrics be compared across networks?
#' @param centrality Type of centrality metrics to test. Can be any of c("all", "closeness", "betweenness", "strength", "expectedInfluence", "bridgeStrength", "bridgeCloseness", "bridgeBetweenness", "bridgeExpectedInfluence")
#' @param nodes Specific nodes for centrality tests. Can be character names or index numbers. Only used if test.centrality=TRUE
#' @param communities Passed to bridge() if computing bridge centrality
#' @param useCommunities Passed to bridge() if computing bridge centrality
#' @param estimator A function that takes data as input and returns a network structure. This can be used for custom estimation algorithms. Note, supplying this function will overwrite the arguments \code{binary.data}, \code{AND}, \code{gamma} and \code{make.positive.definite}.
#' @param estimatorArgs Arguments to the \code{estimator} function.
#' @param verbose Logical: Should some warnings and notes be printed?
#'
#' @return NCT returns a 'NCT' object that contains the following items:
#' 
#' \itemize{
#' \item \code{glstrinv.real } The difference in global strength between the networks of the observed data sets.
#' \item \code{glstrinv.perm } The difference in global strength between the networks of the permutated data sets.
#' \item \code{glstrinv.sep} The global strength values of the individual networks
#' \item \code{glstrinv.pval} The p value resulting from the permutation test concerning difference in global strength.
#' \item \code{nwinv.real} The value of the maximum difference in edge weights of the observed networks.
#' \item \code{nwinv.perm} The values of the maximum difference in edge weights of the permuted networks.
#' \item \code{nwinv.pval} The p value resulting from the permutation test concerning the maximum difference in edge weights.
#' \item \code{einv.pvals} p-values (corrected for multiple testing or not according to 'p.adjust.methods') per edge from the permutation test concerning differences in edges weights. Only returned if test.edges = TRUE.
#' \item \code{einv.real} The value of the difference in edge weight of the observed networks (multiple values if more edges are called to test). Only if test.edges = TRUE.
#' \item \code{einv.perm} The values of the difference in edge weight of the permuted networks. Only if test.edges = TRUE.
#' \item \code{diffcen.real} The values of the difference in centralities of the observed networks. Only if test.centrality = TRUE.
#' \item \code{diffcen.perm} The values of the difference in centralities of the permuted networks. Only if test.centrality = TRUE.
#' \item \code{diffcen.pval} p-values(corrected for multiple testing or not according to 'p.adjust.methods') per node from the permutation test concerning differences in centralities. Only if test.centrality = TRUE.
#' }
#' 
#' @importFrom graphics hist points
#' @importFrom stats cor p.adjust
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Matrix nearPD
#' @importFrom reshape2 melt
#' @importFrom stats na.omit
#' @importFrom methods is
#' 
#' 
#' @export
#'

NCT <- function(data1, data2, 
                gamma, it = 100, binary.data=FALSE, 
                paired=FALSE, weighted=TRUE, AND=TRUE, abs=TRUE,
                test.edges=FALSE, edges="all", 
                progressbar=TRUE, make.positive.definite=TRUE,
                p.adjust.methods= c("none", "holm", "hochberg", "hommel", 
                                    "bonferroni", "BH", "BY", "fdr"), 
                test.centrality=FALSE, 
                centrality=c("strength","expectedInfluence"),nodes="all",
                communities=NULL,useCommunities="all",
                estimator, estimatorArgs = list(), 
                verbose = TRUE){ 
  
  # store function call, including default arguments not explicitly set - credit to Neal Fultz
  match.call.defaults <- function(...) {
    # Extract explicit call arguments
    call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
    # Extract default call arguments
    formals <- evalq(formals(), parent.frame(1))
    
    # When extracting calls, the symbols T or F not properly mapped to TRUE and FALSE
    clean.call.arg <- function(arg){
      if(is.null(arg)){
        return(list(arg))
      }
      if(arg == "T"){
        return(list(TRUE))
      }
      if(arg == "F"){
        return(list(FALSE))
      }
      return(list(arg))
    }
    for(i in 1:length(names(call))){
      call[i] <- clean.call.arg(call[[i]])
    }
    
    # if default argument not explicitly written, add the default value to the saved call
    for(i in setdiff(names(formals), names(call))){
      call[i] <- list(formals[[i]])
    }
    
    match.call(sys.function(sys.parent()), call)
  }
  cl <- match.call.defaults()
  
  p.adjust.methods <- match.arg(p.adjust.methods)
  
  # Fix for networktools example:
  if (missing(edges)) edges <- "all"
  
  # Small test to warn people of arguments that are ignored if the input is a bootnet object:
  # Test for bootnet objects:
  if (is(data1,"bootnetResult") || is(data2,"bootnetResult")){
    mc <- match.call()
    if ("gamma" %in% names(mc)){
      if (verbose) message("Note: Input is a bootnetResult object, argument 'gamma' is ignored.")
    }
    
    if ("binary.data" %in% names(mc)){
      if (verbose) message("Note: Input is a bootnetResult object, argument 'binary.data' is ignored.")
    }
    
    if ("AND" %in% names(mc)){
      if (verbose) message("Note: Input is a bootnetResult object, argument 'AND' is ignored.")
    }
    
    if ("make.positive.definite" %in% names(mc)){
      if (verbose) message("Note: Input is a bootnetResult object, argument 'make.positive.definite' is ignored.")
    }
    
    # Check if estimator function is used:
    if (!missing(estimator)){
      stop("Custom estimator function not supported for bootnet objects.")
    }
  }
  
  # If object 1 is a bootnetResult, extract data, estimator and args:
  if (is(data1,"bootnetResult")){
    if (verbose) message("Note: estimateNetwork object used - estimation method has possibly not been validated.")
    # Estimator:
    # if (missing(estimator)){
    estimator <- data1$estimator
    # }
    # Arguments:
    # if (missing(estimatorArgs)){
    estimatorArgs <- data1$arguments
    estimatorArgs$verbose <- FALSE
    # }
    # Data:
    data1 <- data1$data
  }
  
  # If object 2 is a bootnetResult, extract data, estimator and args:
  if (is(data2,"bootnetResult")){
    # Estimator:
    # if (missing(estimator2)){
    estimator2 <- data2$estimator
    # }
    # Arguments:
    # if (missing(estimatorArgs2)){
    estimatorArgs2 <- data2$arguments
    estimatorArgs2$verbose <- FALSE
    # }
    
    # Test if estimation methods are identical:
    if (!identical(estimator,estimator2)){
      stop("Estimation methods are not identical.")
    }
    
    # Test if arguments are identical:
    # Test if estimation methods are identical:
    if (!identical(estimatorArgs,estimatorArgs2)){
      stop("Estimation arguments are not identical.")
    }
    
    # Data:
    data2 <- data2$data
  }
  
  
  # Gamma (note, not used for bootnet objects):
  if (missing(gamma)){
    if (binary.data){
      gamma <- 0.25
    } else {
      gamma <- 0.5
    }
  } 
  
  # Estimator function:
  if (missing(estimator)){
    
    if (binary.data){
      estimator <- NCT_estimator_Ising
      estimatorArgs$AND <- AND
    } else {
      estimator <- NCT_estimator_GGM
      estimatorArgs$make.positive.definite <- make.positive.definite
    }
    estimatorArgs$gamma <- gamma
  } else {
    # Look at the call if someone also used "binary.data":
    mc <- match.call()
    if ("binary.data" %in% names(mc)){
      if (verbose) message("Note: Both 'estimator' and 'binary.data' arguments used: only the 'estimator' will be used ('binary.data' will be ignored)")
    }
  }
  
  
  if (progressbar==TRUE) pb <- txtProgressBar(max=it, style = 3)
  x1 <- data1
  x2 <- data2
  nobs1 <- nrow(x1)
  nobs2 <- nrow(x2)
  if(is.null(colnames(x1)) && is.null(colnames(x2))){
    colnames(x1) <- colnames(x2) <- paste("var",1:ncol(x1),sep="")
  }
  dataall <- rbind(x1,x2)
  data.list <- list(x1,x2)
  b <- 1:(nobs1+nobs2)
  nvars <- ncol(x1)
  nedges <- nvars*(nvars-1)/2
  nnodes <- ifelse(nodes[1]=="all",nvars,length(nodes))
  nodes <- if(is.numeric(nodes)){colnames(data1)[nodes]} else{nodes}
  if(is.list(edges)){
    edges.tested <- edges
    if(is.character(edges[[1]])){
      whichfun <- function(x){which(colnames(data1)%in%x)}
      edges <- lapply(edges,whichfun)
    }
  }
  
  glstrinv.perm <- glstrinv.real <- nwinv.real <- nwinv.perm <- c()
  diffedges.perm <- matrix(0,it,nedges) 
  einv.perm.all <- array(NA,dim=c(nvars, nvars, it))
  corrpvals.all <- matrix(NA,nvars,nvars)
  edges.pvalmattemp <- matrix(0,nvars,nvars)
  
  
  validCentrality <- c("closeness", "betweenness", 
                       "strength", "expectedInfluence", "bridgeStrength", 
                       "bridgeCloseness", "bridgeBetweenness", "bridgeExpectedInfluence")
  bridgecen <- c("bridgeStrength", "bridgeBetweenness", 
                 "bridgeCloseness", "bridgeExpectedInfluence")
  centrality <- if(centrality[1]=="all") {
    validCentrality
  }else {
    centrality
  }
  diffcen.perm <- matrix(NA, it, nnodes*length(centrality))
  
  #####################################
  ###    procedure for all data     ###
  #####################################
  
  # Estimate the networks:
  nw1 <- do.call(estimator,c(list(x1),estimatorArgs))
  if (is.list(nw1)) nw1 <- nw1$graph
  
  nw2 <- do.call(estimator,c(list(x2),estimatorArgs))
  if (is.list(nw2)) nw2 <- nw2$graph
  
  if(weighted==FALSE){
    nw1=(nw1!=0)*1
    nw2=(nw2!=0)*1
  }
  ##### Invariance measures #####
  
  if(abs){
    ## Global strength invariance
    glstrinv.real <- abs(sum(abs(nw1[upper.tri(nw1)]))-sum(abs(nw2[upper.tri(nw2)])))
    # global strength of individual networks
    glstrinv.sep <- c(sum(abs(nw1[upper.tri(nw1)])), sum(abs(nw2[upper.tri(nw2)])))
  } else {
    ## Global strength invariance
    glstrinv.real <- abs(sum(nw1[upper.tri(nw1)])-sum(nw2[upper.tri(nw2)]))
    # global strength of individual networks
    glstrinv.sep <- c(sum(nw1[upper.tri(nw1)]), sum(nw2[upper.tri(nw2)]))
  }
  
  
  ## Individual edge invariance
  diffedges.real <- abs(nw1-nw2)[upper.tri(abs(nw1-nw2))] 
  diffedges.realmat <- matrix(diffedges.real,it,nedges,byrow=TRUE)
  diffedges.realoutput <- abs(nw1-nw2)
  
  ## Network structure invariance
  nwinv.real <- max(diffedges.real)
  
  ## Centrality invariance
  
  if(test.centrality==TRUE){
    if (!all(centrality %in% validCentrality)) {
      stop(paste0("'centrality' must be one of: ", paste0("'", 
                                                          validCentrality, "'", collapse = ", ")))
    }
    cen1 <- centrality_auto(nw1)$node.centrality
    cen2 <- centrality_auto(nw2)$node.centrality
    names(cen1) <- names(cen2) <- c("betweenness","closeness","strength","expectedInfluence")
    if(TRUE %in% (bridgecen %in% centrality)){
      b1 <- networktools::bridge(nw1, communities=communities, useCommunities=useCommunities)
      b2 <- networktools::bridge(nw2, communities=communities, useCommunities=useCommunities)
      names(b1) <- names(b2) <- c(bridgecen, "bridgeExpectedInfluence2step", 
                                  "communities")
      b1$communities <- b2$communities <- NULL
      cen1 <- data.frame(c(cen1,b1))
      cen2 <- data.frame(c(cen2,b2))
    }
    diffcen.real <- as.matrix(cen1) - as.matrix(cen2)
  }
  
  
  #####################################
  #####     Start permutations    #####
  #####################################
  # warning when paired data is compared
  if(paired==TRUE)
  {
    if (verbose) message("Note: NCT for dependent data has not been validated.")
  }
  
  for (i in 1:it)
  {
    diffedges.permtemp <- matrix(0, nvars, nvars)
    
    # If not paired data
    if(paired==FALSE)
    {
      # Include variance check
      okay <- FALSE
      counter <- 0
      
      if(binary.data) { # if binary data we need to resample the permutation to ensure mininum required variance for glmnet
        while(okay == FALSE) {
          
          # Permute
          s <- sample(1:(nobs1+nobs2),nobs1,replace=FALSE)
          x1perm <- dataall[s,]
          x2perm <- dataall[b[-s], ]
          
          # check glmnet requirement: at least two instances of each category
          ind <- all(apply(x1perm, 2,  function(x) min(c(sum(x==0), sum(x==1)))) > 1) & all(apply(x2perm, 2,  function(x) min(c(sum(x==0), sum(x==1)))) > 1)
          if(ind) okay <- TRUE else counter <- counter + 1
          
        } # end: while
      } else{
        s <- sample(1:(nobs1+nobs2),nobs1,replace=FALSE)
        x1perm <- dataall[s,]
        x2perm <- dataall[b[-s], ]
      }
      
      
      # Estimate the networks:
      r1perm <- do.call(estimator,c(list(x1perm),estimatorArgs))
      if (is.list(r1perm)) r1perm <- r1perm$graph
      
      r2perm <- do.call(estimator,c(list(x2perm),estimatorArgs))
      if (is.list(r2perm)) r2perm <- r2perm$graph
      
      if(weighted==FALSE){
        r1perm=(r1perm!=0)*1
        r2perm=(r2perm!=0)*1
      }
    }
    
    # If paired data
    if(paired==TRUE)
      
    {
      # Include variance check
      okay <- FALSE
      counter <- 0
      
      if(binary.data) { # if binary data we need to resample the permutation to ensure mininum required variance for glmnet
        while(okay == FALSE) {
                    s <- sample(c(1,2),nobs1,replace=TRUE)
          x1perm <- x1[s==1,]
          x1perm <- rbind(x1perm,x2[s==2,])
          x2perm <- x2[s==1,]
          x2perm <- rbind(x2perm,x1[s==2,])
          
          # check glmnet requirement: at least two instances of each category
          ind <- all(apply(x1perm, 2,  function(x) min(c(sum(x==0), sum(x==1)))) > 1) & all(apply(x2perm, 2,  function(x) min(c(sum(x==0), sum(x==1)))) > 1)
          if(ind) okay <- TRUE else counter <- counter + 1
          
        } # end: while
      } else{
        s <- sample(c(1,2),nobs1,replace=TRUE)
        x1perm <- x1[s==1,]
        x1perm <- rbind(x1perm,x2[s==2,])
        x2perm <- x2[s==1,]
        x2perm <- rbind(x2perm,x1[s==2,])
      }
      
      # Estimate the networks:
      r1perm <- do.call(estimator,c(list(x1perm),estimatorArgs))
      if (is.list(r1perm)) r1perm <- r1perm$graph
      
      r2perm <- do.call(estimator,c(list(x2perm),estimatorArgs))
      if (is.list(r2perm)) r2perm <- r2perm$graph
      
      if(weighted==FALSE){
        r1perm=(r1perm!=0)*1
        r2perm=(r2perm!=0)*1
      }
    }
    
    ## Invariance measures for permuted data
    if(abs){
      glstrinv.perm[i] <- abs(sum(abs(r1perm[upper.tri(r1perm)]))-sum(abs(r2perm[upper.tri(r2perm)])))
    } else {
      glstrinv.perm[i] <- abs(sum(r1perm[upper.tri(r1perm)])-sum(r2perm[upper.tri(r2perm)]))
    }
    
    diffedges.perm[i,] <- abs(r1perm-r2perm)[upper.tri(abs(r1perm-r2perm))]
    diffedges.permtemp[upper.tri(diffedges.permtemp, diag=FALSE)] <- diffedges.perm[i,]
    diffedges.permtemp <- diffedges.permtemp + t(diffedges.permtemp)
    einv.perm.all[,,i] <- diffedges.permtemp
    nwinv.perm[i] <- max(diffedges.perm[i,])
    
    
    if(test.centrality==TRUE){
      cen1permtemp <- centrality_auto(r1perm)$node.centrality
      cen2permtemp <- centrality_auto(r2perm)$node.centrality
      names(cen1permtemp) <- names(cen2permtemp) <- c("betweenness","closeness","strength","expectedInfluence")
      if(TRUE %in% (bridgecen %in% centrality)){
        b1permtemp <- networktools::bridge(r1perm, communities=communities, useCommunities=useCommunities)
        b2permtemp <- networktools::bridge(r2perm, communities=communities, useCommunities=useCommunities)
        names(b1permtemp) <- names(b2permtemp) <- c(bridgecen, "bridgeExpectedInfluence2step", 
                                                    "communities")
        b1permtemp$communities <- b2permtemp$communities <- NULL
        cen1permtemp <- data.frame(c(cen1permtemp,b1permtemp))
        cen2permtemp <- data.frame(c(cen2permtemp,b2permtemp))
      }
      diffcen.permtemp <- as.matrix(cen1permtemp) - as.matrix(cen2permtemp)
      if(nodes[1]=="all"){
        diffcen.perm[i,] <- reshape2::melt(diffcen.permtemp[,centrality])$value
      } else {
        diffcen.perm[i,] <- reshape2::melt(diffcen.permtemp[which(nodes%in%colnames(data1)),centrality])$value
      } 
    }
    
    
    if (progressbar==TRUE) setTxtProgressBar(pb, i)
  }
  #####################################
  #####      End permutations     #####
  #####################################
  
  
  #####################################
  #####     Calculate p-values    #####
  #####################################
  if(test.edges==TRUE)
  {
    # vector with uncorrected p values
    edges.pvaltemp <- (colSums(diffedges.perm >= diffedges.realmat) + 1) / (it + 1)
    
    ## If all edges should be tested
    if(is.character(edges))
    {
      # corrected p-values (or not if p.adjust.methods='none')
      corrpvals.all.temp <- p.adjust(edges.pvaltemp, method=p.adjust.methods)
      # matrix with corrected p values
      corrpvals.all
      corrpvals.all[upper.tri(corrpvals.all,diag=FALSE)] <- corrpvals.all.temp 
      rownames(corrpvals.all) <- colnames(corrpvals.all) <- colnames(x1)
      einv.pvals <- melt(corrpvals.all, na.rm=TRUE, value.name = 'p-value')
      einv.perm <- einv.perm.all
      einv.real <- diffedges.realoutput 
      einv.pvals <- cbind(einv.pvals, round(einv.real[upper.tri(einv.real)],8))
      colnames(einv.pvals) <- c('Var1', 'Var2', 'p-value', "Test statistic E")
    }
    
    ## If a selection of edges should be tested
    if(is.list(edges))
    {
      einv.perm <- matrix(NA,it,length(edges))
      colnames(einv.perm) <- edges
      uncorrpvals <- einv.real <- pairs <- c()
      
      # matrix with uncorrected p values
      edges.pvalmattemp[upper.tri(edges.pvalmattemp,diag=FALSE)] <- edges.pvaltemp
      edges.pvalmattemp <- edges.pvalmattemp + t(edges.pvalmattemp)
      
      for(j in 1:length(edges))
      {
        pairs <- rbind(pairs, c(colnames(x1)[edges[[j]][1]], colnames(x1)[edges[[j]][2]]))
        uncorrpvals[j] <- edges.pvalmattemp[edges[[j]][1],edges[[j]][2]]
        einv.real[j] <- diffedges.realoutput[edges[[j]][1],edges[[j]][2]]
        for(l in 1:it){
          einv.perm[l,j] <- einv.perm.all[,,l][edges[[j]][1],edges[[j]][2]]
        }
      }
      corrpvals <- p.adjust(uncorrpvals, method=p.adjust.methods)
      corrpvals_mat <- matrix(NA,length(edges),3)
      corrpvals_mat[,3] <- corrpvals
      corrpvals_mat[,1:2] <- pairs
      einv.pvals <- as.data.frame(corrpvals_mat)
      einv.pvals <- cbind(einv.pvals, einv.real)
      colnames(einv.pvals) <- c('Var1', 'Var2', 'p-value', "Test statistic E")
    }
    
    res <- list(glstrinv.real = glstrinv.real,
                glstrinv.sep = glstrinv.sep,
                glstrinv.pval = (sum(glstrinv.perm >= glstrinv.real) + 1) / (it + 1), 
                glstrinv.perm = glstrinv.perm,
                nwinv.real = nwinv.real,
                nwinv.pval = (sum(nwinv.perm >= nwinv.real) + 1) / (it + 1), 
                nwinv.perm = nwinv.perm,
                einv.real = einv.real,
                einv.pvals = einv.pvals,
                einv.perm = einv.perm, 
                nw1 = nw1,
                nw2 = nw2)
    
  }
  
  if (progressbar==TRUE) close(pb)
  
  if(test.edges==FALSE) 
  {
    res <- list(
      glstrinv.real = glstrinv.real, 
      glstrinv.sep = glstrinv.sep,
      glstrinv.pval = (sum(glstrinv.perm >= glstrinv.real) + 1) / (it + 1), 
      glstrinv.perm = glstrinv.perm,
      nwinv.real = nwinv.real,
      nwinv.pval = (sum(nwinv.perm >= nwinv.real) + 1) / (it + 1),
      nwinv.perm = nwinv.perm, 
      nw1 = nw1,
      nw2 = nw2
    )
  }
  
  if(test.centrality){
    if(nodes[1]=="all"){
      diffcen.real.vec <- reshape2::melt(diffcen.real[,centrality])$value
    } else {
      diffcen.real.vec <- reshape2::melt(diffcen.real[which(nodes%in%colnames(data1)),centrality])$value
    } 
    diffcen.realmat <- matrix(diffcen.real.vec, it, nnodes*length(centrality), 
                              byrow = TRUE)
    diffcen.pvaltemp <- (colSums(abs(diffcen.perm) >= abs(diffcen.realmat)) + 1) / (it + 1)
    diffcen.HBall <- p.adjust(diffcen.pvaltemp, method = p.adjust.methods)
    diffcen.pval <- matrix(diffcen.HBall, nnodes, length(centrality))
    diffcen.real <-  matrix(diffcen.real.vec, nrow=nnodes,ncol=length(centrality))
    colnames(diffcen.pval) <- colnames(diffcen.real) <- centrality
    res[["diffcen.real"]] <- diffcen.real
    res[["diffcen.perm"]] <- diffcen.perm
    res[["diffcen.pval"]] <- diffcen.pval
    
    #Column & row names
    if(nodes[1]=="all"){
      rownames(res[["diffcen.real"]]) <- rownames(res[["diffcen.pval"]]) <- colnames(data1)
      colnames(res[["diffcen.perm"]]) <- apply(expand.grid(colnames(data1), centrality), 1, paste, collapse=".")
    }
    else {
      rownames(res[["diffcen.real"]]) <- rownames(res[["diffcen.pval"]]) <- nodes
    }
  }
  
  res$info$call <- cl
  
  class(res) <- "NCT"
  return(res)
}


## Methods:

## TODO: add option for global expected influence

