NCT <- function(data1, 
                data2, 
                gamma, 
                it = 100, 
                binary.data=FALSE, 
                paired=FALSE, 
                weighted=TRUE, 
                AND=TRUE, 
                test.edges=FALSE, 
                edges, 
                progressbar=TRUE,
                make.positive.definite=TRUE,
                p.adjust.methods= c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"), # Edit from Payton
                estimator, # Assign this a function with an estimator (e.g., obtained from bootnet). This will OVERWRITE binary.data when used!
                estimatorArgs = list(), # Arguments sent to the estimator function
                
                # Optionally set a different estimator for group 2 (if missing overwritten with estimator and estimatorArgs)
                # estimator2, # Assign this a function with an estimator (e.g., obtained from bootnet). This will OVERWRITE binary.data when used!
                # estimatorArgs2 = list(), # Arguments sent to the estimator function
                verbose = TRUE
){ 
  p.adjust.methods <- match.arg(p.adjust.methods)
  
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
  dataall <- rbind(x1,x2)
  data.list <- list(x1,x2)
  b <- 1:(nobs1+nobs2)
  nvars <- ncol(x1)
  nedges <- nvars*(nvars-1)/2
  
  glstrinv.perm <- glstrinv.real <- nwinv.real <- nwinv.perm <- c()
  diffedges.perm <- matrix(0,it,nedges) 
  einv.perm.all <- array(NA,dim=c(nvars, nvars, it))
  corrpvals.all <- matrix(NA,nvars,nvars)
  edges.pvalmattemp <- matrix(0,nvars,nvars)
  
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
  
  ## Global strength invariance
  glstrinv.real <- abs(sum(abs(nw1[upper.tri(nw1)]))-sum(abs(nw2[upper.tri(nw2)])))
  # global strength of individual networks
  glstrinv.sep <- c(sum(abs(nw1[upper.tri(nw1)])), sum(abs(nw2[upper.tri(nw2)])))
  
  ## Individual edge invariance
  diffedges.real <- abs(nw1-nw2)[upper.tri(abs(nw1-nw2))] 
  diffedges.realmat <- matrix(diffedges.real,it,nedges,byrow=TRUE)
  diffedges.realoutput <- abs(nw1-nw2)
  
  ## Network structure invariance
  nwinv.real <- max(diffedges.real)
  
  
  #####################################
  #####     Start permutations    #####
  #####################################
  for (i in 1:it)
  {
    diffedges.permtemp <- matrix(0, nvars, nvars)
    
    # If not paired data
    if(paired==FALSE)
    {
      s <- sample(1:(nobs1+nobs2),nobs1,replace=FALSE)
      x1perm <- dataall[s,]
      x2perm <- dataall[b[-s],]
      
      
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
      s <- sample(c(1,2),nobs1,replace=TRUE)
      x1perm <- x1[s==1,]
      x1perm <- rbind(x1perm,x2[s==2,])
      x2perm <- x2[s==1,]
      x2perm <- rbind(x2perm,x1[s==2,])
      
      
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
    glstrinv.perm[i] <- abs(sum(abs(r1perm[upper.tri(r1perm)]))-sum(abs(r2perm[upper.tri(r2perm)])))
    diffedges.perm[i,] <- abs(r1perm-r2perm)[upper.tri(abs(r1perm-r2perm))]
    diffedges.permtemp[upper.tri(diffedges.permtemp, diag=FALSE)] <- diffedges.perm[i,]
    diffedges.permtemp <- diffedges.permtemp + t(diffedges.permtemp)
    einv.perm.all[,,i] <- diffedges.permtemp
    nwinv.perm[i] <- max(diffedges.perm[i,])
    
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
    edges.pvaltemp <- colSums(diffedges.perm >= diffedges.realmat)/it 
    
    ## If all edges should be tested
    if(is.character(edges))
    {
      # corrected p-values (or not if p.adjust.methods='none')
      corrpvals.all.temp <- round(p.adjust(edges.pvaltemp, method=p.adjust.methods),3)
      # matrix with corrected p values
      corrpvals.all
      corrpvals.all[upper.tri(corrpvals.all,diag=FALSE)] <- corrpvals.all.temp 
      rownames(corrpvals.all) <- colnames(corrpvals.all) <- colnames(x1)
      einv.pvals <- melt(corrpvals.all, na.rm=TRUE, value.name = 'p-value')
      einv.perm <- einv.perm.all
      einv.real <- diffedges.realoutput
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
      colnames(einv.pvals) <- c('Var1', 'Var2', 'p-value')
    }
    
    edges.tested <- colnames(einv.perm)
    
    res <- list(glstrinv.real = glstrinv.real,
                glstrinv.sep = glstrinv.sep,
                glstrinv.pval = sum(glstrinv.perm >= glstrinv.real)/it, 
                glstrinv.perm = glstrinv.perm,
                nwinv.real = nwinv.real,
                nwinv.pval = sum(nwinv.perm >= nwinv.real)/it, 
                nwinv.perm = nwinv.perm,
                edges.tested = edges.tested,
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
      glstrinv.pval = sum(glstrinv.perm >= glstrinv.real)/it, 
      glstrinv.perm = glstrinv.perm,
      nwinv.real = nwinv.real,
      nwinv.pval = sum(nwinv.perm >= nwinv.real)/it,
      nwinv.perm = nwinv.perm, 
      nw1 = nw1,
      nw2 = nw2
    )
  }
  
  # }
  
  # This is no obsolete (handled in section above)
  # 
  #   #################################
  #   ### procedure for binary data ###
  #   #################################
  #   ## Real data
  #   if(binary.data==TRUE) 
  #   {
  #     IF1 <- IsingFit(x1, AND = AND, gamma=gamma, plot=FALSE, progressbar=FALSE)
  #     IF2 <- IsingFit(x2, AND = AND, gamma=gamma, plot=FALSE, progressbar=FALSE)
  #     nw1 <- IF1$weiadj
  #     nw2 <- IF2$weiadj
  #     if(weighted==FALSE){
  #       nw1=(nw1!=0)*1
  #       nw2=(nw2!=0)*1
  #     }
  #     
  #     ##### Invariance measures #####
  #     
  #     ## Global strength invariance
  #     glstrinv.real <- abs(sum(abs(nw1[upper.tri(nw1)]))-sum(abs(nw2[upper.tri(nw2)])))
  #     # global strength of individual networks
  #     glstrinv.sep <- c(sum(abs(nw1[upper.tri(nw1)])), sum(abs(nw2[upper.tri(nw2)])))
  #     
  #     ## Individual edge invariance
  #     diffedges.real <- abs(nw1-nw2)[upper.tri(abs(nw1-nw2))] 
  #     diffedges.realmat <- matrix(diffedges.real,it,nedges,byrow=TRUE)
  #     diffedges.realoutput <- abs(nw1-nw2)
  #     
  #     ## Network structure invariance
  #     nwinv.real <- max(diffedges.real)
  #     
  #     #####################################
  #     #####     Start permutations    #####
  #     #####################################
  #     for (i in 1:it)
  #     {
  #       diffedges.permtemp <- matrix(0, nvars, nvars)
  #       if(paired==FALSE)
  #       {
  #         # a check to prevent an error of lognet(): a variable with only one 0 or 1 is not allowed
  #         checkN=0
  #         while(checkN==0){
  #           s <- sample(1:(nobs1+nobs2),nobs1,replace=FALSE)
  #           x1perm <- dataall[s,]
  #           x2perm <- dataall[b[-s],]
  #           
  #           cm1 <- colMeans(x1perm)
  #           cm2 <- colMeans(x2perm)
  #           checkN=ifelse(any(cm1<(1/nobs1))|any(cm1>(nobs1-1)/nobs1)|any(cm2<(1/nobs2))|any(cm2>(nobs2-1)/nobs2),0,1) 
  #         }
  #         IF1perm <- IsingFit(x1perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
  #         IF2perm <- IsingFit(x2perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
  #         r1perm <- IF1perm$weiadj
  #         r2perm <- IF2perm$weiadj
  #         if(weighted==FALSE){
  #           r1perm=(r1perm!=0)*1
  #           r2perm=(r2perm!=0)*1
  #         }
  #       }
  #       
  #       if(paired==TRUE)
  #       {
  #         checkN=0
  #         while(checkN==0)
  #         {
  #           s <- sample(c(1,2),nobs1,replace=TRUE)
  #           x1perm <- x1[s==1,]
  #           x1perm <- rbind(x1perm,x2[s==2,])
  #           x2perm <- x2[s==1,]
  #           x2perm <- rbind(x2perm,x1[s==2,])
  #           
  #           cm1 <- colMeans(x1perm)
  #           cm2 <- colMeans(x2perm)
  #           checkN=ifelse(any(cm1<(1/nobs1))|any(cm1>(nobs1-1)/nobs1)|any(cm2<(1/nobs2))|any(cm2>(nobs2-1)/nobs2),0,1) 
  #         }
  #         
  #         IF1perm <- IsingFit(x1perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
  #         IF2perm <- IsingFit(x2perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
  #         r1perm <- IF1perm$weiadj
  #         r2perm <- IF2perm$weiadj
  #         if(weighted==FALSE)
  #         {
  #           r1perm=(r1perm!=0)*1
  #           r2perm=(r2perm!=0)*1
  #         }
  #       }
  #       
  #       ## Invariance measures for permuted data
  #       glstrinv.perm[i] <- abs(sum(abs(r1perm[upper.tri(r1perm)]))-sum(abs(r2perm[upper.tri(r2perm)])))
  #       diffedges.perm[i,] <- abs(r1perm-r2perm)[upper.tri(abs(r1perm-r2perm))]
  #       diffedges.permtemp[upper.tri(diffedges.permtemp, diag=FALSE)] <- diffedges.perm[i,]
  #       diffedges.permtemp <- diffedges.permtemp + t(diffedges.permtemp)
  #       einv.perm.all[,,i] <- diffedges.permtemp
  #       nwinv.perm[i] <- max(diffedges.perm[i,])
  #       
  #       if (progressbar==TRUE) setTxtProgressBar(pb, i)
  #     }
  #     #####################################
  #     #####      End permutations     #####
  #     #####################################
  #     
  #     
  #     #####################################
  #     #####     Calculate p-values    #####
  #     #####################################
  #     if(test.edges==TRUE)
  #     {
  #       # vector with uncorrected p values
  #       edges.pvaltemp <- colSums(diffedges.perm >= diffedges.realmat)/it 
  #       # # matrix with uncorrected p values
  #       edges.pvalmattemp[upper.tri(edges.pvalmattemp,diag=FALSE)] <- edges.pvaltemp
  #       edges.pvalmattemp <- edges.pvalmattemp + t(edges.pvalmattemp)
  #       
  #       ## If all edges should be tested
  #       if(is.character(edges))
  #       {
  #         # corrected p-values (or not if p.adjust.methods='none')
  #         corrpvals.all.temp <- round(p.adjust(edges.pvaltemp, method=p.adjust.methods),3)
  #         # matrix with corrected p values
  #         corrpvals.all[upper.tri(corrpvals.all,diag=FALSE)] <- corrpvals.all.temp 
  #         rownames(corrpvals.all) <- colnames(corrpvals.all) <- colnames(x1)
  #         einv.pvals <- melt(corrpvals.all, na.rm=TRUE, value.name = 'p-value')
  #         einv.perm <- einv.perm.all
  #         einv.real <- diffedges.realoutput
  #       }
  #       
  #       ## If a selection of edges should be tested
  #       if(is.list(edges))
  #       {
  #         einv.perm <- matrix(NA,it,length(edges))
  #         colnames(einv.perm) <- edges
  #         uncorrpvals <- einv.real <- pairs <- c()
  #         for(j in 1:length(edges))
  #         {
  #           pairs <- rbind(pairs, c(colnames(x1)[edges[[j]][1]], colnames(x1)[edges[[j]][2]]))
  #           uncorrpvals[j] <- edges.pvalmattemp[edges[[j]][1],edges[[j]][2]]
  #           einv.real[j] <- diffedges.realoutput[edges[[j]][1],edges[[j]][2]]
  #           for(l in 1:it){
  #             einv.perm[l,j] <- einv.perm.all[,,l][edges[[j]][1],edges[[j]][2]]
  #           }
  #         }
  #         corrpvals <- p.adjust(uncorrpvals, method=p.adjust.methods)
  #         corrpvals_mat <- matrix(NA,length(edges),3)
  #         corrpvals_mat[,3] <- corrpvals
  #         
  #         corrpvals_mat[,1:2] <- pairs
  #         einv.pvals <- as.data.frame(corrpvals_mat)
  #         colnames(einv.pvals) <- c('Var1', 'Var2', 'p-value')
  #       }
  #       
  #       edges.tested <- colnames(einv.perm)
  #       
  #       res <- list(nw1 = nw1,
  #                   nw2 = nw2,
  #                   glstrinv.real = glstrinv.real,
  #                   glstrinv.sep = glstrinv.sep,
  #                   glstrinv.pval = sum(glstrinv.perm >= glstrinv.real)/it, 
  #                   glstrinv.perm = glstrinv.perm,
  #                   nwinv.real = nwinv.real,
  #                   nwinv.pval = sum(nwinv.perm >= nwinv.real)/it, 
  #                   nwinv.perm = nwinv.perm,
  #                   edges.tested = edges.tested,
  #                   einv.real = einv.real,
  #                   einv.pvals = einv.pvals,
  #                   einv.perm = einv.perm)
  #     }
  #     
  #     
  #     if (progressbar==TRUE) close(pb)
  #     
  #     if(test.edges==FALSE) 
  #     {
  #       res <- list(nw1 = nw1,
  #                   nw2 = nw2,
  #                   glstrinv.real = glstrinv.real, 
  #                   glstrinv.sep = glstrinv.sep,
  #                   glstrinv.pval = sum(glstrinv.perm >= glstrinv.real)/it,
  #                   glstrinv.perm = glstrinv.perm,
  #                   nwinv.real = nwinv.real,
  #                   nwinv.pval = sum(nwinv.perm >= nwinv.real)/it,
  #                   nwinv.perm = nwinv.perm
  #       )
  #     }
  #   }
  
  class(res) <- "NCT"
  return(res)
}


## Methods:

summary.NCT <- function(x,...){
  
  cat("\n NETWORK INVARIANCE TEST
      Test statistic M: ", x$nwinv.real,
      "\n p-value", x$nwinv.pval,
      "\n\n GLOBAL STRENGTH INVARIANCE TEST
      Global strength per group: ", x$glstrinv.sep,
      "\n Test statistic S: ", x$glstrinv.real,
      "\n p-value", x$glstrinv.pval,
      "\n\n EDGE INVARIANCE TEST
      Edges tested: ", x$edges.tested,
      "\n Test statistic E: ", x$einv.real,
      "\n p-value", x$einv.pvals
  )
  
}

plot.NCT <- function(x,what = c("strength","network","edge"),...){
  
  what <- match.arg(what)
  
  ## Plot results of global strength invariance test (not reliable with only 10 permutations!):
  if (what == "strength"){
    hist(x$glstrinv.perm, main=paste('p =',x$glstrinv.pval),xlab='Difference in global strength',xlim=c(0,max(x$glstrinv.real,x$glstrinv.perm)))
    points(x$glstrinv.real,0,col='red',pch=17)
  } 
  
  if (what == "network"){
    
    ## Plot results of the network invariance test (not reliable with only 10 permutations!):
    hist(x$nwinv.perm, main=paste('p =',x$nwinv.pval),xlab='Maximum of difference',xlim=c(0,max(x$nwinv.real,x$nwinv.perm)))
    points(x$nwinv.real,0,col='red',pch=17)
  } 
  
  if (what == "edge"){
    
    ## Plot results of the maximum difference in edge weights (not reliable with only 10 permutations)
    nedgetests <- ncol(x$einv.perm)
    for(i in 1:nedgetests){
      hist(x$einv.perm[,i], main=paste('p =',x$einv.pval[i,3]),xlab=paste('Difference in edge strength ', colnames(x$einv.perm)[i]),xlim=c(0,max(x$einv.real,x$einv.perm)))
      points(x$einv.real[i],0,col='red',pch=17)
    }
    
  } #else stop("Method not implemented yet.")
  
}