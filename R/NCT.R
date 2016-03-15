NCT <- function(data1, data2, gamma, it, binary.data=FALSE, paired=FALSE, weighted=TRUE, AND=TRUE, test.edges=FALSE, edges, progressbar=TRUE, ...){ 
  
  # weighted: are the connection weights taken into account?
  # paired: default is FALSE, so applicable to independent samples. May be unequal sample sizes. When TRUE: samples are dependent. Sample sizes must be equal. Relabeling within each pair of (dependent) observations.
  # AND: default is TRUE. When using binary data, IsingFit will apply the AND-rule for determining whether there is a connection or not.
  # test.edges: default is FALSE and no differences in individual edges are tested. When TRUE the edges provided in argument 'edges' are tested.
  # edges: a list should be provided with the index (indices) of the edge(s) to be tested. A Holm-Bonferroni correction is applied to correct for multiple testing.
  # progressbar: default is TRUE. Shows the progression of the analysis.
  
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
  
  #####################################
  ### procedure for non-binary data ###
  #####################################
  # Real data
  if(binary.data==FALSE) 
  {
    glstr.perm <- glstr.real <- max.real <- max.perm <- c()
    diffedges.perm <- matrix(NA,it,nedges) 
    edges.pval.HBall <- matrix(NA,nvars,nvars)
    edges.pvalmattemp <- matrix(0,nvars,nvars)
    
    # real data
    nw1 <- EBICglasso(cor(x1),nrow(x1),gamma=gamma)
    nw2 <- EBICglasso(cor(x2),nrow(x2),gamma=gamma)
    if(weighted==FALSE){
      nw1=(nw1!=0)*1
      nw2=(nw2!=0)*1
    }
    
    # Difference measures
    # difference in global strength
    glstr.real <- abs(sum(abs(nw1[upper.tri(nw1)]))-sum(abs(nw2[upper.tri(nw2)])))
    # global strength of individual networks
    glstr.sep <- c(sum(abs(nw1[upper.tri(nw1)])), sum(abs(nw2[upper.tri(nw2)])))
    # difference in weights by edge
    diffedges.real <- abs(nw1-nw2)[upper.tri(abs(nw1-nw2))] 
    # maximum norm
    max.real <- max(diffedges.real)
    
    # permuted data
    for (i in 1:it)
    {
      if(paired==FALSE)
      {
        s <- sample(1:(nobs1+nobs2),nobs1,replace=FALSE)
        x1perm <- dataall[s,]
        x2perm <- dataall[b[-s],]
        r1perm <- EBICglasso(cor(x1perm),nrow(x1perm),gamma=gamma)
        r2perm <- EBICglasso(cor(x2perm),nrow(x2perm),gamma=gamma)
        if(weighted==FALSE){
          r1perm=(r1perm!=0)*1
          r2perm=(r2perm!=0)*1
        }
      }
      
      if(paired==TRUE)
      {
        s <- sample(c(1,2),nobs1,replace=TRUE)
        x1perm <- x1[s==1,]
        x1perm <- rbind(x1perm,x2[s==2,])
        x2perm <- x2[s==1,]
        x2perm <- rbind(x2perm,x1[s==2,])
        r1perm <- EBICglasso(cor(x1perm),nrow(x1perm),gamma=gamma)
        r2perm <- EBICglasso(cor(x2perm),nrow(x2perm),gamma=gamma)
        if(weighted==FALSE){
          r1perm=(r1perm!=0)*1
          r2perm=(r2perm!=0)*1
        }
      }
      
      # Difference measures for permuted data
      glstr.perm[i] <- abs(sum(abs(r1perm[upper.tri(r1perm)]))-sum(abs(r2perm[upper.tri(r2perm)])))
      diffedges.perm[i,] <- abs(r1perm-r2perm)[upper.tri(abs(r1perm-r2perm))]
      max.perm[i] <- max(diffedges.perm[i,])
      
      if (progressbar==TRUE) setTxtProgressBar(pb, i)
    }
    
    diffedges.realmat <- matrix(diffedges.real,it,nedges,byrow=TRUE)
    # vector with uncorrected p values
    edges.pvaltemp <- colSums(diffedges.perm >= diffedges.realmat)/it 
    # matrix with uncorrected p values
    edges.pvalmattemp[upper.tri(edges.pvalmattemp,diag=FALSE)] <- edges.pvaltemp
    edges.pvalmattemp <- edges.pvalmattemp + t(edges.pvalmattemp)
    
    if(test.edges==TRUE)
    {
      if(is.character(edges))
      {
        ept.HBall <- p.adjust(edges.pvaltemp, method='holm')
        # matrix with corrected p values
        edges.pval.HBall[upper.tri(edges.pval.HBall,diag=FALSE)] <- ept.HBall 
        rownames(edges.pval.HBall) <- colnames(edges.pval.HBall) <- colnames(data1)
        el.pvals <- melt(edges.pval.HBall, na.rm=TRUE, value.name = 'p-value')
      }
      
      if(is.list(edges))
      {
        uncorrpvals <- c()
        for(i in 1:length(edges))
        {
          uncorrpvals[i] <- edges.pvalmattemp[edges[[i]][1],edges[[i]][2]]
        }
        HBcorrpvals <- p.adjust(uncorrpvals, method='holm')
        el.pvals <- HBcorrpvals
      }
      res <- list(glstr.real = glstr.real, 
                  glstr.perm = glstr.perm,
                  glstr.sep = glstr.sep,
                  glstr.pval = sum(glstr.perm >= glstr.real)/it, 
                  max.real = max.real,
                  max.perm = max.perm,
                  max.pval = sum(max.perm >= max.real)/it, 
                  el.pvals = el.pvals,
                  nw1 = nw1,
                  nw2 = nw2)
    }
    
    if (progressbar==TRUE) close(pb)
    
    if(test.edges==FALSE) 
    {
      res <- list(
        glstr.real = glstr.real, 
        glstr.perm = glstr.perm,
        glstr.sep = glstr.sep,
        glstr.pval = sum(glstr.perm >= glstr.real)/it, 
        max.real = max.real,
        max.perm = max.perm,
        max.pval = sum(max.perm >= max.real)/it, 
        nw1 = nw1,
        nw2 = nw2)
    }
    
  }
  
  #################################
  ### procedure for binary data ###
  #################################
  # Real data
  if(binary.data==TRUE) 
  {
    glstr.perm <- glstr.real <- max.real <- max.perm <- c()
    diffedges.perm <- matrix(NA,it,nedges) 
    edges.pval.HBall <- matrix(NA,nvars,nvars)
    edges.pvalmattemp <- matrix(0,nvars,nvars)
    # real data
    IF1 <- IsingFit(x1, AND = AND, gamma=gamma, plot=FALSE, progressbar=FALSE)
    IF2 <- IsingFit(x2, AND = AND, gamma=gamma, plot=FALSE, progressbar=FALSE)
    nw1 <- IF1$weiadj
    nw2 <- IF2$weiadj
    if(weighted==FALSE){
      nw1=(nw1!=0)*1
      nw2=(nw2!=0)*1
    }
    
    # Difference measures
    # difference in global strength
    glstr.real <- abs(sum(abs(nw1[upper.tri(nw1)]))-sum(abs(nw2[upper.tri(nw2)])))
    # global strength of individual networks
    glstr.sep <- c(sum(abs(nw1[upper.tri(nw1)])), sum(abs(nw2[upper.tri(nw2)])))
    # difference in weights by edge
    diffedges.real <- abs(nw1-nw2)[upper.tri(abs(nw1-nw2))] 
    # maximum norm
    max.real <- max(diffedges.real)
    
    # permuted data
    for (i in 1:it)
    {
      if(paired==FALSE)
      {
        # a check to prevent an error of lognet(): a variable with only one 0 or 1 is not allowed
        checkN=0
        while(checkN==0){
          s <- sample(1:(nobs1+nobs2),nobs1,replace=FALSE)
          x1perm <- dataall[s,]
          x2perm <- dataall[b[-s],]
          
          cm1 <- colMeans(x1perm)
          cm2 <- colMeans(x2perm)
          checkN=ifelse(any(cm1<(1/nobs1))|any(cm1>(nobs1-1)/nobs1)|any(cm2<(1/nobs2))|any(cm2>(nobs2-1)/nobs2),0,1) 
        }
        IF1perm <- IsingFit(x1perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
        IF2perm <- IsingFit(x2perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
        r1perm <- IF1perm$weiadj
        r2perm <- IF2perm$weiadj
        if(weighted==FALSE){
          r1perm=(r1perm!=0)*1
          r2perm=(r2perm!=0)*1
        }
      }
      
      if(paired==TRUE)
      {
        checkN=0
        while(checkN==0)
        {
          s <- sample(c(1,2),nobs1,replace=TRUE)
          x1perm <- x1[s==1,]
          x1perm <- rbind(x1perm,x2[s==2,])
          x2perm <- x2[s==1,]
          x2perm <- rbind(x2perm,x1[s==2,])
          
          cm1 <- colMeans(x1perm)
          cm2 <- colMeans(x2perm)
          checkN=ifelse(any(cm1<(1/nobs1))|any(cm1>(nobs1-1)/nobs1)|any(cm2<(1/nobs2))|any(cm2>(nobs2-1)/nobs2),0,1) 
        }
        
        IF1perm <- IsingFit(x1perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
        IF2perm <- IsingFit(x2perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
        r1perm <- IF1perm$weiadj
        r2perm <- IF2perm$weiadj
        if(weighted==FALSE)
        {
          r1perm=(r1perm!=0)*1
          r2perm=(r2perm!=0)*1
        }
      }
      glstr.perm[i] <- abs(sum(abs(r1perm[upper.tri(r1perm)]))-sum(abs(r2perm[upper.tri(r2perm)])))
      diffedges.perm[i,] <- abs(r1perm-r2perm)[upper.tri(abs(r1perm-r2perm))]
      max.perm[i] <- max(diffedges.perm[i,])
      
      if (progressbar==TRUE) setTxtProgressBar(pb, i)
    }
    
    
    diffedges.realmat <- matrix(diffedges.real,it,nedges,byrow=TRUE)
    # vector with uncorrected p values
    edges.pvaltemp <- colSums(diffedges.perm >= diffedges.realmat)/it 
    # matrix with uncorrected p values
    edges.pvalmattemp[upper.tri(edges.pvalmattemp,diag=FALSE)] <- edges.pvaltemp
    edges.pvalmattemp <- edges.pvalmattemp + t(edges.pvalmattemp)
    
    if(test.edges==TRUE)
    {
      if(is.character(edges))
      {
        ept.HBall <- p.adjust(edges.pvaltemp, method='holm')
        # matrix with corrected p values
        edges.pval.HBall[upper.tri(edges.pval.HBall,diag=FALSE)] <- ept.HBall 
        rownames(edges.pval.HBall) <- colnames(edges.pval.HBall) <- colnames(data1)
        el.pvals <- melt(edges.pval.HBall, na.rm=TRUE, value.name = 'p-value')
      }
      
      if(is.list(edges))
      {
        uncorrpvals <- c()
        for(i in 1:length(edges))
        {
          uncorrpvals[i] <- edges.pvalmattemp[edges[[i]][1],edges[[i]][2]]
        }
        HBcorrpvals <- p.adjust(uncorrpvals, method='holm')
        el.pvals <- HBcorrpvals
      }
      res <- list(glstr.real = glstr.real, 
                  glstr.perm = glstr.perm,
                  glstr.sep = glstr.sep,
                  glstr.pval = sum(glstr.perm >= glstr.real)/it, 
                  max.real = max.real,
                  max.perm = max.perm,
                  max.pval = sum(max.perm >= max.real)/it, 
                  el.pvals = el.pvals,
                  nw1 = nw1,
                  nw2 = nw2)
    }
    
    
    
    if (progressbar==TRUE) close(pb)
    
    if(test.edges==FALSE) 
    {
      res <- list(
        glstr.real = glstr.real, 
        glstr.perm = glstr.perm,
        glstr.sep = glstr.sep,
        glstr.pval = sum(glstr.perm >= glstr.real)/it, 
        max.real = max.real,
        max.perm = max.perm,
        max.pval = sum(max.perm >= max.real)/it, 
        nw1 = nw1,
        nw2 = nw2)
    }
  }
  
  class(res) <- "NCT"
  return(res)
}