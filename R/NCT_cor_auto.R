NCT_cor_auto <- function(data1, 
                         data2, 
                         gamma, 
                         it = 100, 
                         paired=FALSE, 
                         weighted=TRUE, 
                         test.edges=FALSE, 
                         edges, 
                         progressbar=TRUE,
                         make.positive.definite=TRUE,
                         p.adjust.methods="none",
                         model = c("lasso","TMFG")
){ 
  
  #### June 13 2017
  #### Note that the use of cor_auto in combination with NCT is not validated
  #### When treating the data as ordinal, you assume that the data comes from a normal distribution and is thresholded into ordinal categories. When using that in NCT, you make the additional assumption that the thresholding is similar in both groups. It is the question whether this is a plausible assumption. If you do not want to assume this, you might want to treat the data as normal in the ordinary NCT. 
  
  if (missing(gamma)) gamma <- 0.5
  
  if(missing(model)){
      model <- "lasso"
    } else {
      model <- match.arg(model)
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
  #####        Real data          #####
  #####################################
  
    cor_x1 <- cor_auto(x1, verbose=FALSE)
    cor_x2 <- cor_auto(x2, verbose=FALSE)
    
    if(make.positive.definite){
      cor_x1 <- matrix(nearPD(cor_x1, corr=TRUE)$mat, ncol = nvars)
      cor_x1 <- (cor_x1 + t(cor_x1)) / 2 # make symmetric
      cor_x2 <- matrix(nearPD(cor_x2, corr=TRUE)$mat, ncol = nvars)
      cor_x2 <- (cor_x2 + t(cor_x2)) / 2 # make symmetric
    }
    
  if(model=="lasso")
    {
      nw1 <- EBICglasso(cor_x1, nrow(x1),gamma=gamma)
      nw2 <- EBICglasso(cor_x2, nrow(x2),gamma=gamma)
    }else if(model=="TMFG")
    {
      nw1 <- NetworkToolbox::TMFG(cor_x1)$A
      nw2 <- NetworkToolbox::TMFG(cor_x2)$A
    }
    
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
        
        cor_x1 <- cor_auto(x1perm, verbose=FALSE)
        cor_x2 <- cor_auto(x2perm, verbose=FALSE)
        
        if(make.positive.definite){
          cor_x1 <- matrix(nearPD(cor_x1, corr=TRUE)$mat, ncol = nvars)
          cor_x1 <- (cor_x1 + t(cor_x1)) / 2 # make symmetric
          cor_x2 <- matrix(nearPD(cor_x2, corr=TRUE)$mat, ncol = nvars)
          cor_x2 <- (cor_x2 + t(cor_x2)) / 2 # make symmetric
        }
        
        if(model=="lasso")
          {
            r1perm <- EBICglasso(cor_x1, nrow(x1perm),gamma=gamma)
            r2perm <- EBICglasso(cor_x2, nrow(x2perm),gamma=gamma)
          }else if(model=="TMFG")
          {
            r1perm <- NetworkToolbox::TMFG(cor_x1)$A
            r2perm <- NetworkToolbox::TMFG(cor_x2)$A
          }
          
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
        
        cor_x1 <- cor_auto(x1perm, verbose=FALSE)
        cor_x2 <- cor_auto(x2perm, verbose=FALSE)
        
        if(make.positive.definite){
          cor_x1 <- matrix(nearPD(cor_x1, corr=TRUE)$mat, ncol = nvars)
          cor_x1 <- (cor_x1 + t(cor_x1)) / 2 # make symmetric
          cor_x2 <- matrix(nearPD(cor_x2, corr=TRUE)$mat, ncol = nvars)
          cor_x2 <- (cor_x2 + t(cor_x2)) / 2 # make symmetric
        }
        
        if(model=="lasso")
          {
            r1perm <- EBICglasso(cor_x1,nrow(x1perm),gamma=gamma)
            r2perm <- EBICglasso(cor_x2,nrow(x2perm),gamma=gamma)
          }else if(model=="TMFG")
          {
            r1perm <- NetworkToolbox::TMFG(cor_x1)$A
            r2perm <- NetworkToolbox::TMFG(cor_x2)$A
          }
          
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
    #####################################
    #####      End p-value calc     #####
    #####################################
    
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

  class(res) <- "NCT"
  return(res)
}
