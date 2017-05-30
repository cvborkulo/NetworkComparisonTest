NCT <- function(data1, data2, gamma, it = 100, binary.data=FALSE, paired=FALSE, weighted=TRUE, AND=TRUE, test.edges=FALSE, edges, 
                progressbar=TRUE, method = c("permute", "bootstrap"), bootcut=c("cutEqual", "none")){ 
  
  if (missing(gamma)){
    if (binary.data){
      gamma <- 0.25
    } else {
      gamma <- 0.5
    }
  }
  
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
  x1_grps <- x1 ## Create a copy of x1
  x1_grps$group <- rep(1, nobs1) ## Add a "group" column to x1
  x2_grps <- x2 ## Create a copy of x2
  x2_grps$group <- rep(2, nobs2) ## Add a "group" column to x2
  dataall_grps <- rbind(x1_grps, x2_grps) ## Create a combined dataset that also contains a "group" column
  
  glstrinv.perm <- glstrinv.real <- nwinv.real <- nwinv.perm <- c()
  diffedges.perm <- matrix(0,it,nedges) 
  diffedges.permtemp <- matrix(0, nvars, nvars)
  einv.perm.all <- array(NA,dim=c(nvars, nvars, it))
  edges.pval.HBall <- matrix(NA,nvars,nvars)
  edges.pvalmattemp <- matrix(0,nvars,nvars)
  
  edgelist <- matrix(NA,nvars,nvars)                          ## The next 7 lines generate a named edgelist vector
  edgelist[upper.tri(edgelist, diag=FALSE)] <- rep(1, nedges)
  rownames(edgelist) <- colnames(edgelist) <- colnames(data1)
  edgelist <- melt(edgelist, na.rm=TRUE, value.name= "Edge Strength")
  edgelistVec <- vector()
  for(i in 1:nedges) {
    edgelistVec[i] <- paste(edgelist[i, 1], edgelist[i, 2], sep = " -- ")
  }
  
  
  if(method == "permute") {
    
    #####################################
    ###    Calculate real values      ###
    #####################################
    ## NONBINARY DATA
    if(binary.data==FALSE) 
    {
      nw1 <- EBICglasso(cor(x1),nrow(x1),gamma=gamma)
      nw2 <- EBICglasso(cor(x2),nrow(x2),gamma=gamma)
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
    }
    
    ## BINARY DATA
    ## Real data
    if(binary.data==TRUE) 
    {
      IF1 <- IsingFit(x1, AND = AND, gamma=gamma, plot=FALSE, progressbar=FALSE)
      IF2 <- IsingFit(x2, AND = AND, gamma=gamma, plot=FALSE, progressbar=FALSE)
      nw1 <- IF1$weiadj
      nw2 <- IF2$weiadj
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
      
    }
    
    #####################################
    ###    Nonbinary Permutation      ###
    #####################################
    
    if(binary.data == FALSE) {
      ## permuted data
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
        
        ## Invariance measures for permuted data
        glstrinv.perm[i] <- abs(sum(abs(r1perm[upper.tri(r1perm)]))-sum(abs(r2perm[upper.tri(r2perm)])))
        diffedges.perm[i,] <- abs(r1perm-r2perm)[upper.tri(abs(r1perm-r2perm))]
        diffedges.permtemp[upper.tri(diffedges.permtemp, diag=FALSE)] <- diffedges.perm[i,]
        diffedges.permtemp <- diffedges.permtemp + t(diffedges.permtemp)
        einv.perm.all[,,i] <- diffedges.permtemp
        nwinv.perm[i] <- max(diffedges.perm[i,])
        
        if (progressbar==TRUE) setTxtProgressBar(pb, i)
      }
      
      if(test.edges==TRUE)
      {
        # vector with uncorrected p values
        edges.pvaltemp <- colSums(diffedges.perm >= diffedges.realmat)/it 
        # matrix with uncorrected p values
        edges.pvalmattemp[upper.tri(edges.pvalmattemp,diag=FALSE)] <- edges.pvaltemp
        edges.pvalmattemp <- edges.pvalmattemp + t(edges.pvalmattemp)
        
        if(is.character(edges))
        {
          ept.HBall <- p.adjust(edges.pvaltemp)
          # matrix with corrected p values
          edges.pval.HBall[upper.tri(edges.pval.HBall,diag=FALSE)] <- ept.HBall 
          rownames(edges.pval.HBall) <- colnames(edges.pval.HBall) <- colnames(data1)
          einv.pvals <- melt(edges.pval.HBall, na.rm=TRUE, value.name = 'p-value')
          einv.perm <- einv.perm.all
          einv.real <- diffedges.realoutput
        }
        
        if(is.list(edges))
        {
          einv.perm <- matrix(NA,it,length(edges))
          colnames(einv.perm) <- edges
          uncorrpvals <- einv.real <- c()
          for(j in 1:length(edges))
          {
            uncorrpvals[j] <- edges.pvalmattemp[edges[[j]][1],edges[[j]][2]]
            einv.real[j] <- diffedges.realoutput[edges[[j]][1],edges[[j]][2]]
            for(l in 1:it){
              einv.perm[l,j] <- einv.perm.all[,,l][edges[[j]][1],edges[[j]][2]]
            }
          }
          HBcorrpvals <- p.adjust(uncorrpvals)
          einv.pvals <- HBcorrpvals
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
                    method = "permute")
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
          method = "permute"
        )
      }
      
      
    }
    
    #################################
    ###    Binary Permutation     ###
    #################################
    if(binary.data == TRUE) {
      
      ## permuted data
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
        
        ## Invariance measures for permuted data
        glstrinv.perm[i] <- abs(sum(abs(r1perm[upper.tri(r1perm)]))-sum(abs(r2perm[upper.tri(r2perm)])))
        diffedges.perm[i,] <- abs(r1perm-r2perm)[upper.tri(abs(r1perm-r2perm))]
        diffedges.permtemp[upper.tri(diffedges.permtemp, diag=FALSE)] <- diffedges.perm[i,]
        diffedges.permtemp <- diffedges.permtemp + t(diffedges.permtemp)
        einv.perm.all[,,i] <- diffedges.permtemp
        nwinv.perm[i] <- max(diffedges.perm[i,])
        
        if (progressbar==TRUE) setTxtProgressBar(pb, i)
      }
      
      if(test.edges==TRUE)
      {
        # vector with uncorrected p values
        edges.pvaltemp <- colSums(diffedges.perm >= diffedges.realmat)/it 
        # matrix with uncorrected p values
        edges.pvalmattemp[upper.tri(edges.pvalmattemp,diag=FALSE)] <- edges.pvaltemp
        edges.pvalmattemp <- edges.pvalmattemp + t(edges.pvalmattemp)
        
        if(is.character(edges))
        {
          ept.HBall <- p.adjust(edges.pvaltemp)
          # matrix with corrected p values
          edges.pval.HBall[upper.tri(edges.pval.HBall,diag=FALSE)] <- ept.HBall 
          rownames(edges.pval.HBall) <- colnames(edges.pval.HBall) <- colnames(data1)
          einv.pvals <- melt(edges.pval.HBall, na.rm=TRUE, value.name = 'p-value')
          einv.perm <- einv.perm.all
          einv.real <- diffedges.realoutput
        }
        
        
        if(is.list(edges))
        {
          einv.perm <- matrix(NA,it,length(edges))
          colnames(einv.perm) <- edges
          uncorrpvals <- einv.real <- c()
          for(j in 1:length(edges))
          {
            uncorrpvals[j] <- edges.pvalmattemp[edges[[j]][1],edges[[j]][2]]
            einv.real[j] <- diffedges.realoutput[edges[[j]][1],edges[[j]][2]]
            for(l in 1:it){
              einv.perm[l,j] <- einv.perm.all[,,l][edges[[j]][1],edges[[j]][2]]
            }
          }
          HBcorrpvals <- p.adjust(uncorrpvals)
          einv.pvals <- HBcorrpvals
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
                    method = "permute")
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
          method = "permute"
        )
      }
    }
    
    class(res) <- "NCT"
    return(res)
    
    
  } #End of permute
  
  if (method == "bootstrap") {
    
    glstrinv <- list()
    nwinv <- list()
    einv <- list()
    
    glstrinv$t <- numeric(0) ## Empty vector to insert bootstrapped values
    nwinv$t <- numeric(0) ## Empy vector to insert bootstrapped values
    einv$t <- data.frame(matrix(NA, it, length(edgelistVec))) ## Empty data frame. Each column is an empty vector representing an edge
    colnames(einv$t) <- edgelistVec ## Name the columns by edge in my empty data frame
    einv$ci <- data.frame(matrix(NA, nedges, 3)) ## Create an empty data frame for confidence intervals of edges
    einv$ci[,1] <- edgelistVec
    colnames(einv$ci) <- c("Edge", "2.5%", "97.5%")
    einv$est <- data.frame(matrix(NA, nedges, 2)) ## Creating an estimate data frame for edges. This will keep output consistent within the bootstrap, rather than using a matrix like in the permutation output
    einv$est[,1] <- edgelistVec
    colnames(einv$est) <- c("Edge", "Estimate")
    
    
    #################################################
    ###            Nonbinary Bootstrap            ###
    #################################################
    if(binary.data==FALSE) {
      
      ##### Run the bootstrap #####
      for(i in 1:it) {
        x <- dataall_grps[sample(nrow(dataall_grps), dim(dataall_grps)[1], replace = TRUE),]      ## Create a bootstrapped sample of the combined data (both x1 and x2 at once). 
        ## Because of the "group" column we added, we will know which data points belong to each group 
        division <- split(x, x$group)                       ## Split the bootstrapped sample by the "group" column 
        div1 <- division[[1]][,ncolall]                     ## Define the x1 half as div1, and cut off the "group" column
        div2 <- division[[2]][,ncolall]                     ## Define the x2 half as div2, and cut off the "group" column
        if(match.arg(bootcut)=="cutEqual") {
          if(dim(div1)[1]<dim(div2)[1]){div2 <- div2[-sample(1:dim(div2)[1], abs(dim(div1)[1]-dim(div2)[1])),]}
          if(dim(div1)[1]>dim(div2)[1]){div1 <- div1[-sample(1:dim(div1)[1], abs(dim(div1)[1]-dim(div2)[1])),]}
        } else {}
        nw1 <- EBICglasso(cor(div1),nrow(div1),gamma=gamma) ## Generate a network for div1
        nw2 <- EBICglasso(cor(div2),nrow(div2),gamma=gamma) ## Generate a network for div2
        if(weighted==FALSE){                                ## If unweighted, reassign values as 0 or 1
          nw1=(nw1!=0)*1
          nw2=(nw2!=0)*1 
        }
        glstrinv$t[i] <- sum(abs(nw1[upper.tri(nw1)]))-sum(abs(nw2[upper.tri(nw2)])) ## Calculate the global strength invariance for this bootstrapped sample
        einv$t[i,] <- (nw1-nw2)[upper.tri((nw1-nw2), diag=FALSE)]  ## Compute a vector of all edge invariances for this sample
        nwinv$t[i] <- max(abs(einv$t[i,])) ## Compute network structure invariance for this sample
        
        if (progressbar==TRUE) setTxtProgressBar(pb, i)
      }
      
      
    } # End of nonbinary
    
    #################################################
    ###             Binary Boostrap               ###
    #################################################
    if(binary.data==TRUE) {
      
      ##### Run the bootstrap #####
      for(i in 1:it) {
        x <- dataall_grps[sample(nrow(dataall_grps), dim(dataall_grps)[1], replace = TRUE),]      ## Create a bootstrapped sample of the combined data (both x1 and x2 at once).
        division <- split(x, x$group)                       ## Split the bootstrapped sample by the "group" column 
        div1 <- division[[1]][,ncolall]                     ## Define the x1 half as div1, and cut off the "group" column
        div2 <- division[[2]][,ncolall]                     ## Define the x2 half as div2, and cut off the "group" column
        if(match.arg(bootcut)=="cutEqual") {
          if(dim(div1)[1]<dim(div2)[1]){div2 <- div2[-sample(1:dim(div2)[1], abs(dim(div1)[1]-dim(div2)[1])),]}
          if(dim(div1)[1]>dim(div2)[1]){div1 <- div1[-sample(1:dim(div1)[1], abs(dim(div1)[1]-dim(div2)[1])),]}
        } else {}
        IF1 <- IsingFit(x1, AND = AND, gamma=gamma, plot=FALSE, progressbar=FALSE) ## Generate networks
        IF2 <- IsingFit(x2, AND = AND, gamma=gamma, plot=FALSE, progressbar=FALSE)
        nw1 <- IF1$weiadj
        nw2 <- IF2$weiadj
        if(weighted==FALSE){                                ## If unweighted, reassign values as 0 or 1
          nw1=(nw1!=0)*1
          nw2=(nw2!=0)*1 
        }
        glstrinv$t[i] <- sum(abs(nw1[upper.tri(nw1)]))-sum(abs(nw2[upper.tri(nw2)])) ## Calculate the global strength invariance for this bootstrapped sample. Note: NOT an absolute value
        einv$t[i,] <- (nw1-nw2)[upper.tri((nw1-nw2), diag=FALSE)]  ## Compute a vector of all edge invariances for this sample. Note: NOT absolute values
        nwinv$t[i] <- max(abs(einv$t[i,])) ## Compute network structure invariance for this sample
        
        if (progressbar==TRUE) setTxtProgressBar(pb, i)
      }
      
      
    } # End of nonbinary
    
    ##### Calculate confidence intervals #####
    
    glstrinv$ci <- quantile(glstrinv$t, probs = c(0.025, 0.975))
    nwinv$ci <- quantile(nwinv$t, probs = c(0.025, 0.975))
    for(i in 1:nedges) {
      einv$ci[i, 2:3] <- quantile(einv$t[,i], probs = c(0.025, 0.975))[1:2]
    }
    
    ##### Calculate estimated true values #####
    
    glstrinv.est <- quantile(glstrinv$t, probs = 0.5)
    nwinv.est <- quantile(nwinv$t, probs = 0.5)
    for(i in 1:nedges) {
      einv$est[i,2] <- quantile(einv$t[,i], probs = 0.5)
    }
    
    ## Printing instructions
    if(test.edges) {
      res <- list()
      res$glstrinv.est <- glstrinv.est
      res$glstrinv.ci <- glstrinv$ci
      res$glstrinv.t <- glstrinv$t
      res$nwinv.est <- nwinv.est
      res$nwinv.ci <- nwinv$ci
      res$nwinv.t <- nwinv$t
      res$edges.tested <- test.edges
      res$einv <- data.frame(cbind(einv$est, einv$ci[,2:3]))
      res$method <- "bootstrap"
    }
    if(test.edges == FALSE) {
      res <- list()
      res$glstrinv.real <- glstrinv.real
      res$glstrinv.ci <- glstrinv$ci
      res$glstrinv.t <- glstrinv$t
      res$nwinv.est <- nwinv.est
      res$nwinv.ci <- nwinv$ci
      res$nwinv.t <- nwinv$t
      res$edges.tested <- test.edges
      res$method <- "bootstrap"
    }
    class(res) <- "NCT"
    return(res)
    
  } #End of bootstrap
  
} # End of NCT function
