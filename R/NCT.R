NCT <- function(data1, data2, gamma, it, binary.data=FALSE, paired=FALSE, weighted=TRUE, AND=TRUE, progressbar=TRUE, ...){ 
  
  # weighted: are the connection weights taken into account?
  # paired: default is FALSE, so independent samples. May be unequal sample sizes.
  # When TRUE: sample sizes must be equal. Relabeling within each pair of (dependent) observations
  
  if (progressbar==TRUE) pb <- txtProgressBar(max=it, style = 3)
  x1 <- data1
  x2 <- data2
  nobs1 <- nrow(x1)
  nobs2 <- nrow(x2)
  dataall <- rbind(x1,x2)
  data.list <- list(x1,x2)
  b <- 1:(nobs1+nobs2)  
  
  #####################################
  ### procedure for non-binary data ###
  #####################################
  if(binary.data==FALSE) {
    ## test at single gamma
    if(is.numeric(gamma)){
      diffperm <- diffreal <- c()
      # real data
      res1 <- EBICglasso(cor(x1),nrow(x1),gamma=gamma)
      res2 <- EBICglasso(cor(x2),nrow(x2),gamma=gamma)
      if(weighted==FALSE){
        res1=(res1!=0)*1
        res2=(res2!=0)*1
      }
      diffreal <- abs(sum(abs(res1[upper.tri(res1)]))-sum(abs(res2[upper.tri(res2)])))
      # permuted data
      for (i in 1:it){
        if(paired==FALSE){
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
        if(paired==TRUE){
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
        diffperm[i] <- abs(sum(abs(r1perm[upper.tri(r1perm)]))-sum(abs(r2perm[upper.tri(r2perm)])))
        if (progressbar==TRUE) setTxtProgressBar(pb, i)
      }
    }
    
    ### test across whole range of gamma's (0-1)
    if(is.character(gamma)){
      gseq <- seq(0,1,by=.1)
      diffreal <- diffperm <- c()
      diffrealtemp <- matrix(NA,length(gseq),2)
      diffpermtemp <- matrix(NA,it,length(gseq))
      # real data
      for(k in 1:length(gseq)){
        gamma=gseq[k]
        res1 <- EBICglasso(cor(x1),nrow(x1),gamma=gamma)
        res2 <- EBICglasso(cor(x2),nrow(x2),gamma=gamma)
        if(weighted==FALSE){
          res1=(res1!=0)*1
          res2=(res2!=0)*1
        }
        diffrealtemp[k,] <- c(sum(abs(res1[upper.tri(res1)])),sum(abs(res2[upper.tri(res2)])))
      }
      diffreal <- sum(abs(diffrealtemp[,1]-diffrealtemp[,2]))
      # permuted data
      for (i in 1:it){
        if(paired==FALSE){
          s <- sample(1:(nobs1+nobs2),nobs1,replace=FALSE)
          x1perm <- dataall[s,]
          x2perm <- dataall[b[-s],]
          for(j in 1:length(gseq)){
            gamma=gseq[j]
            r1perm <- EBICglasso(cor(x1perm),nrow(x1perm),gamma=gamma)
            r2perm <- EBICglasso(cor(x2perm),nrow(x2perm),gamma=gamma)
            if(weighted==FALSE){
              r1perm=(r1perm!=0)*1
              r2perm=(r2perm!=0)*1
            }
            # make this suited for plotting area under curve
            temp1 <- sum(abs(r1perm[upper.tri(r1perm)]))
            temp2 <- sum(abs(r2perm[upper.tri(r2perm)]))
            diffpermtemp[i,j] <- abs(temp1-temp2)
          }
        }
        if(paired==TRUE){
          s <- sample(c(1,2),nobs1,replace=TRUE)
          x1perm <- x1[s==1,]
          x1perm <- rbind(x1perm,x2[s==2,])
          x2perm <- x2[s==1,]
          x2perm <- rbind(x2perm,x1[s==2,])
          for(j in 1:length(gseq)){
            gamma=gseq[j]
            r1perm <- EBICglasso(cor(x1perm),nrow(x1perm),gamma=gamma)
            r2perm <- EBICglasso(cor(x2perm),nrow(x2perm),gamma=gamma)
            if(weighted==FALSE){
              r1perm=(r1perm!=0)*1
              r2perm=(r2perm!=0)*1
            }
            # make this suited for plotting area under curve
            temp1 <- sum(abs(r1perm[upper.tri(r1perm)]))
            temp2 <- sum(abs(r2perm[upper.tri(r2perm)]))
            diffpermtemp[i,j] <- abs(temp1-temp2)
          }
        }
        if (progressbar==TRUE) setTxtProgressBar(pb, i)
      }
      diffperm <- abs(rowSums(diffpermtemp))
    }
    if (progressbar==TRUE) close(pb)
    res <- list(diffreal = diffreal, 
                diffperm = diffperm,
                pval = sum(diffperm >= diffreal)/it)
  }
  
  #################################
  ### procedure for binary data ###
  #################################
  if(binary.data==TRUE) {
    ## test at single gamma
    if(is.numeric(gamma)){
      diffperm <- diffreal <- c()
      # real data
      IF1 <- IsingFit(x1,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
      IF2 <- IsingFit(x2,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
      res1 <- IF1$weiadj
      res2 <- IF2$weiadj
      if(weighted==FALSE){
        res1=(res1!=0)*1
        res2=(res2!=0)*1
      }
      diffreal <- abs(sum(abs(res1[upper.tri(res1)]))-sum(abs(res2[upper.tri(res2)])))
      # permuted data
      for (i in 1:it){
        if(paired==FALSE){
          s <- sample(1:(nobs1+nobs2),nobs1,replace=FALSE)
          x1perm <- dataall[s,]
          x2perm <- dataall[b[-s],]
          IF1perm <- IsingFit(x1perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
          IF2perm <- IsingFit(x2perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
          r1perm <- IF1perm$weiadj
          r2perm <- IF2perm$weiadj
          if(weighted==FALSE){
            r1perm=(r1perm!=0)*1
            r2perm=(r2perm!=0)*1
          }
        }
        if(paired==TRUE){
          s <- sample(c(1,2),nobs1,replace=TRUE)
          x1perm <- x1[s==1,]
          x1perm <- rbind(x1perm,x2[s==2,])
          x2perm <- x2[s==1,]
          x2perm <- rbind(x2perm,x1[s==2,])
          IF1perm <- IsingFit(x1perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
          IF2perm <- IsingFit(x2perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
          r1perm <- IF1perm$weiadj
          r2perm <- IF2perm$weiadj
          if(weighted==FALSE){
            r1perm=(r1perm!=0)*1
            r2perm=(r2perm!=0)*1
          }
        }
        diffperm[i] <- abs(sum(abs(r1perm[upper.tri(r1perm)]))-sum(abs(r2perm[upper.tri(r2perm)])))
        if (progressbar==TRUE) setTxtProgressBar(pb, i)
      }
    }
    
    ### test across whole range of gamma's (0-1)
    if(is.character(gamma)){
      gseq <- seq(0,1,by=.1)
      diffreal <- diffperm <- c()
      diffrealtemp <- matrix(NA,length(gseq),2)
      diffpermtemp <- matrix(NA,it,length(gseq))
      # real data
      for(k in 1:length(gseq)){
        gamma=gseq[k]
        IF1 <- IsingFit(x1,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
        IF2 <- IsingFit(x2,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
        res1 <- IF1$weiadj
        res2 <- IF2$weiadj
        if(weighted==FALSE){
          res1=(res1!=0)*1
          res2=(res2!=0)*1
        }
        diffrealtemp[k,] <- c(sum(abs(res1[upper.tri(res1)])),sum(abs(res2[upper.tri(res2)])))
      }
      diffreal <- sum(abs(diffrealtemp[,1]-diffrealtemp[,2]))
      # permuted data
      for (i in 1:it){
        if(paired==FALSE){
          s <- sample(1:(nobs1+nobs2),nobs1,replace=FALSE)
          x1perm <- dataall[s,]
          x2perm <- dataall[b[-s],]
          for(j in 1:length(gseq)){
            gamma=gseq[j]
            IF1perm <- IsingFit(x1perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
            IF2perm <- IsingFit(x2perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
            r1perm <- IF1perm$weiadj
            r2perm <- IF2perm$weiadj
            if(weighted==FALSE){
              r1perm=(r1perm!=0)*1
              r2perm=(r2perm!=0)*1
            }
          }
        }
        if(paired==TRUE){
          s <- sample(c(1,2),nobs1,replace=TRUE)
          x1perm <- x1[s==1,]
          x1perm <- rbind(x1perm,x2[s==2,])
          x2perm <- x2[s==1,]
          x2perm <- rbind(x2perm,x1[s==2,])
          for(j in 1:length(gseq)){
            gamma=gseq[j]
            IF1perm <- IsingFit(x1perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
            IF2perm <- IsingFit(x2perm,AND = AND, gamma=gamma,plot=FALSE,progressbar=FALSE)
            r1perm <- IF1perm$weiadj
            r2perm <- IF2perm$weiadj
            if(weighted==FALSE){
              r1perm=(r1perm!=0)*1
              r2perm=(r2perm!=0)*1
            }
          }
        }
        if (progressbar==TRUE) setTxtProgressBar(pb, i)
      }
      diffperm <- abs(rowSums(diffpermtemp))
    }
    if (progressbar==TRUE) close(pb)
    res <- list(diffreal = diffreal, 
                diffperm = diffperm,
                pval = sum(diffperm >= diffreal)/it)
  }
  
  class(res) <- "NCT"
  return(res)
}