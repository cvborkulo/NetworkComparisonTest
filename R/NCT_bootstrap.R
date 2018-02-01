NCT_bootstrap <- function(data1, data2, nBoots = 500, 
                          default=c("association", "concentration", "EBICglasso", "IsingFit", "custom"), 
                          paired=FALSE, weighted=TRUE, progressbar=TRUE, 
                          bootcut=c("none", "cutEqual"), custom_func, AND=TRUE, global.strength=c("raw", "absolute_value")){
  
  if (progressbar==TRUE) pb <- txtProgressBar(max=nBoots, style = 3)
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
  if(match.arg(global.strength) == "aboslute_value"){
    gl.str <- function(x){
      return(abs(sum(abs(x[upper.tri(x)]))))
    } 
  } else if (match.arg(global.strength) == "raw") {
    gl.str <- function(x){
      return(sum(x[upper.tri(x)]))
    } 
  }

  ## The next 7 lines generate an empty edgelist  
  ## We will store edge invariances in edgelist format, for ease with CIs
  edgelist <- matrix(NA,nvars,nvars)                          
  edgelist[upper.tri(edgelist, diag=FALSE)] <- rep(1, nedges)
  rownames(edgelist) <- colnames(edgelist) <- colnames(data1)
  edgelist <- reshape2::melt(edgelist, na.rm=TRUE, value.name= "Edge Strength")
  edgelistVec <- vector()
  for(i in 1:nedges) {
    edgelistVec[i] <- paste(edgelist[i, 1], edgelist[i, 2], sep = " -- ")
  }
  
  glstrinv <- list()
  einv <- list()
  
  glstrinv$t <- numeric(0) ## Empty vector to insert bootstrapped values
  einv$t <- data.frame(matrix(NA, nBoots, length(edgelistVec))) ## Empty data frame. Each column is an empty vector representing an edge
  colnames(einv$t) <- edgelistVec ## Name the columns by edge in my empty data frame
  einv$boot <- data.frame(matrix(NA, nedges, 7)) ## Create an empty data frame for confidence intervals of edges
  einv$boot[,1] <- edgelistVec
  colnames(einv$boot) <- c("Edge","Network1", "Network2", "RealInv", "Estimate", "2.5% CI", "97.5% CI")
  
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
  
  ##### Calculate real values #####
  nw1.real <- fun(x1) 
  nw2.real <- fun(x2)
  glstrinv.real <- gl.str(nw1.real) - gl.str(nw2.real)
  glstrinv.sep <- c(gl.str(nw1.real), gl.str(nw2.real))
  diffedges.real <- abs(nw1.real - nw2.real)[upper.tri(abs(nw1.real - nw2.real))]
  diffedges.realmat <- matrix(diffedges.real, nBoots, nedges, 
                              byrow = TRUE)
  diffedges.realoutput <- abs(nw1.real - nw2.real)
  einv$boot[,2] <- (nw1.real)[upper.tri((nw1.real), diag=FALSE)]
  einv$boot[,3] <- (nw1.real)[upper.tri((nw1.real), diag=FALSE)]
  einv$boot[,4] <- (nw1.real-nw2.real)[upper.tri((nw1.real-nw2.real), diag=FALSE)]
  
  ##### Run the bootstrap #####
  for(i in 1:nBoots) {
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

    if (progressbar==TRUE) setTxtProgressBar(pb, i)
  }
  
  ##### Calculate confidence intervals #####
  
  glstrinv$ci <- quantile(glstrinv$t, probs = c(0.025, 0.975))
  for(i in 1:nedges) {
    einv$boot[i, 6:7] <- quantile(einv$t[,i], probs = c(0.025, 0.975))[1:2]
  }
  
  ##### Calculate estimated true values #####
  
  glstrinv.est <- quantile(glstrinv$t, probs = 0.5)
  for(i in 1:nedges) {
    einv$boot[i, 5] <- quantile(einv$t[,i], probs = 0.5)
  }
  

  ## Return instructions
  res <- list()
  res$glstrinv.real <- glstrinv.real
  res$glstrinv.sep <- glstrinv.sep
  res$glstrinv.est <- glstrinv.est
  res$glstrinv.ci <- glstrinv$ci
  res$glstrinv.t <- glstrinv$t
  res$edgeinv.mat <- einv$boot
  res$edgeinv.t <- einv$t
  res$method <- "bootstrap"
  
  class(res) <- "NCT"
  return(res)
}

