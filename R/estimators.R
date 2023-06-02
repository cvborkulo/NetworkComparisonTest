# These are the two estimator functions based on code exactly as in the original NCT:
NCT_estimator_Ising <- function(x, gamma = 0.25, AND = TRUE){
  IF <- IsingFit::IsingFit(x, AND = AND, gamma=gamma, plot=FALSE, progressbar=FALSE)
  IF$weiadj
}


NCT_estimator_GGM <- function(x, make.positive.definite = TRUE, gamma = 0.5, corMethod = c("cor","cor_auto"), verbose=FALSE){

  corMethod <- match.arg(corMethod)
  
  if (corMethod == "cor"){
    cor_x <- cor(x)  
  } else if (corMethod == "cor_auto") {
    cor_x <- cor_auto(x, verbose = FALSE)
  }
  
  
  if(make.positive.definite){
    cor_x <- matrix(nearPD(cor_x, corr=TRUE)$mat, ncol = ncol(cor_x))
    cor_x <- (cor_x + t(cor_x)) / 2 # make symmetric
  }
  

  if(verbose){
    nw <- EBICglasso(cor_x,nrow(x),gamma=gamma)
  } else {
    nw <- suppressWarnings(suppressMessages(EBICglasso(cor_x,nrow(x),gamma=gamma)))
  }
  
  return(nw)
}
