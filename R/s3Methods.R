summary.NCT <- function(object,...){
  
  cat("\n NETWORK INVARIANCE TEST
 Test statistic M: ", object$nwinv.real,
      "\n p-value", object$nwinv.pval,
      "\n\n GLOBAL STRENGTH INVARIANCE TEST
 Global strength per group: ", object$glstrinv.sep,
      "\n Test statistic S: ", object$glstrinv.real,
      "\n p-value", object$glstrinv.pval,
      "\n\n EDGE INVARIANCE TEST
 Edges tested: ", object$edges.tested,
      "\n Test statistic E: ", object$einv.real,
      "\n p-value", object$einv.pvals
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
      hist(x$einv.perm[,i], main=paste('p =',x$einv.pval[i]),xlab=paste('Difference in edge strength ', colnames(x$einv.perm)[i]),xlim=c(0,max(x$einv.real,x$einv.perm)))
      points(x$einv.real[i],0,col='red',pch=17)
    }

  } #else stop("Method not implemented yet.")
  
}