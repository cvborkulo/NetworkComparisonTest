summary.NCT <- function(x,...){
  if(x$method == "permute") {
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
    )}
  if(x$method == "bootstrap") {
    cat("\n NETWORK INVARIANCE TEST
        estimate: ", x$nwinv.est,
        "\n Confidence interval", x$nwinv.ci,
        "\n\n GLOBAL STRENGTH INVARIANCE TEST
        estimate: ", x$glstrinv.est,
        "\n Confidence interval: ", x$glstrinv.ci,
        "\n\n EDGE INVARIANCE TEST
        Edges tested?: ", x$edges.tested,
        "\n Edge invariance summary: ",
        head(x$edgeinv.mat),
        "\n"
    )
    return(format.data.frame(x$einv, digits = 3))
  }
}

print.NCT <- function(x,...){
  class(x) <- NULL
  if(x$method == "permute") {
    print(x)
  }
  if(x$method=="bootstrap"){
    badItems <- c("glstrinv.t", "edgeinv.t")
    x[which(names(x) %in% badItems)] <- NULL
    x$edgeinv.mat <- head(x$edgeinv.mat)
    cat("\n")
    print(x)
  } 
}

plot.NCT <- function(x,what = c("strength","network","edge"),...){
  
  what <- match.arg(what)
  
  if(x$method == "permute") {
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
  if(x$method == "bootstrap"){
    if (what == "strength"){
      ## Plot results of the global strength invariance test:
      hist(x$glstrinv.t, main = "Global Strength Invariance", 
           xlab = paste('t0=', round(x$glstrinv.est, 3), 'c.i=', round(x$glstrinv.ci[1], 3), 'to', round(x$glstrinv.ci[2], 3)),
           xlim = c(min(x$glstrinv.t), max(x$glstrinv.t)))
      abline(v = x$glstrinv.est, col = "grey", lwd = 2)
      abline(v = x$glstrinv.ci[1], col = "blue", lty = 2, lwd = 2)
      abline(v = x$glstrinv.ci[2], col = "blue", lty = 2, lwd = 2)
    } 
    
    if (what == "edge"){
      print("Method not implemented yet")
      
    }
  }
}