print.NCT <- function(x,...){
  cat("\n NETWORK INVARIANCE TEST \n Test statistic M: ", x$nwinv.real,
      "\n p-value", x$nwinv.pval,
      "\n\n GLOBAL STRENGTH INVARIANCE TEST \n Global strength per group: ", x$glstrinv.sep,
      "\n Test statistic S: ", x$glstrinv.real,
      "\n p-value", x$glstrinv.pval,
      "\n\n EDGE INVARIANCE TEST \n\n")
  print(x$einv.pvals)
  cat("\n CENTRALITY INVARIANCE TEST \n \n")
  print(x$diffcen.pval)
}

summary.NCT <- function(object,...){
  cenTest <- if(is.null(object$diffcen.real)){
    c(NA,NA)
  } else {
    c(reshape2::melt(object$diffcen.real)$value,
      reshape2::melt(object$diffcen.pval)$value)
  }
  
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
      "\n p-value", object$einv.pvals$`p-value`,
      "\n\n CENTRALITY INVARIANCE TEST
      Nodes tested:", rownames(object$diffcen.pval),
      "\n Centralities tested:", colnames(object$diffcen.pval),
      "\n Test statistic C: ", na.omit(cenTest[1]),
      "\n p-value", na.omit(cenTest[2]))
}

plot.NCT <- function(x,what = c("strength","network","edge","centrality"),...){
  
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
    ## Plot results of the difference in edge weight (not reliable with only 10 permutations)
    if(length(dim(x$einv.perm))>2){
      extractLowTri <- function(einv.perm){
        grid <- expand.grid(1:dim(einv.perm)[2], 1:dim(einv.perm)[2])
        grid <- grid[grid[,1]!=grid[,2],]
        grid <- grid[1:(nrow(grid)/2),]
        out <- matrix(NA, nrow=dim(einv.perm)[3], ncol=nrow(grid))
        colnames(out) <- 1:nrow(grid)
        for(i in 1:nrow(grid)){
          out[,i] <- einv.perm[grid[i,1],grid[i,2],]
          colnames(out)[i] <- paste(grid[i,1],grid[i,2],sep="-")
        }
        return(out)
      }
      x$einv.perm <- extractLowTri(x$einv.perm)
    }
    
    nedgetests <- ncol(x$einv.perm)
    for(i in 1:nedgetests){
      hist(x$einv.perm[,i], 
           main=paste('p =',x$einv.pval[i,3]),
           xlab=paste('Difference in edge strength:', 
                      colnames(x$einv.perm)[i]),
           xlim=c(0,max(x$einv.real,x$einv.perm)*1.1))
      points(x$einv.real[i],0,col='red',pch=17)
    }
  }
  
  
  
  if (what == "centrality"){
    ncentests <- ncol(x$diffcen.perm)
    for(i in 1:ncentests){
      hist(x$diffcen.perm[,i], 
           main=paste('p =',reshape2::melt(x$diffcen.pval)$value[i]),
           xlab=paste('Difference in \'', 
                      reshape2::melt(x$diffcen.pval)$Var2[i],
                      '\' for node \'',
                      reshape2::melt(x$diffcen.pval)$Var1[i],
                      '\'',
                      sep=""),
           xlim=c(0,max(x$diffcen.perm[,i],abs(reshape2::melt(x$diffcen.real)$value[i]))))
      points(abs(reshape2::melt(x$diffcen.real)$value[i]),0,col='red',pch=17)
    }
  } 
}

