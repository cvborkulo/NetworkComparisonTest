#' Print method for NCT
#'
#' @description Print method, prints the NCT output, plot method plots the output, summary method returns a summary of the output.
#'
#' @name plot.NCT
#'
#' @param x output of NCT
#' @param ... additional arguments 
#'
#' @export
#'

print.NCT <- function(x, ...){
  if(x$info$call$abs){
    global_stat_message = "\n\n GLOBAL STRENGTH INVARIANCE TEST \n Global strength per group: "
  } else {
    global_stat_message = "\n\n GLOBAL EXPECTED INFLUENCE INVARIANCE TEST \n Global EI per group: "
  }
  cat("\n NETWORK INVARIANCE TEST \n Test statistic M: \n", x$nwinv.real,
      "\n p-value", x$nwinv.pval,
      global_stat_message, x$glstrinv.sep,
      "\n Test statistic S: ", x$glstrinv.real,
      "\n p-value", x$glstrinv.pval)
  if(x$info$call$test.edges){
    cat("\n\n EDGE INVARIANCE TEST \n")
    print(x$einv.pvals) 
  }
  if(x$info$call$test.centrality){
    cat("\n CENTRALITY INVARIANCE TEST p-value\n")
    print(x$diffcen.pval)
  }
}

#' Summary method for NCT
#'
#' @param object output of NCT
#' @param ... additional arguments 
#'
#' @export
#'
#' 
summary.NCT <- function(object,...){
  if(object$info$call$abs){
    global_stat_message = "\n\n GLOBAL STRENGTH INVARIANCE TEST \n Global strength per group: "
  } else {
    global_stat_message = "\n\n GLOBAL EXPECTED INFLUENCE INVARIANCE TEST \n Global EI per group: "
  }
  
  if(object$info$call$paired){
    paired <- "DEPENDENT GROUPS"
  } else {
    paired <- "INDEPENDENT GROUPS"
  }
  
  if(object$info$call$binary.data){
    data_type <- "BINARY"
  } else {
    data_type <- "GAUSSIAN"
  }
  
  cat("", paired, data_type, "NETWORK COMPARISON TEST \n")
  p_adjust <- ifelse(is.language(object$info$call$p.adjust.methods), "none", object$info$call$p.adjust.methods)
  cat("\n P-VALUE CORRECTION:", p_adjust, "\n")
  
  cat("\n NETWORK INVARIANCE TEST \n Test statistic M:", object$nwinv.real,
      "\n p-value", object$nwinv.pval,
      global_stat_message, object$glstrinv.sep,
      "\n Test statistic S: ", object$glstrinv.real,
      "\n p-value", object$glstrinv.pval)
  if(object$info$call$test.edges){
    cat("\n\n EDGE INVARIANCE TEST \n")
    print(object$einv.pvals) 
  }
  if(object$info$call$test.centrality){
    cenTest <- c(reshape2::melt(object$diffcen.real)$value,
                 reshape2::melt(object$diffcen.pval)$value)  
    cat("\n\n CENTRALITY INVARIANCE TEST \n Nodes tested:", rownames(object$diffcen.pval),
        "\n Centralities tested:", colnames(object$diffcen.pval))
    cat("\n Test statistics C: \n")
    print(object$diffcen.real)
    cat("\n p-values: \n")
    print(object$diffcen.pval)
  }
}

#' Plot method for NCT
#'
#' @param x output of NCT
#' @param what defines what has to be plotted: results pertaining to test on invariance of global strength ("strength"), network structure ("network"), edge strength ("edge"), or specific centrality measure ("centrality")
#' @param ... additional arguments
#'
#' @export
#'
#' 
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

