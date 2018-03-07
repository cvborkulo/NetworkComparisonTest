\name{NetworkComparisonTest}
\alias{NetworkComparisonTest}
\alias{NCT}
\title{
Statistical Comparison of Two Networks Based on Three Invariance Measures
}
\description{
This permutation based hypothesis test, suited for gaussian and binary data, assesses the difference between two networks based on several invariance measures (network structure invariance, global strength invariance, edge invariance). Network structures are estimated with l1-regularized partial correlations (gaussian data) or with l1-regularized logistic regression (eLasso, binary data). Suited for comparison of independent and dependent samples (currently, only for one group measured twice).
}
\usage{
NCT(data1, data2, gamma, it, binary.data=FALSE, paired=FALSE, 
    weighted=TRUE, AND=TRUE, test.edges=FALSE, edges, 
    progressbar=TRUE, make.positive.definite=TRUE)
}

\arguments{
  \item{data1}{
One of two datasets. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables. Must be cross-sectional data.
}
  \item{data2}{
The other of two datasets. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables. Must be cross-sectional data.
}
  \item{gamma}{
A single value between 0 and 1. When not entered, gamma is set to 0.25 for binary data and 0.50 for gaussian data. Networks are estimated with this value for hyperparameter gamma in the extended BIC.
}
  \item{it}{
The number of iterations (permutations).
}
  \item{binary.data}{
Logical. Can be TRUE or FALSE to indicate whether the data is binary or not. If binary.data is FALSE, the data is regarded gaussian.
}
  \item{paired}{
Logical. Can be TRUE of FALSE to indicate whether the samples are dependent or not. If paired is TRUE, relabeling is performed within each pair of observations. If paired is FALSE, relabeling is not restricted to pairs of observations. Note that, currently, dependent data is assumed to entail one group measured twice.
}
  \item{weighted}{
Logical. Can be TRUE of FALSE to indicate whether the networks to be compared should be weighted of not. If not, the estimated networks are dichotomized. Defaults to TRUE.
}
  \item{AND}{
Logical. Can be TRUE of FALSE to indicate whether the AND-rule or the OR-rule should be used to define the edges in the network. Defaults to TRUE. Only necessary for binary data.
}
  \item{test.edges}{
Logical. Can be TRUE of FALSE to indicate whether or not differences in individual edges should be tested.
}
  \item{edges}{
Character or list. When 'all', differences between all individual edges are tested. When provided a list with one or more pairs of indices referring to variables, the provided edges are tested. A Holm-Bonferroni correction is applied to control for multiple testing.
}
  \item{progressbar}{
Logical. Should the pbar be plotted in order to see the progress of the estimation procedure? Defaults to TRUE.
}
  \item{make.positive.definite}{
If \code{make.positive.definite = TRUE}, the covariance matrices used for the glasso are projected to the nearest positive definite matrices, if they are not yet positive definite. This is useful for small n, for which it is very likely that at least one of the bootstrap comparisons involves a covariance matrix that is not positive definite.}


}



\value{
NCT returns a 'NCT' object that contains the following items:
\item{glstrinv.real }{The difference in global strength between the networks of the observed data sets.}
\item{glstrinv.perm }{The difference in global strength between the networks of the permutated data sets.}
\item{glstrinv.sep}{The global strength values of the individual networks}
\item{glstrinv.pval }{The p value resulting from the permutation test concerning difference in global strength.}
\item{nwinv.real}{The value of the maximum difference in edge weights of the observed networks}
\item{nwinv.perm}{The values of the maximum difference in edge weights of the permuted networks}
\item{nwinv.pval }{The p value resulting from the permutation test concerning the maximum difference in edge weights.}
\item{einv.pvals}{The Holm-Bonferroni corrected p values per edge from the permutation test concerning differences in edges weights. Only if test.edges = TRUE.}
\item{edges.tested}{The pairs of variables between which the edges are called to be tested. Only if test.edges = TRUE.}
\item{einv.real}{The value of the difference in edge weight of the observed networks (multiple values if more edges are called to test). Only if test.edges = TRUE.}
\item{einv.perm}{The values of the difference in edge weight of the permuted networks. Only if test.edges = TRUE.}
}

\references{
Ernst, M.D. Permutation methods: A basis for exact inference. Stat Sci. 2004;19(4):676-685.

Good, P.I. Permutation, parametric and bootstrap tests of hypotheses. Vol. 3. New York:: Springer, 2005.

van Borkulo, C. D., Boschloo, L., Borsboom, D., Penninx, B. W. J. H., Waldorp, L. J., & Schoevers, R.A. (2015). Association of symptom network structure with the course of depression. JAMA Psychiatry. 2015;72(12). doi:10.1001/jamapsychiatry.2015.2079

van Borkulo, C. D., Boschloo, Kossakowski, J., Tio, P., L., Schoevers, R.A., Borsboom, D., & , Waldorp, L. J. (2016). Comparing network structures on three aspects: A permutation test. Manuscript submitted for publication. 
}
\author{
Claudia D. van Borkulo, with contributions from Sacha Epskamp and Alex Millner

Maintainer: Claudia D. van Borkulo <cvborkulo@gmail.com>
}
\note{
See also my website: http://cvborkulo.com
}

\examples{
library("IsingSampler")
library("IsingFit")

### Simulate binary datasets under null hypothesis:
### underlying network structures have the same strength 
# Input:
N <- 6 # Number of nodes
nSample <- 500 # Number of samples

# Ising parameters:
Graph <- matrix(sample(0:1,N^2,TRUE,prob = c(0.8, 0.2)),N,N) * runif(N^2,0.5,2)
Graph <- pmax(Graph,t(Graph))
diag(Graph) <- 0
Thresh <- -rowSums(Graph) / 2

# Simulate:
data1 <- IsingSampler(nSample, Graph, Thresh)
data2 <- IsingSampler(nSample, Graph, Thresh)

### Compare networks of data sets using NCT ###
# with gamma = 0. Iterations set to 10 to save time. Should be 1000 at least.

# Testing on all three aspects
# 2 edges are tested here: between variable 1 and 2, 
# and between 3 and 6 (can be list(c(2,1),c(6,3)) as well)
Res_1 <- NCT(data1, data2, gamma=0, it=10, binary.data = TRUE, 
test.edges=TRUE, edges=list(c(1,2),c(3,6)))

## Plotting of NCT results
## See the help file of plot.NCT for more information about the plotting function and its arguments

# Plot results of the network structure invariance test (not reliable with only 10 permutations!):
plot(Res_1, what="network")

# Plot results of global strength invariance test (not reliable with only 10 permutations!):
plot(Res_1, what="strength")

# Plot results of the edge invariance test (not reliable with only 10 permutations!):
# Note that two distributions are plotted
plot(Res_1, what="edge")

# Without testing for (an) individual edge(s)
# The arguments 'test.edges' and 'edges' don't need to be specified
# Not run
# Res_0 <- NCT(data1, data2, gamma=0, it=10, binary.data = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
