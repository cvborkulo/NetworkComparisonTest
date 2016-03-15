\name{NetworkComparisonTest}
\alias{NetworkComparisonTest}
\alias{NCT}
\title{
Statistical comparison of two networks based on the difference in strength
}
\description{
This permutation based hypothesis test, suited for gaussian and binary data, assesses the difference between two networks based on several measures (difference in global strength, the largest difference between edges, individual edges). Network structures are estimated with l1-regularized partial correlations (gaussian data) or with l1-regularized logistic regression (eLasso, binary data). Suited for comparison of independent and dependent samples (currently, only for one group measured twice).
}
\usage{
NCT(data1, data2, gamma, it, binary.data=FALSE, paired=FALSE, weighted=TRUE, AND=TRUE, test.edges=FALSE, edges, progressbar=TRUE, ...)
}

\arguments{
  \item{data1}{
One of two datasets. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables. Must be cross-sectional data.
}
  \item{data2}{
The other of two datasets. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables. Must be cross-sectional data.
}
  \item{gamma}{
A single value between 0 and 1. Networks are estimated with this value for hyperparameter gamma in the extended BIC.
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
Logical. Can be TRUE of FALSE to indicate whether the networks to be compared should be weigthed of not. If not, the estimated networks are dichotomized. Defaults to TRUE.
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
 \item{\dots}{ not used. }
}

\value{
NCT returns a 'NCT' object that contains the following items:
\item{glstr.real }{The difference in global strength between the networks of the observed data sets.}
\item{glstr.perm }{The difference in global strength between the networks of the permutated data sets.}
\item{glstr.sep}{The global strength values of the individual networks}
\item{glstr.pval }{The p value resulting from the permutation test concerning difference in global strength.}
\item{max.real}{The value of the maximum difference in edge weights of the observed networks}
\item{max.perm}{The values of the maximum difference in edge weights of the permuted networks}
\item{max.pval }{The p value resulting from the permutation test concerning the maximum difference in edge weights.}
\item{el.pvals}{The Holm-Bonferroni corrected p values per edge from the permutation test concerning differences in edges weights. Only if test.edges = TRUE.}
\item{nw1}{The weighted adjacency matrix of the observed network of data1}
\item{nw2}{The weighted adjacency matrix of the observed network of data2}
}

\references{
Ernst MD. Permutation methods: A basis for exact inference. Stat Sci. 2004;19(4):676-685.

Good PI. Permutation, parametric and bootstrap tests of hypotheses. Vol. 3. New York:: Springer, 2005.

van Borkulo, C. D., Boschloo, L., Borsboom, D., Penninx, B. W. J. H., Waldorp, L. J., & Schoevers, R.A. (2015). Association of symptom network structure with the course of depression. In press.

van Borkulo, C. D., Waldorp, L. J., Boschloo, L., Schoevers, R.A., & Borsboom, D. (2015). Distinguishing networks: A permutation test for comparing network structures. Manuscript in preparation. 
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

### Simulate binary datasets under null hypothesis:
### underlying network structures have the same strength 
# Input:
N <- 6 # Number of nodes
nSample <- 1000 # Number of samples

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
Res_0 <- NCT(data1, data2, gamma=0, it=10, binary.data = TRUE)

# Plot results of global strength (not reliable with only 10 permutations!):
hist(Res_0$glstr.perm, main=paste('p =',Res_0$glstr.pval),xlab='Difference in global strength')
points(Res_0$glstr.real,col='red')

# Plot results of the maximum difference in edge weights (not reliable with only 10 permutations!):
hist(Res_0$max.perm, main=paste('p =',Res_0$max.pval),xlab='Maximum')
points(Res_0$max.real,col='red')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
