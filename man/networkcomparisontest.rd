\name{NetworkComparisonTest}
\alias{NetworkComparisonTest}
\alias{NCT}
\title{
Statistical comparison of two networks based on the difference in strength
}
\description{
This permutation based hypothesis test, suited for gaussian and binary data, assesses the difference in strength (weighted sum of connections) of two networks. Network structures are estimated with l1-regularized partial correlations (gaussian data) or with l1-regularized logistic regression (eLasso, binary data).
}
\usage{
NCT(data1, data2, gamma, it, binary.data, weighted = TRUE, 
AND = TRUE, progressbar = TRUE, ...)
}

\arguments{
  \item{data1}{
One of two datasets. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables. Must be cross-sectional data.
}
  \item{data2}{
The other of two datasets. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables. Must be cross-sectional data.
}
  \item{gamma}{
A single value between 0 and 1 or 'range'. When a single value is provided, networks are estimated with this value for hyperparameter gamma in the extended BIC. When 'range' is provided, networks are estimated across the entire range of gamma's (0, 0.1, 0.2, ..., 1).
}
  \item{it}{
The number of iterations (permutations).
}
  \item{binary.data}{
Logical. Can be TRUE of FALSE to indicate whether the data is binary or not. If binary.data is FALSE, the data is regarded gaussian.
}
  \item{weighted}{
Logical. Can be TRUE of FALSE to indicate whether the networks to be compared should be weigthed of not. If not, the estimated networks are dichotomized. Defaults to TRUE.
}
  \item{AND}{
Logical. Can be TRUE of FALSE to indicate whether the AND-rule or the OR-rule should be used to define the edges in the network. Defaults to TRUE. Only necessary for binary data.
}
  \item{progressbar}{
Logical. Should the pbar be plotted in order to see the progress of the estimation procedure? Defaults to TRUE.
}
 \item{\dots}{ not used. }
}

\value{
NCT returns a 'NCT' object that contains the following items:
\item{diffreal }{The difference in strength between the networks of the observed data sets.}
\item{diffperm }{The difference in strength between the networks of the permutated data sets.}
\item{pval }{The p value resulting from the permutation test.}
}

\references{
Ernst MD. Permutation methods: A basis for exact inference. Stat Sci. 2004;19(4):676-685.

Good PI. Permutation, parametric and bootstrap tests of hypotheses. Vol. 3. New York:: Springer, 2005.

van Borkulo, C. D., Waldorp, L. J., Boschloo, L., Schoevers, R.A., & Borsboom, D. (2015). Distinguishing networks: A permutation test for comparing network structures. Manuscript in preparation. 
}
\author{
Claudia D. van Borkulo

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
# across entire range of gamma. Warning: this takes a while!
# Res_range <- NCT(data1, data2, gamma='range', it=10, binary.data = TRUE)

# Plot results (not reliable with 10 permutations!):
hist(Res_0$diffperm, main=paste('p =',Res_0$pval),xlab='Difference in strength')
points(Res_0$diffreal,col='red')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
