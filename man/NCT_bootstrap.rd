\name{NCT_bootstrap}
\alias{NetworkComparisonTest_bootstrap}
\title{
Statistical Comparison of Two Networks Using Bootstrapping
}
\description{
This bootstrapping methodology, suited for gaussian and binary data, assesses the effect size of the difference between two networks based on two invariance measures (global strength invariance, edge invariance). Network structures are estimated with Pearson correlations, partial correlations, l1-regularized partial correlations, or with l1-regularized logistic regression (eLasso, binary data). Suited for comparison of independent and dependent samples.
}
\usage{
NCT_bootstrap <- function(data1, data2, nBoots = 500, default=c("association", "concentration", "EBICglasso", "IsingFit", "custom"), 
                          paired=FALSE, weighted=TRUE, progressbar=TRUE, 
                          bootcut=c("none", "cutEqual"), custom_func=NULL, AND=TRUE)
}

\arguments{
  \item{data1}{
One of two datasets. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables. Must be cross-sectional data.
}
  \item{data2}{
The other of two datasets. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables. Must be cross-sectional data.
}
  \item{nBoots}{
The number of bootstrapped samples to draw upon.
}
  \item{default}{
Which method should be used to calculate networks? If "custom" is used, a "custom_func" argument must be specified
}
  \item{paired}{
Logical. Can be TRUE of FALSE to indicate whether the samples are dependent or not. If paired is TRUE, relabeling is performed within each pair of observations. If paired is FALSE, relabeling is not restricted to pairs of observations. Note that, currently, dependent data is assumed to entail one group measured twice.
}
  \item{weighted}{
Logical. Can be TRUE of FALSE to indicate whether the networks to be compared should be weighted of not. If not, the estimated networks are dichotomized. Defaults to TRUE.
}
  \item{progressbar}{
Logical. Should the pbar be plotted in order to see the progress of the estimation procedure? Defaults to TRUE.
}
  \item{bootcut}{
If "none" is specified, each group is resampled according to their own sample size. 
If "cutEqual" is specified, each group is resampled according to the sample size of the smaller group.
}
  \item{custom_func}{
A custom function if default is set to "custom". Should take participant data (e.g., x1) as the only input, and should output an adjacency matrix.
}
}

\value{
NCT_bootstrap returns a 'NCT' object that contains the following items:
\item{glstrinv.real }{The difference in global strength between the networks of the observed data sets.}
\item{glstrinv.sep}{The global strength values of the individual networks}
\item{glstrinv.est}{The estimated difference in global strength according to the bootstrapped samples}
\item{glstrinv.ci}{The 95\% confidence interval according to the bootstrapped samples}
\item{glstrinv.t}{A vector containing the global strength invariance in each bootstrapped sample}
\item{edgeinv.mat}{A matrix which contains an edgelist with the following columns: edge value in first network, edge value in second network,
edge invariance in observed datasets, estimated edge invariance based on bootstrapping, confidence interval based on bootstrapping}
\item{edgeinv.t}{A dataframe. Each column represents an edge in the network, and each row represents the invariance of that edge for
a bootstrapped sample}
\item{method}{NCT method: "perm" or "bootstrap"}

}

\references{
}
\author{
Payton J. Jones

Maintainers: Payton J. Jones <payton_jones@g.harvard.edu>
Claudia D. van Borkulo <cvborkulo@gmail.com>
}
\note{
See also the website: http://cvborkulo.com
}

\examples{
library("networktools")

# Use the two random halves of the "depression" dataset
s1 <- sample(1:1000, 500)
s2 <- c(1:1000)[-s1]
data1 <- depression[s1,]
data2 <- depression[s2,]

Res_1 <- NCT_bootstrap(data1, data2, nBoots=10, default="IsingFit")

## Plotting of NCT results
## See the help file of plot.NCT for more information about the plotting function and its arguments

# Plot results of global strength invariance test (not reliable with only 10 permutations!):
plot(Res_1, what="strength")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
