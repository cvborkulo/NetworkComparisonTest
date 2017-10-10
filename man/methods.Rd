\name{NCT-methods}
\alias{plot.NCT}
\alias{summary.NCT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for NCT objects
}
\description{
Plot method plots the NCT output, summary method returns a summary of the output.
}
\usage{
\method{plot}{NCT}(x, what = c("strength", "network", "edge"), \dots)
\method{summary}{NCT}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
output of \code{\link{NCT}}
}
\item{what}{
defines what has to be plotted: results pertaining to test on invariance of network structure ("network"), edge strength ("edge"), or global strength ("strength")
}
  \item{object}{
  output of \code{\link{NCT}}
  }
  \item{\dots}{
  Arguments only used in plot method
  }
}

\author{
Claudia D. van Borkulo, with contributions from Sacha Epskamp and Alex Millner

Maintainer: Claudia D. van Borkulo <cvborkulo@gmail.com>
}
