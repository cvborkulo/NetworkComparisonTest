\name{NCT-methods}
\alias{plot.NCT}
\alias{summary.NCT}
\alias{print.NCT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for NCT objects
}
\description{
Print method, prints the NCT output, plot method plots the output, summary method returns a summary of the output.
}
\usage{
\method{print}{NCT}(x,\dots)
\method{summary}{NCT}(object,\dots)
\method{plot}{NCT}(x, what = c("strength","network","edge","centrality"),\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
output of \code{\link{NCT}}
}
  \item{object}{
output of \code{\link{NCT}}
}
  \item{what}{
defines what has to be plotted: results pertaining to test on invariance of global strength ("strength"), network structure ("network"), edge strength ("edge"), or specific centrality measure ("centrality")
}
  \item{\dots}{
  Arguments only used in plot method
  }
}

\author{
Claudia D. van Borkulo, with contributions from Jonas Haslbeck, Sacha Epskamp, Payton Jones and Alex Millner

Maintainer: Claudia D. van Borkulo <cvborkulo@gmail.com>
}
