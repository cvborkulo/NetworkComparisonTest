% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3Methods.R
\name{plot.NCT}
\alias{plot.NCT}
\alias{print.NCT}
\title{Print method for NCT}
\usage{
\method{print}{NCT}(x, ...)

\method{plot}{NCT}(x, what = c("strength", "network", "edge", "centrality"), ...)
}
\arguments{
\item{x}{output of NCT}

\item{...}{for now ignored}

\item{what}{defines what has to be plotted: results pertaining to test on invariance of global strength ("strength"), network structure ("network"), edge strength ("edge"), or specific centrality measure ("centrality")}
}
\description{
Print method, prints the NCT output, plot method plots the output, summary method returns a summary of the output.
}
