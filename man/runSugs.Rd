% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runSugs.R
\name{runSugs}
\alias{runSugs}
\title{The SUGS clustering algorithm by Wang and Dunson (2011)}
\usage{
runSugs(iter, mydata, Model, mu_0 = NULL, lambda_0 = 0.01, nu_0 = NULL,
  S_0 = NULL, betaHat = c(1, 5, 15, 30, 50, 100), a = 10, b = 1,
  BPPARAM = bpparam())
}
\arguments{
\item{iter}{The number of random orderings of the observations for which to run
sugs.}

\item{mydata}{Data matrix with observations as rows}

\item{Model}{Character string indicating whether to use PML, ML or
both for model selection. "Both" defaults to PML.}

\item{mu_0}{The mean hyperparameter, default is the column means of the data matrix.}

\item{lambda_0}{The variance of the Guassian mean prior, the dafault value is \code{0.01}.}

\item{nu_0}{The degrees of freedom hyperparameter, the default value is \code{2 * (D + 2)}, where \code{D} is the number of variables.}

\item{S_0}{The scale hyperparamter, the deault value is a tenth of the column variance of the data matrix.}

\item{betaHat}{A grid of hyperparameters for the dirichlet concentration parameter, the default is \code{c(1, 5, 15, 30, 50, 100)}.}

\item{a}{The scale of the gamma prior for the dirichlet concentration parameter, the dafault value is \code{10}.}

\item{b}{The rate of the gamma prior for the dirichlet concentration parameter, the default value is \code{1}.}

\item{BPPARAM}{Support for parallel processing using the
\code{BiocParallel} infrastructure. When missing (default),
    the default registered \code{BiocParallelParam} parameters are
    used. Alternatively, one can pass a valid
    \code{BiocParallelParam} parameter instance: \code{SnowParam},
    \code{MulticoreParam}, \code{DoparParam}, \ldots see the
    \code{BiocParallel} package for details. To revert to the
    origianl serial implementation, use \code{serialParam}.}
}
\value{
A matrix of cluster allocation, K the number of clusters, a matrix
indicating the number of observations allocated to each cluster. The value of
the model selection criteria either log PML, log ML or both and the random orderings used
}
\description{
The SUGS clustering algorithm by Wang and Dunson (2011)
}
