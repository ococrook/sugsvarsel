% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AddStatsDiag.R
\name{addStatsDiag}
\alias{addStatsDiag}
\title{Function to sequentially update statistics and posterior in SUGS (Gaussian mixtures)}
\usage{
addStatsDiag(x, K, clustnew, n, x_bar, SCL, m, lambda, S, nu, beta, S_0,
  lambda_0, mu_0, nu_0)
}
\arguments{
\item{x}{The observation currently under considetion}

\item{K}{The number of currently occupied clusters}

\item{clustnew}{The new proposed cluster to add new statistics}

\item{n}{The vector indicating the number of observations in each cluster}

\item{x_bar}{The sample mean statistics}

\item{SCL}{The sample variable statistics}

\item{m}{The current posterior mean}

\item{lambda}{The current posterior mean variance}

\item{S}{The current posterior scale vector}

\item{nu}{The current posterior degrees of freedom}

\item{beta}{The current posterior concentration parameter (not currently used in this implementation)}

\item{S_0}{The prior scale vector}

\item{lambda_0}{The prior mean variance}

\item{mu_0}{The prior mean}

\item{nu_0}{The prior degrees of freedom}
}
\value{
Returns the updated statistics and posteriors.
}
\description{
Function to sequentially update statistics and posterior in SUGS (Gaussian mixtures)
}
