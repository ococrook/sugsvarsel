% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{sugsclusterProbcpp}
\alias{sugsclusterProbcpp}
\title{A C++ accelerated version of sugsclusterProb}
\usage{
sugsclusterProbcpp(x, K, i, D, n, betaHat, phi, m, nu, S, lambda)
}
\arguments{
\item{x}{A numeric vector containing data of the observation}

\item{K}{An integer specifiying the number of clusters}

\item{i}{An integer specifying the current iteration of SUGS algorithm}

\item{D}{An integer specifiying the number of variables}

\item{n}{A numeric vector contain the number of people already allocated to each cluster}

\item{betaHat}{A numeric vector containing the grid of hyperpriors for the dirichlet
concentration parameter}

\item{phi}{A numeric matrix containing the weights for the dirichlet hyperpiors}

\item{m}{A numeric matrix containing the means for each component}

\item{nu}{A numeric vector containg the degrees of freedom for each component}

\item{S}{A numeric matrix containing the scale prior for each component}

\item{lambda}{A numeric vector containing the mean variance hyperparamter for each component}
}
\value{
An arma::vec containg the probability that observation \code{x} belongs to each component
}
\description{
A C++ accelerated version of sugsclusterProb
}
