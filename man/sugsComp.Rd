% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sugsComp.R
\name{sugsComp}
\alias{sugsComp}
\title{The computational function for the SUGS with variable selection algorithm.}
\usage{
sugsComp(mydata, intfeature, Model, mu_0 = NULL, lambda_0 = 0.01,
  nu_0 = NULL, S_0 = NULL, betaHat = c(0.01, 0.1, 1, 5, 10, 15, 30, 50,
  100), a = 10, b = 1, w = c(0.5, 0.5))
}
\arguments{
\item{mydata}{Data matrix with observations as rows.}

\item{intfeature}{A binary vector of feature which are parition as irrelevant (0)
or relevant (1).}

\item{Model}{A character string sating whether PML, ML or both are used for feature selection
and returned.}

\item{mu_0}{The mean hyperparameter, default is the column means of the data matrix.}

\item{lambda_0}{The variance of the Guassian mean prior, the dafault value is \code{0.01}.}

\item{nu_0}{The degrees of freedom hyperparameter, the default value is \code{D}, where \code{D} is the number of variables.}

\item{S_0}{The scale hyperparamter, the deault value is a fifth of the column variance of the data matrix.}

\item{betaHat}{A grid of hyperparameters for the dirichlet concentration parameter, the default is \code{c(0.01, 0.1, 1, 5, 10, 15, 30, 50, 100)}.}

\item{a}{The scale of the gamma prior for the dirichlet concentration parameter, the dafault value is \code{10}.}

\item{b}{The rate of the gamma prior for the dirichlet concentration parameter, the default value is \code{1}.}

\item{w}{The prior probability of a variable belong to the irrelevant or relevant partition. The vector must
contain two entries the first entry being the probabiliy of being irreleavnt and the second being the probability of being relevant
The default value is \code{c(0.5,0.5)}.}
}
\value{
An allocation vector, called member, of observation to clusters, K the number of clusters,
n the number of observation in each cluster, the log PML or log ML or both as appropriate, a binary
vector of features that have been classified as irrelevant (0) or relevant (1).
}
\description{
The computational function for the SUGS with variable selection algorithm.
}
