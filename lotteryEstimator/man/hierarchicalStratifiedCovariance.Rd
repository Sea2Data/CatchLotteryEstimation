% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.R
\name{hierarchicalStratifiedCovariance}
\alias{hierarchicalStratifiedCovariance}
\title{Hierarchical stratified estimate}
\usage{
hierarchicalStratifiedCovariance(sample, partitionId, subEstimator)
}
\arguments{
\item{sample}{\code{\link[data.table]{data.table}} with sample data}

\item{partitionId}{character() identifying the columns in 'sample' that identify the strata to estimate from}

\item{subEstimator}{\code{\link[lotteryEstimator]{ParameterizedCovarianceEstimator}} for estimating covariances for each strata}
}
\value{
matrix of covariance estimates
}
\description{
\code{\link[lotteryEstimator]{HierarchicalCovarianceEstimator}} for estimating covariances for a stratified sample
}
