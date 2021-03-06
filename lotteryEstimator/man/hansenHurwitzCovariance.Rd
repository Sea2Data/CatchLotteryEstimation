% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.R
\name{hansenHurwitzCovariance}
\alias{hansenHurwitzCovariance}
\title{Hansen-Hurwitz covariance estimator}
\usage{
hansenHurwitzCovariance(sampleTotals, selectionProbabilities)
}
\arguments{
\item{sampleTotals}{list() of numeric() vectors, each corresponding to a sample of the population parameters to be estimated.}

\item{selectionProbabilities}{numeric() vector of selection probabilites for the samples listed in 'sampleTotals'.}
}
\value{
matrix() representing a symmetric matrix with the estimated covariances.
}
\description{
Estimator of the covariance of a Hansen-Hurwitz estimate of a set of population parameters, for single-stage sampling.

This also estimates the first terms of the covariance of a multi-stage Hansen-Hurwitz estimate.
That is the covariance due to variation between sampling units.
}
\details{
For multi-stage sampling, this function estimates the covariance due to sampling variation between sampling untis.
The within-unit variance can be estimated with \code{\link[lotteryEstimator]{hansenHurwitzIntra}}.
Combine these estimators for total covariance in hierarchical sampling.
}
