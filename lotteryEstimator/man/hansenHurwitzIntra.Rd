% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.R
\name{hansenHurwitzIntra}
\alias{hansenHurwitzIntra}
\title{Intra-unit covariance of Hansen-Hurwitz estimator}
\usage{
hansenHurwitzIntra(sampleCovariances, selectionProbabilities)
}
\arguments{
\item{sampleCovariances}{list() of matrix() matrices, with the covariances of sample estimates.}

\item{selectionProbabilities}{numeric() vector of selection probabilites for the samples whoose covariance are listed in 'sampleCovariances'.}
}
\value{
data.table() representing a symmetric matrix with the estimated covariances.

matrix() representing a symmetric matrix with the estimated covariances. intra-unit covariances.
}
\description{
Estimator of the last term of the covariance of a Hansen-Hurwitz estimate for a multi-stage design.
That is the covariance due to variation within sampling units.
}
\details{
This function estimates the covariance due to sampling variation within sampling untis.
The between-unit variance can be estimated with \code{\link[lotteryEstimator]{hansenHurwitzCovariance}}.
Combine these estimators for hierarchical sampling.

For single-stage sampling, the sample totals are known, and there is no within-unit sampling variance.
}
