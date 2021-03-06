% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampleStatistics.R
\name{calculateSampleProportionCovariance}
\alias{calculateSampleProportionCovariance}
\title{Covariance of sample proprotions}
\usage{
calculateSampleProportionCovariance(proportions)
}
\arguments{
\item{proportions}{vector of size 'n' containging sampled proportions for a categorical variable with 'n' levels}
}
\value{
matrix() 'n'x'n' matrix with proportion covariances.
}
\description{
Covariance of sample proprotions
}
\examples{
 ages <- c(4,5,5,6,10)
 categories <- 1:20
 sampleCounts <- countCategorical(ages, categories)
 sampleProrportions <- sampleCounts / sum(sampleCounts)
 covariance <- calculateSampleProportionCovariance(sampleProrportions)
}
