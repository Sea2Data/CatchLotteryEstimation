% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap.R
\name{resample}
\alias{resample}
\title{resample data}
\usage{
resample(
  samples,
  hierarchy,
  nSamples = rep(NA, length(hierarchy)),
  replacement = rep(T, length(hierarchy)),
  popSize = rep(NA, length(hierarchy)),
  prefix = ""
)
}
\arguments{
\item{samples}{\code{\link[data.table]{data.table}} with samples}

\item{hierarchy}{character() vector describing the sampling units in hierarchical order}

\item{nSamples}{character() vector corresponding to 'hierarchy', indicating columns for that specify how many samples should be included in the resampling at the corresponding stage. NAs signify that the number of samples that is available in 'samples' should be included.}

\item{replacement}{logical() vector corresponding to 'hierarchy', indicating whether the corresponding sampling units should be resampled with replacement}

\item{popSize}{character() vector corresponding to 'hierarchy', indicating columns for population sizes, NAs signify that the number of samples that is available in 'samples' should be considered to be the population.}

\item{prefix}{character() a string to prefix used to rename sampling units in resampled data}
}
\value{
resampled data, with sampling units renamed. Samples are sorted by the order of selection at each level in the hiearchy (first row is the sample first selected)
}
\description{
Resamples data from hierarchically clustered data,
such as those ariseing from multi-stage sampling with equal selection probabilities at each stage.
}
\details{
The parameter 'samples' should have one row for each set of measurements / observations,
and have have columns identifying the sampling units that observations were sampled from.
The parameter 'hierarchy' enumaretes the columns that identify the sampling units in hierarchical order,
primary sampling unit first.

Equal probability selection is assumed for all stages of sampling.
}
\examples{
 # The data set contains SSUs where no replicates where sampled
 data(NSSH2019)

 # set the number of replicates desired for SSU
 NSSH2019$SSUrepl <- 2

 # Mimicking replicate SSU sampling as selection without replacement at SSU level and
 # selection with replacement as if the fish sampled in an SSU is resampled from the enitre haul
 rs <- resample(NSSH2019, c("PSUid", "SSUid","FishId"),
                          replacement = c(FALSE,FALSE,TRUE),
                          nSamples = c(NA,"SSUrepl",NA),
                          popSize = c(NA,"nSSU",NA))

 # resampling only the SSUs,
 # keeping the PSUs and fish as is, leads to exact duplication of samples
 # since there is only one SSU pr PSU in this example
 rs <- resample(NSSH2019, c("PSUid", "SSUid"),
                          replacement = c(FALSE,TRUE),
                          nSamples = c(NA,"SSUrepl"))
 nrow(rs) == 2*nrow(NSSH2019)
 sum(rs$age) == 2*sum(NSSH2019$age)

 # replicate sampling at SSU level can be mimicked
 # by resampling within SSU.
 # Note that the number of rows are still exactly twice the original samples
 # because the SSUs are resampled within their PSU.
 rs <- resample(NSSH2019, c("PSUid", "SSUid","FishId"),
                          replacement = c(FALSE,FALSE,TRUE),
                          nSamples = c(NA,"SSUrepl",NA), popSize = c(NA,"nSSU",NA))
 nrow(rs) == 2*nrow(NSSH2019)
 sum(rs$age) != 2*sum(NSSH2019$age)

 # resampling without replacement at all levels,
 # and population sizes equal to sample sizes (default)
 # simply reorders the data
 rs <- resample(NSSH2019, c("PSUid", "SSUid", "FishId"),
                          replacement = c(FALSE,FALSE,FALSE))
 nrow(rs) == nrow(NSSH2019)
 sum(rs$age) == sum(NSSH2019$age)

}
