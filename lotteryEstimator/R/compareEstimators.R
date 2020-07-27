
#' Sample example population
#' @details Performs two-stage sampling without replacement from an example population
#' @param nPSU the number of primary sampling units to sample
#' @param nSSU the maximal number of secondary sampling untis to sample
#' @param population \code{\link[data.table]{data.table}} with the entire popoulation to sample. The population must be organised in clusters for two stage sampling, identfied by 'PSU' and 'SSU'.
#' @param PSU character() specifying the column in 'population' that identifies the primary sampling units
#' @param SSU character() specifying the column in 'population' that identifies the secondary sampling units
#' @param addSamplingProbs logical() whether sampling probabilites and population totals will be appended
#' @return \code{\link[data.table]{data.table}} formatted as 'population' containing only the sampled elements. If addSamplingProbs is True, the following additional columns will be added:
#'  \describe{
#'   \item{Npsu}{total number of PSUs in the population}
#'   \item{Nssu}{total number of SSUs within the given PSU in the population}
#'   \item{nPSU}{total number of PSUs selected}
#'   \item{nSSU}{total number of SSUs selected withing the given PSU.}
#'   \item{ipPSU}{inclusion probability of PSU}
#'   \item{spPSU}{selection probability of PSU}
#'   \item{ipSSU}{inclusion probability of SSU, given that PSU is selected}
#'   \item{spSSU}{selection probability of SSU, given that PSU is selected}
#'  }
#' @examples
#'  sample <- sampleExamplePopulation(10,10)
#' @export
sampleExamplePopulation <- function(nPSU, nSSU, population=lotteryEstimator::longlinerPopulation, PSU="vesselId", SSU="opid", addSamplingProbs=T){

  if (!data.table::is.data.table(population)){
    stop("population must be a data table.")
  }

  newcols <- c("Npsu", "Nssu", "nPSU", "nSSU", "ipPSU", "spPSU", "ipSSU", "spSSU")
  if (any(newcols %in% names(population))){
    existing <- newcols[newcols %in% names(population)]
    stop(paste("The following columns alread exisits in data table 'population':", paste(existing, collapse=",")))
  }

  #annotate total PSUs and SSUs
  ssuTotal <- population[,.(Nssu=length(unique(get(SSU)))), by=PSU]

  population <- merge(population, ssuTotal)
  population$Npsu <- length(unique(population[[PSU]]))

  #annotate sample sizes
  population$nPSU <- nPSU
  population$nSSU <- pmin(population$Nssu, nSSU)

  #annotate sampling probabilities
  population$ipPSU <- population$nPSU / population$Npsu
  population$ipSSU <- population$nSSU / population$Nssu

  sample <- resample(population, hierarchy=c(PSU, SSU), nSamples=c("nPSU", "nSSU"), replacement=c(F,F))

  PSUselectionOrder <- match(sample[[PSU]], sample[[PSU]][!duplicated(sample[[PSU]])])
  SSUselectionOrder <- match(sample[[SSU]], sample[[SSU]][!duplicated(sample[[SSU]])])
  SSUselectionOrder <- SSUselectionOrder - SSUselectionOrder[!duplicated(sample[[PSU]])][match(sample[[PSU]],sample[[PSU]][!duplicated(sample[[PSU]])])] + 1

  sample$spPSU <- 1 / (sample$Npsu - PSUselectionOrder + 1)
  sample$spSSU <- 1 / (sample$Nssu - SSUselectionOrder + 1)

  if (!addSamplingProbs){
    sample[,newcols] <- NULL
  }

  return(sample)
}

#' Check estimator
#' @description checks an estimator against a know population by repeated sampling and estmation
#' @details The analysis is restricted to 2-stage sampling using \code{\link[lotteryEstimator]{sampleExamplePopulation}}. Parameters 'nPSU', 'nSSU', 'population', 'PSU' and 'SSU' are simply passed to this function.
#' @param estimator function that takes a single parameter (the output from \code{\link[lotteryEstimator]{sampleExamplePopulation}}) and returns an estimate (vector or scalar).
#' @param popParameter The value of the population parameter to compare the estimate to (vector or scalar)
#' @param iterations the number of iterations to run
#' @param nPSU parameter passed to \code{\link[lotteryEstimator]{sampleExamplePopulation}}
#' @param nSSU parameter passed to \code{\link[lotteryEstimator]{sampleExamplePopulation}}
#' @param population parameter passed to \code{\link[lotteryEstimator]{sampleExamplePopulation}}
#' @param PSU parameter passed to \code{\link[lotteryEstimator]{sampleExamplePopulation}}
#' @param SSU parameter passed to \code{\link[lotteryEstimator]{sampleExamplePopulation}}
#' @return list with the following values (vector or scalar):
#'  \describe{
#'   \item{iterations}{The number of estimations run}
#'   \item{simulated.bias}{The difference between the arithmetric mean of estimates and 'popParameter'}
#'   \item{simulated.relative.bias}{simulated.bias divided by 'popParameter'}
#'   \item{mean.sq.error}{The mean squared error of estimates (MSE)}
#'   \item{root.mean.sq.error}{The square root of the mean squared error of estimates (RMSE)}
#'  }
#' @examples
#'  data(longlinerPopulation)
#'
#'  #total COD in example population
#'  total <- sum(longlinerPopulation$wholeWeightKg[longlinerPopulation$speciesFAO=="COD"])
#'
#'  #Naive estimator treat clustered sampling as unclustered
#'  est <- function(sample){
#'     mean(sample$wholeWeightKg[sample$speciesFAO=="COD"])*
#'      sum(longlinerPopulation$speciesFAO=="COD")}
#'
#'  #evaluate naive estimator
#'  \dontrun{checkEstimatorslonglinerPopulation(est, total, 1000, 20, 1)}
#'
#'  #confirm that naive estimator is correct when entire population is sampled
#'  checkEstimatorslonglinerPopulation(est, total, 2, 72, 365)
#'
#' @export
checkEstimatorsExamplePopulation <- function(estimator, popParameter, iterations, nPSU, nSSU, population=lotteryEstimator::longlinerPopulation, PSU="vesselId", SSU="opid"){

  estimates <- list()
  sq.error <- list()
  stopifnot(iterations>0)

  for (i in 1:iterations){
    sample <- sampleExamplePopulation(nPSU, nSSU, population = population, PSU=PSU, SSU=SSU)
    est <- estimator(sample)
    diff <- est - popParameter
    estimates[[i]] <- est
    sq.error[[i]] <- diff*diff
  }

  result <- list()
  result$iterations <- iterations
  result$simulated.bias <- (Reduce("+", estimates) / iterations ) - popParameter
  result$simulated.relative.bias <- result$simulated.bias / popParameter
  result$mean.sq.error <- (Reduce("+", sq.error) / iterations )
  result$root.mean.sq.error <- sqrt(result$mean.sq.error)

  return(result)
}
