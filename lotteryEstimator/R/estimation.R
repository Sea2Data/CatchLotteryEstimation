
#' Hansen-Hurwitz estimator
#' @description
#'  Hansen-Hurwitz estimator of the population total for a set of population parameters.
#' @details
#'  Implement estimation from hierarchical sampling by successive application of this estimator.
#' @param sampleTotals list() of numeric() vectors, each corresponding to a sample of the population parameters to be estimated.
#' @param selectionProbabilities numeric() vector of selection probabilites corresponding to the rows in 'sampleTotals'.
#' @return numeric() vector containing the estimated population totals.
#' @export
hansenHurwitz <- function(sampleTotals, selectionProbabilities){
  if (length(sampleTotals) != length(selectionProbabilities)){
    stop("selectionProbabilities does not correspond to the listed sample totals")
  }
  if (sum(selectionProbabilities) > 1 | sum(selectionProbabilities) <= 0){
    stop("sum of selectionProbabilities must be in [0,1>")
  }
  if (any(selectionProbabilities > 1) | any(selectionProbabilities <= 0)){
    stop("all selectionProbabilities must be in [0,1>")
  }

  sumSamples <- sampleTotals[[1]] / selectionProbabilities[1]

  if (length(sampleTotals) > 1){
    for (i in 2:length(sampleTotals)){
      sumSamples <- sumSamples + ( sampleTotals[[i]] / selectionProbabilities[i] )
    }
  }

  return(sumSamples / length(sampleTotals))

}

#' Hierarchical Hansen-Hurwitz
#' @description
#'  \code{\link[lotteryEstimator]{HierarchicalEstimator}} for applying
#'  \code{\link[lotteryEstimator]{hansenHurwitz}} in hierarchical implementations.
#' @param sample \code{\link[data.table]{data.table}} with sample data
#' @param sampleId character() identifying the columns in 'sample' that identify the sampling units to estimate from
#' @param subEstimator function, function for estimating totals for each sample
#' @param selectionProbabilities character() indentifying the columns in 'sample' with selection probabilites for the sampling units.
#' @return numeric() estimate of total
#' @examples
#'  numAtAgeSample <- function(sample){countCategorical(sample$age, 2:20)}
#'  numAtAgeHaul <- function(sample){hierarchicalHorwitzThompson(sample, "SSUid",
#'                                   numAtAgeSample, "SSUinclusionProb")}
#'  exampleSamples <- lotteryEstimator::NSSH2019
#'  exampleSamples$SSUinclusionProb <- exampleSamples$SSUselectionProb * exampleSamples$nSSU
#'  numAtAgeTotal <- hierarchicalHansenHurwitz(exampleSamples, "PSUid",
#'                                             numAtAgeHaul, "PSUselectionProb")
#' @export
hierarchicalHansenHurwitz <- function(sample, sampleId, subEstimator, selectionProbabilities){

  sampleTotals <- list()
  sProb <- list()
  for (id in unique(sample[[sampleId]])){
    sampleUnitData <- sample[sample[[sampleId]] == id,]
    sProb[[id]] <- sampleUnitData[[selectionProbabilities]][1]
    sampleTotals[[id]] <- subEstimator(sampleUnitData)
  }
  return(hansenHurwitz(sampleTotals, unlist(sProb)))

}

#' Hansen-Hurwitz covariance estimator
#' @description
#'  Estimator of the covariance of a Hansen-Hurwitz estimate of a set of population parameters, for single-stage sampling.
#'
#'  This also estimates the first terms of the covariance of a multi-stage Hansen-Hurwitz estimate.
#'  That is the covariance due to variation between sampling units.
#' @details
#'  For multi-stage sampling, this function estimates the covariance due to sampling variation between sampling untis.
#'  The within-unit variance can be estimated with \code{\link[lotteryEstimator]{hansenHurwitzIntra}}.
#'  Combine these estimators for total covariance in hierarchical sampling.
#'
#' @param sampleTotals list() of numeric() vectors, each corresponding to a sample of the population parameters to be estimated.
#' @param selectionProbabilities numeric() vector of selection probabilites for the samples listed in 'sampleTotals'.
#' @return matrix() representing a symmetric matrix with the estimated covariances.
#' @export
hansenHurwitzCovariance <- function(sampleTotals, selectionProbabilities){

  if (length(sampleTotals) != length(selectionProbabilities)){
    stop("selectionProbabilities does not correspond to the listed sample totals")
  }
  if (sum(selectionProbabilities) > 1 | sum(selectionProbabilities) <= 0){
    stop("sum of selectionProbabilities must be in [0,1>")
  }
  if (any(selectionProbabilities > 1) | any(selectionProbabilities <= 0)){
    stop("all selectionProbabilities must be in [0,1>")
  }
  if (length(sampleTotals) < 2){
    stop("estimator not defined for less than two samples.")
  }

  point <- hansenHurwitz(sampleTotals, selectionProbabilities)

  sumSamples <- outer(sampleTotals[[1]], sampleTotals[[1]])*0

  for (i in 1:length(sampleTotals)){
      dev <- (sampleTotals[[i]] / selectionProbabilities[i]) - point
      sumSamples <- sumSamples + outer(dev, dev)
  }


  n <- length(sampleTotals)
  return(sumSamples / (n*(n-1)))

}

#' Intra-unit covariance of Hansen-Hurwitz estimator
#' @description
#'  Estimator of the last term of the covariance of a Hansen-Hurwitz estimate for a multi-stage design.
#'  That is the covariance due to variation within sampling units.
#' @details
#'  This function estimates the covariance due to sampling variation within sampling untis.
#'  The between-unit variance can be estimated with \code{\link[lotteryEstimator]{hansenHurwitzCovariance}}.
#'  Combine these estimators for hierarchical sampling.
#'
#'  For single-stage sampling, the sample totals are known, and there is no within-unit sampling variance.
#' @param sampleCovariances list() of matrix() matrices, with the covariances of sample estimates.
#' @param selectionProbabilities numeric() vector of selection probabilites for the samples whoose covariance are listed in 'sampleCovariances'.
#' @return data.table() representing a symmetric matrix with the estimated covariances.
#' @return matrix() representing a symmetric matrix with the estimated covariances. intra-unit covariances.
#' @export
hansenHurwitzIntra <- function(sampleCovariances, selectionProbabilities){
  if (length(sampleCovariances) != length(selectionProbabilities)){
    stop("selectionProbabilites does not correspond to the listed sample covariances")
  }
  if (sum(selectionProbabilities) > 1 | sum(selectionProbabilities) <= 0){
    stop("sum of selectionProbabilities must be in [0,1>")
  }
  if (any(selectionProbabilities > 1) | any(selectionProbabilities <= 0)){
    stop("all selectionProbabilities must be in [0,1>")
  }

  sumSamples <- sampleCovariances[[1]] / selectionProbabilities[1]**2

  if (length(sampleCovariances) > 1){
    for (i in 2:length(sampleCovariances)){
      sumSamples <- sumSamples + (sampleCovariances[[i]] / selectionProbabilities[i]**2)
    }
  }

  n <- length(sampleCovariances)
  return(sumSamples / n**2)

}

#' Hierarchical Hansen-Hurwitz covariance
#' @description
#'  \code{\link[lotteryEstimator]{HierarchicalCovarianceEstimator}} for estimating
#'  covariance of \code{\link[lotteryEstimator]{hansenHurwitz}} in hierarchical implementations.
#' @param sample \code{\link[data.table]{data.table}} with sample data
#' @param sampleId character() identifying the columns in 'sample' that identify the sampling units to estimate from
#' @param subEstimator function, function for estimating totals for each sample
#' @param subCovarianceEstimator function, function for estimating covariances for each sample
#' @param selectionProbabilities character() indentifying the columns in 'sample' with selection probabilites for the sampling units.
#' @return numeric() estimate of total
#' @examples
#'  numAtAgeSample <- function(sample){countCategorical(sample$age, 2:20)}
#'  numAtAgeHaul <- function(sample){hierarchicalHorwitzThompson(sample, "SSUid",
#'                                   numAtAgeSample, "SSUinclusionProb")}
#'  exampleSamples <- lotteryEstimator::NSSH2019
#'  exampleSamples$SSUinclusionProb <- exampleSamples$SSUselectionProb * exampleSamples$nSSU
#'  covariance <- hierarchicalHansenHurwitzCovariance(exampleSamples, "PSUid", numAtAgeHaul, function(x){0}, "PSUselectionProb")
#' @export
hierarchicalHansenHurwitzCovariance <- function(sample, sampleId, subEstimator, subCovarianceEstimator, selectionProbabilities){
  sampleTotals <- list()
  sampleCovariances <- list()
  sProb <- list()
  for (id in unique(sample[[sampleId]])){
    sampleUnitData <- sample[sample[[sampleId]] == id,]
    sProb[[id]] <- sampleUnitData[[selectionProbabilities]][1]
    sampleTotals[[id]] <- subEstimator(sampleUnitData)
    sampleCovariances[[id]] <- subCovarianceEstimator(sampleUnitData)
  }

  intra <- hansenHurwitzIntra(sampleCovariances, unlist(sProb))
  inter <- hansenHurwitzCovariance(sampleTotals, unlist(sProb))

  return(intra + inter)

}

#' Horvitz-Thompson estimator
#' @description
#'  Horvitz-Thompson estimator of the population total for a set of population parameters.
#' @details
#'  Implement estimation from hierarchical sampling by successive application of this estimator.
#' @param sampleTotals list() of numeric() vectors, each corresponding to a sample of the population parameters to be estimated.
#' @param inclusionProbabilities numeric() vector of inclusion probabilites corresponding to the rows in 'sampleTotals'.
#' @return numeric() vector containing the estimated population totals.
#' @export
horvitzThompson <- function(sampleTotals, inclusionProbabilities){
  if (length(sampleTotals) != length(inclusionProbabilities)){
    stop("inclusionProbabilities does not correspond to the listed sample totals")
  }
  if (any(inclusionProbabilities > 1) | any(inclusionProbabilities <= 0)){
    stop("all inclusionProbabilities must be in [0,1>")
  }

  sumSamples <- sampleTotals[[1]] / inclusionProbabilities[1]

  if (length(sampleTotals) > 1){
    for (i in 2:length(sampleTotals)){
      sumSamples <- sumSamples + ( sampleTotals[[i]] / inclusionProbabilities[i])
    }
  }

  return(sumSamples)

}

#' Hierarchical Howritz Thompson
#' @description
#'  \code{\link[lotteryEstimator]{HierarchicalEstimator}} for applying
#'  \code{\link[lotteryEstimator]{horvitzThompson}} in hierarchical implementations.
#' @details
#'  Example provided in the documentation for \code{\link[lotteryEstimator]{hierarchicalHansenHurwitz}}
#' @param sample \code{\link[data.table]{data.table}} with sample data
#' @param sampleId character() identifying the columns in 'sample' that identify the sampling units to estimate from
#' @param subEstimator function, function for estimating totals for each sample
#' @param inclusionProbabilities character() indentifying the columns in 'sample' with inclusion probabilites for the sampling units.
#' @return numeric() estimate of total
#' @export
hierarchicalHorwitzThompson <- function(sample, sampleId, subEstimator, inclusionProbabilities){
  sampleTotals <- list()
  iProb <- list()
  for (id in unique(sample[[sampleId]])){
    sampleUnitData <- sample[sample[[sampleId]] == id,]
    iProb[[id]] <- sampleUnitData[[inclusionProbabilities]][1]
    sampleTotals[[id]] <- subEstimator(sampleUnitData)
  }
  return(horvitzThompson(sampleTotals, unlist(iProb)))

}

#' Covariance of Horvitz-Thompson estimator
#' @description
#'  Estimator of the covariance of a Horvitz-Thompson estimate for single stage sampling without replacement.
#'  That is the covariance due to variation between sampling units.
#' @details
#'  For multi-stage sampling, this function estimates the covariance due to sampling variation between sampling untis.
#'  The within-unit variance can be estimated with \code{\link[lotteryEstimator]{horvitzThompsonIntra}}.
#'  Combine these estimators for total covariance in hierarchical sampling.
#'
#' @param sampleTotals list() of numeric() vectors, each corresponding to a sample of the population parameters to be estimated.
#' @param inclusionProbabilities numeric() vector of inclusion probabilites for the samples whoose covariance are listed in 'sampleCovariances'.
#' @param coInclusionProbabilities matrix()  of co-inclusion probabilites (joint-inclusionprobabilities) for the samples whoose covariance are listed in 'sampleCovariances'.
#' @return matrix() representing a symmetric matrix with the estimated covariances.
#' @noRd
horvitzThompsonCovariance <- function(sampleTotals, inclusionProbabilities, coInclusionProbabilities){
  if (length(sampleTotals) != length(inclusionProbabilities)){
    stop("inclusionProbabilities does not correspond to the listed sample covariances")
  }
  if (any(inclusionProbabilities > 1) | any(inclusionProbabilities <= 0)){
    stop("all inclusionProbabilities must be in [0,1>")
  }
  warning("Not tested")
  sumSamples <- (1-inclusionProbabilities[1]) * outer(sampleTotals[[1]], sampleTotals[[1]]) / inclusionProbabilities[1]**2

  if (length(sampleTotals) > 1){
    for (i in 2:length(sampleTotals)){
      sumSamples <- sumSamples + (1-inclusionProbabilities[i]) * outer(sampleTotals[[i]], sampleTotals[[i]]) / inclusionProbabilities[i]**2
    }
  }

  for (i in 1:length(sampleTotals)){
    for (j in i:length(sampleTotals)){
      if (i != j){
        sumSamples <- sumSamples + (coInclusionProbabilities[i,j]-inclusionProbabilities[i]*inclusionProbabilities[j]) * outer(sampleTotals[[i]], sampleTotals[[j]]) / (coInclusionProbabilities[i,j] *inclusionProbabilities[i] * inclusionProbabilities[j])
      }
    }
  }

  return(sumSamples)

}

#' Intra-unit covariance of Horvitz-Thompson estimator
#' @description
#'  Estimator of the last term of the covariance of a Horvitz-Thompson estimate for multi-stage sampling without replacement.
#'  That is the covariance due to variation within sampling units.
#' @details
#'  This function estimates the covariance due to sampling variation within sampling untis.
#'  The between-unit variance can be estimated with \code{\link[lotteryEstimator]{horvitzThompsonCovariance}}.
#'  Combine these estimators for hierarchical sampling.
#'
#'  For single-stage sampling, the sample totals are known, and there is no within-unit sampling variance.
#' @param sampleCovariances list() of matrix() matrices, with the covariances of sample estimates.
#' @param inclusionProbabilities numeric() vector of inclusion probabilites for the samples whoose covariance are listed in 'sampleCovariances'.
#' @return data.table() representing a symmetric matrix with the estimated covariances.
#' @return matrix() representing a symmetric matrix with the estimated covariances. intra-unit covariances.
#' @noRd
horvitzThompsonIntra <- function(sampleCovariances, inclusionProbabilities){
  if (length(sampleCovariances) != length(inclusionProbabilities)){
    stop("inclusionProbabilities does not correspond to the listed sample covariances")
  }
  if (any(inclusionProbabilities > 1) | any(inclusionProbabilities <= 0)){
    stop("all inclusionProbabilities must be in [0,1>")
  }

  warning("Not tested")

  sumSamples <- sampleCovariances[[1]] / inclusionProbabilities[1]

  if (length(sampleCovariances) > 1){
    for (i in 2:length(sampleCovariances)){
      sumSamples <- sumSamples + (sampleCovariances[[i]] / inclusionProbabilities[i])
    }
  }

  return(sumSamples)

}

#' Estimate total
#' @description
#'  Estimates total catch at age from estiamtes in strata.
#' @param stratifiedCatchAtAge estimates in each strata, formatted as \code{\link[lotteryEstimator]{stratifiedCatchAtAge}}
#' @return estimates of total catch at age, formatted as \code{\link[lotteryEstimator]{catchAtAge}}
#' @export
estimateFromStrataTotals <- function(stratifiedCatchAtAge){
  estimate <- list()
  estimate$catchAtAge <- Reduce('+', lapply(stratifiedCatchAtAge, FUN=function(x){x$catchAtAge}))
  estimate$covariance <- Reduce('+', lapply(stratifiedCatchAtAge, FUN=function(x){x$covariance}))
  return(estimate)
}


#
# Some datatype conventions
#

#' Hierarchical estimator
#'
#' Function contract for hierarchical estimators.
#' Hierarchical estimators are estimators that depend on other estimates being made for each sampling unit.
#' A hierachical estimator will divide the sample into subsets with data for each sampling unit
#' and call another function to calculate estimates from each of these subsets.
#'
#' These functions take three arguments:
#' \describe{
#'  \item{sample}{\code{\link{data.table}{data.table}} with the sample to estimate from}
#'  \item{sampleId}{character() identifying column in 'sample' encoding the sampling units to estimate from}
#'  \item{subEstimator}{that will provide an estimate for each sampling unit}
#' }
#'
#' @name HierarchicalEstimator
#'
NULL

#' Hierarchical covariance estimator
#'
#' Function contract for hierarchical covariance estimators.
#' Hierarchical covariance estimators are estimators that depend on other estimates being made for each sampling unit.
#' A hierachical covariance estimator will divide the sample into subsets with data for each sampling unit
#' and call another function to calculate estimates from each of these subsets.
#'
#' These functions take three arguments:
#' \describe{
#'  \item{sample}{\code{\link{data.table}{data.table}} with the sample to estimate from}
#'  \item{sampleId}{character() identifying column in 'sample' encoding the sampling units to estimate from}
#'  \item{subEstimator}{that will provide a point estimate for each sampling unit}
#'  \item{subCovarianceEstimator}{that will provide a covariance estimate for each sampling unit}
#' }
#'
#' @name HierarchicalCovarianceEstimator
#'
NULL

#' Catch at age
#'
#' Estimate of catch at age in numbers
#'
#' @details
#'  list with two memebers:
#'  \describe{
#'   \item{catchAtAge}{numeric() named vector with catch at age in numbers for each age group}
#'   \item{covariance}{matrix() with covariance of catch at age between age groups}
#'  }
#'
#' @name catchAtAge
#'
NULL

#' Stratified catch at age
#'
#' Stratified Estimate of catch at age in numbers
#'
#' @details
#'  list with three members:
#'  \describe{
#'     \item{estimates}{named list of \code{\link[lotteryEstimator]{catchAtAge}} estimates, with names identifying strata the estimate was done for.}
#'  }
#'
#'
#' @name stratifiedCatchAtAge
#'
NULL
