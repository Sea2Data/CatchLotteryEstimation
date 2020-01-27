
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

#
# Some datatype conventions
#

#' Catch at age
#'
#' Estiamte of catch at age in numbers
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

