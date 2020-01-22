
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
      sumSamples <- sumSamples + ( sampleTotals[[i]] / selectionProbabilities[i])
    }
  }

  return(sumSamples / length(sampleTotals))

}

#' Inter-unit covariance of Hansen-Hurwitz estimator
#' @description
#'  Estimator of the first terms of the covariance of a Hansen-Hurwitz estimate.
#'  That is the covariance due to variation between sampling units.
#' @details
#'  This function estimates the covariance due to sampling variation between sampling untis.
#'  The within-unit variance can be estimated with \code{\link[lotteryEstimator]{HansenHurwitzIntra}}.
#'  Combine these estimators for hierarchical sampling.
#'
#'  For single-stage sampling, the sample totals are known, and there is no within-unit sampling variance.
#' @param sampleTotals list() of numeric() vectors, each corresponding to a sample of the population parameters to be estimated.
#' @param selectionProbabilities numeric() vector of selection probabilites for the samples listed in 'sampleTotals'.
#' @return matrix() representing a symmetric matrix with the estimated inter-unit covariances.
#' @export
hansenHurwitzInter <- function(sampleTotals, selectionProbabilities){

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
  deviations <- lapply(sampleTotals, FUN=function(x){ x - point})

  sumSamples <- (deviations[[1]] %*% t(deviations[[1]])) / selectionProbabilities[1]

  if (length(deviations) > 1){
    for (i in 2:length(deviations)){
      sumSamples <- sumSamples + ( (deviations[[i]] %*% t(deviations[[i]])) / selectionProbabilities[i] )
    }
  }

  n <- length(deviations)
  return(sumSamples / (n*(n-1)))

}

#' Intra-unit covariance of Hansen-Hurwitz estimator
#' @description
#'  Estimator of the last term of the covariance of a Hansen-Hurwitz estimate.
#'  That is the covariance due to variation within sampling units.
#' @details
#'  This function estimates the covariance due to sampling variation within sampling untis.
#'  The between-unit variance can be estimated with \code{\link[lotteryEstimator]{HansenHurwitzInter}}.
#'  Combine these estimators for hierarchical sampling.
#'
#'  For single-stage sampling, the sample totals are known, and there is no within-unit sampling variance.
#' @param sampleCovariances list() of matrix() matrices, with the covariances of sample estimates.
#' @param selectionProbabilities numeric() vector of selection probabilites for the samples whoose covariance are listed in 'sampleCovariances'.
#' @return data.table() representing a symmetric matrix with the estimated covariances.
#' @return matrix() representing a symmetric matrix with the estimated covariances. intra-unit covariances.
#' @export
hansenHurwitzIntra <- function(sampleCovariances, selectionProbabilities){
  if (length(sampleTotals) != length(selectionProbabilites)){
    stop("selectionProbabilites does not correspond to the listed sample covariances")
  }
  if (sum(selectionProbabilities) > 1 | sum(selectionProbabilities) <= 0){
    stop("sum of selectionProbabilities must be in [0,1>")
  }
  if (any(selectionProbabilities > 1) | any(selectionProbabilities <= 0)){
    stop("all selectionProbabilities must be in [0,1>")
  }

}




