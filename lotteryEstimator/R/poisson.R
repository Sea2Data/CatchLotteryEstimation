#' Poisson inclusion probability
#' @description calculates inclusion probability for Poisson sampling, based on selection probability and expected sample size.
#' @param selectionProbability numeric: selection probability for sample
#' @param sampleSize numeric: expected sample size
#' @examples
#'  poissonInclusionProbability(.005,300)
#'  poissonInclusionProbability(c(.0005,0.0001),300)
#'
#'  # selection proportional to fraction of total in sample
#'  # e.g. catch relative to quota
#'  inSample <- 200
#'  total <- 400000
#'  sampleSize <- 300
#'  poissonInclusionProbability(inSample/total,sampleSize)
#'
#' @export
poissonInclusionProbability <- function(selectionProbability, sampleSize){

  if (any(selectionProbability <= 0) | any(selectionProbability > 1)){
    stop("Selection probabilites must be in <0,1].")
  }

  probNotinSample <- (1-selectionProbability)**sampleSize
  return(1-probNotinSample)
}

#' Poisson selection probability
#' @description calculates selection probability for Poisson sampling, based on inclusion probability and expected sample size.
#' @param inclusionProbability numeric: inclusion probability for sample
#' @param sampleSize numeric: expected sample size
#' @examples
#'  # selection proportional to fraction of total in sample
#'  # e.g. catch relative to quota
#'  inSample <- 200
#'  total <- 400000
#'  sampleSize <- 300
#'  inclusionProbability <- poissonInclusionProbability(inSample/total,sampleSize)
#'  selectionProbability <- poissonSelectionProbability(inclusionProbability, sampleSize)
#'
#' @export
poissonSelectionProbability <- function(inclusionProbability, sampleSize){
  if (any(inclusionProbability < 0) | any(inclusionProbability > 1)){
    stop("Inclusion probabilites must be in <0,1].")
  }
  return(1-(1-inclusionProbability)**(1/sampleSize))
}

#' Poission joint inclusion probability
#' @description calculates the pariwise joint inclusion probability for Poisson sampling, based on respective inclusion probabilities
#' @param inclusionProbabilities inclusion probabilities for a sample
#' @return the pairwise joint inclusion probability matrix
#' @export
poissonJointInclusionProbabilityMatrix <- function(inclusionProbabilities){
  coinc <- outer(inclusionProbabilities,inclusionProbabilities)
  diag(coinc) <- NA
  return(coinc)
}
