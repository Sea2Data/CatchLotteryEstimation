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

#' Possion sampling
#' @description
#'  Perform Poisson sampling
#' @details
#'  Sampling is performed on vectors. For sampling a data frame, perform sampling on indecies.
#'
#'  Note the definition of 'sampleSize'. In poisson sampling the sample size is not fixed and may vary from draw to draw.
#'
#' @param population numeric: vector representing the population to sample from
#' @param inclusionProbabilities numeric: vector representing the inclusion probabilites of the populatin
#' @param sampleSize numeric: expected sample size
#' @return numeric a subset of 'population'
#' @examples
#'  data("longlinerPopulation")
#'  codset <- longlinerPopulation[longlinerPopulation$speciesFAO == "COD",]
#'
#'  #select cod catches with inclusion probabilites proportional to their weight
#'  codselection <- codset[poissionSample(1:nrow(codset), codset$wholeWeightKg/sum(codset$wholeWeightKg), 100)]
#'  nrow(codselection)
#' @export
poissionSample <- function(population, inclusionProbabilities, sampleSize){

  if (any(inclusionProbabilities>1) | any(inclusionProbabilities<0)){
    stop("Inclusion probabilites must be in [0,1]")
  }

  rn <- runif(length(population))
  selection <- (rn <= inclusionProbabilities*sampleSize)
  return(population[selection])
}

#' Random sample loss correction
#' @description
#'  Corrects inclusion probabilities due to random loss of samples.
#'  Apprioriate if selected sampling units are not sampled due to random inteference.
#'  The random interference should affect different selections with equal probability
#'  (not related to the fact that the selections may have been made with unequal probability).
#' @param inclusionProbability numeric: vector inclusion probabilites for sampled sampling units.
#' @param selected The number of sampling units that was selected
#' @param sampled The number of sampling units that was actually sampled and recorded
#' @return numeric: vector of corrected inclusion probabilites.
#' @examples
#'  # correct a vector of inclusion probabilities
#'  # 4 samples selected, two missing
#'  randomNonResponseCorrection(c(.4,.1), 4)
#'
#'  # correct a single inclusion probability
#'  # 100 samples selected. 20 missing
#'  randomNonResponseCorrection(.1, 100, 80)
#'
#' @export
randomNonResponseCorrection <- function(inclusionProbabilities, selected, sampled=length(inclusionProbabilities)){
  lost <- selected - sampled
  ex <- 1 - (lost / selected)
  return(1 - (1-inclusionProbabilities)**ex)
}


