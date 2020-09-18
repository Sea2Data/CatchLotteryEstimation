#' Random sample loss correction
#' @description
#'  Calulcates factor for correcting inclusion probabilities due to random loss of samples.
#'  Appropriate if selected sampling units are not sampled due to random inteference.
#'  The random interference should affect different selections with equal probability
#'  (not related to the fact that the selections may have been made with unequal probability).
#' @param selectedInclusionProbabilities numeric: vector with inclusion probabilites for selected sampling units.
#' @param sampledInclusionProbabilities numeric: vector with inclusion probabilites for sampled sampling units ('selectedInclusionProbabilities' excluding non-response).
#' @return scalar for correcting inclusionprobabilites
#' @examples
#'  # correct a vector of inclusion probabilities
#'  # 4 samples selected, two missing
#'  inclusionProbabilities <- c(.4,.1,.2,.1)
#'  nonresponse <- c(F,F,T,T)
#'  sampled <- inclusionProbabilities[!nonresponse]
#'  correctionFactor <- randomNonResponseCorrectionFactor(inclusionProbabilities, sampled)
#'  correctedInclusionProbabilites <- sampled * correctionFactor
#'
#'  sampleTotals <- c(3,4)
#'
#'  #estimate
#'  preEst <- horvitzThompson(sampleTotals, correctedInclusionProbabilites)
#'
#'  #correct post-estimate
#'  postEst <- horvitzThompson(sampleTotals, inclusionProbabilities[!nonresponse])/correctionFactor
#'
#'  preEst - postEst
#'
#' @export
randomNonResponseCorrectionFactor <- function(selectedInclusionProbabilities, sampledInclusionProbabilities){
  correction <- sum(1/sampledInclusionProbabilities) / (sum(1/selectedInclusionProbabilities))
  return(correction)

}

nrc <- function(inclusionProbabilities, selected, sampled=length(inclusionProbabilities)){
  lost <- selected - sampled
  ex <- 1 - (lost / selected)
  return(1 - (1-inclusionProbabilities)**ex)
}

