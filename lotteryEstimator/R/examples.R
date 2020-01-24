
#' Simple estimator for NSSH
#' @description
#'  Estimator for 2-stage lottery sampling, assuming with replacement haul-selection and ignoring intra-haul variance.
#' @details
#'  This estimator is applicable for two stage sampling based on the lottery system.
#'  PSU: haul, assumed UEPWR, Hansen-Hurwitz estimator
#'  SSU: bucket, assumed SRSWOR, Horvitz-Thompson estimator
#' @param samples samples formatted as \code{\link[lotteryEstimator]{NSSH2019}}
#' @param minAge lowest age to produce estimate for
#' @param maxAge highest age to produce estimate for
#' @return estimates, formatted as \code{\link[lotteryEstimator]{catchAtAge}}
#' @examples
#'  data(NSSH2019)
#'  estimates <- simpleNsshEstimator(NSSH2019, 1,20)
#' @export
simpleNsshEstimator <- function(samples, minAge=1, maxAge=20){


  # get inclusion probabilities for SSUs
  NSSU <- stats::aggregate(list(NSSU = samples$SSUid), by=list(PSUid=samples$PSUid), function(x){length(unique(x))})
  samples <- merge(samples, NSSU)
  samples$SSUinclusionProb <- samples$SSUselectionProb * samples$NSSU

  PSUtotals <- list()
  PSUids <- unique(samples$PSUid)
  PSUselectionProbabilities <- c()
  for (i in 1:length(PSUids)){
    PSU <- PSUids[i]
    SSUtotals <- list()
    SSUinclusionProbabilities <- c()
    SSUids <- unique(samples$SSUid[samples$PSUid == PSU])
    for (j in 1:length(SSUids)){
      SSU <- SSUids[j]
      SSUtotals[[j]] <- countCategorical(samples$age[samples$SSUid == SSU], minAge:maxAge)
      SSUinclusionProbabilities <- c(SSUinclusionProbabilities, samples$SSUinclusionProb[match(SSU, samples$SSUid)])
    }
    PSUtotals[[i]] <- horvitzThompson(SSUtotals, SSUinclusionProbabilities)
    PSUselectionProbabilities <- c(PSUselectionProbabilities, samples$PSUselectionProb[match(PSU, samples$PSUid)])
  }

  total <- hansenHurwitz(PSUtotals, PSUselectionProbabilities)
  variance <- hansenHurwitzCovariance(PSUtotals, PSUselectionProbabilities)

  result <- list()
  result$catchAtAge <- total
  result$covariance <- variance

  # standardize age range and estimate
  return(result)
}
