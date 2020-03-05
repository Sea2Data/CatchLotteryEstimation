
#' Reference implementation. Kept for testing purposes
#' @noRd
simpleNsshEstimatorReference <- function(samples, minAge=1, maxAge=20){

  sampledSSUs <- stats::aggregate(list(sampledSSUs=samples$SSUid), by=list(PSUid=samples$PSUid), FUN=function(x){length(unique(x))})
  samples <- merge(samples, sampledSSUs)

  # calculate inclusion probability for SSUs
  samples$SSUinclusionProb <- samples$SSUselectionProb * samples$sampledSSUs

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

  sampledSSUs <- stats::aggregate(list(sampledSSUs=samples$SSUid), by=list(PSUid=samples$PSUid), FUN=function(x){length(unique(x))})
  samples <- merge(samples, sampledSSUs)

  # calculate inclusion probability for SSUs
  samples$SSUinclusionProb <- samples$SSUselectionProb * samples$sampledSSUs

  # cencus in each bucket
  numAtAgeSample <- function(sample){countCategorical(sample$age, minAge:maxAge)}

  # Horwitz Thompson for haul (estimating from SSUs)
  numAtAgeHaul <- function(sample){hierarchicalHorwitzThompson(sample, "SSUid", numAtAgeSample, "SSUinclusionProb")}

  # Hansen Hurwitz for total (estimating from PSUS)
  numAtAgeTotal <- hierarchicalHansenHurwitz(samples, "PSUid", numAtAgeHaul, "PSUselectionProb")

  # Hansen Hurwitz for covairance, assuming 0 intra-haul covariance
  covariance <- hierarchicalHansenHurwitzCovariance(samples, "PSUid", numAtAgeHaul, function(x){0}, "PSUselectionProb")

  result <- list()
  result$catchAtAge <- numAtAgeTotal
  result$covariance <- covariance

  return(result)
}

#' Estimator for NSSH
#' @description
#'  Estimator for 2-stage lottery sampling, assuming with replacement haul-selection and estimating intra-haul variance from sample covariances.
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
#'  estimates <- twoStageNsshEstimator(NSSH2019, 1,20)
#' @export
twoStageNsshEstimator <- function(samples, minAge=1, maxAge=20){

  # get number of fish in each haul
  samples$nFishHaul <- samples$nFish * samples$nSSU

  # estimate for each haul
  numAtAgeHaul <- function(sample){proportionCategorical(sample$age, minAge:maxAge) * sample$nFishHaul[1]}

  # Hansen Hurwitz for total (estimating from PSUS)
  numAtAgeTotal <- hierarchicalHansenHurwitz(samples, "PSUid", numAtAgeHaul, "PSUselectionProb")

  # sample covariance for for each haul
  intraHaulCovariance <- function(sample){calculateSampleProportionCovariance(proportionCategorical(sample$age, minAge:maxAge)) * (sample$nFishHaul[1]**2 / ((sample$nFish[1])*(sample$nFish[1]-1)))}

  # Hansen Hurwitz for covariance for total
  covariance <- hierarchicalHansenHurwitzCovariance(samples, "PSUid", numAtAgeHaul, intraHaulCovariance, "PSUselectionProb")

  result <- list()
  result$catchAtAge <- numAtAgeTotal
  result$covariance <- covariance
  return(result)
}

#' Estimator for NSSH
#' @description
#'  Estimator for 2-stage length-stratified lottery sampling, assuming with replacement haul-selection, and no intra-haul covariance.
#' @details
#'  This estimator is applicable for two stage sampling based on the lottery system.
#'  PSU: haul, assumed UEPWR, Hansen-Hurwitz estimator
#'  SSU: bucket, assumed SRSWOR
#' @param samples samples formatted as \code{\link[lotteryEstimator]{NSSH2019}}
#' @param minAge lowest age to produce estimate for
#' @param maxAge highest age to produce estimate for
#' @return estimates, formatted as \code{\link[lotteryEstimator]{catchAtAge}}
#' @examples
#'  data(NSSH2019)
#'  estimates <- twoStageNsshEstimator(NSSH2019, 1,20)
#' @export
twoStageNsshEstimatorStratified <- function(samples, minAge=1, maxAge=20){

  lengthStrata <- unique(samples$lengthStrata)

  # get number of fish in each haul
  samples$nFishHaul <- samples$nFish * samples$nSSU

  # estimate for each haul
  numAtAgeHaul <- function(sample){
    m <- sample$nFishHaul[1] * calculateAgeLengthKey(sample$age, sample$lengthStrata, minAge:maxAge, lengthStrata) %*% proportionCategorical(sample$lengthStrata, lengthStrata)
    v <- as.vector(m)
    names(v) <- rownames(m)
    return(v)}

  # Hansen Hurwitz for total (estimating from PSUS)
  numAtAgeTotal <- hierarchicalHansenHurwitz(samples, "PSUid", numAtAgeHaul, "PSUselectionProb")

  # Hansen Hurwitz for covariance for total
  covariance <- hierarchicalHansenHurwitzCovariance(samples, "PSUid", numAtAgeHaul, function(x){0}, "PSUselectionProb")

  result <- list()
  result$catchAtAge <- numAtAgeTotal
  result$covariance <- covariance
  return(result)
}

#' Estimator for NSSH
#' @description
#'  Estimator for 2-stage lottery sampling, assuming with replacement haul-selection and estimating intra-haul variance from bootstrap.
#' @details
#'  As \code{\link[lotteryEstimator]{twoStageNsshEstimator}}, but with intra-haul variance bootstrapped
#'  by selection without replacement.
#' @param samples samples formatted as \code{\link[lotteryEstimator]{NSSH2019}}
#' @param minAge lowest age to produce estimate for
#' @param maxAge highest age to produce estimate for
#' @param iterations number of bootstrap iterations to run for each haul.
#' @return estimates, formatted as \code{\link[lotteryEstimator]{catchAtAge}}
#' @examples
#'  data(NSSH2019)
#'  estimates <- twoStageNsshEstimatorBootstrapped(NSSH2019, 1,20)
#' @export
twoStageNsshEstimatorBootstrapped <- function(samples, minAge=1, maxAge=20, iterations=10){

  # get number of fish in each haul
  samples$nFishHaul <- samples$nFish * samples$nSSU

  # estimate for each haul
  numAtAgeHaul <- function(sample){proportionCategorical(sample$age, minAge:maxAge) * sample$nFishHaul[1]}

  # Hansen Hurwitz for total (estimating from PSUS)
  numAtAgeTotal <- hierarchicalHansenHurwitz(samples, "PSUid", numAtAgeHaul, "PSUselectionProb")

  # bootstrapped WOR covariance for for each haul
  intraHaulCovariance <- function(sample){bootstrap(sample, numAtAgeHaul, iterations, c("FishId"), replacement = c(F), popSize = c("nFishHaul"))$covariances}

  # Hansen Hurwitz for covariance for total
  covariance <- hierarchicalHansenHurwitzCovariance(samples, "PSUid", numAtAgeHaul, intraHaulCovariance, "PSUselectionProb")

  result <- list()
  result$catchAtAge <- numAtAgeTotal
  result$covariance <- covariance
  return(result)
}



defaultLengthStrata <- as.list(c(0,seq(.25,.39,.02),41))
names(defaultLengthStrata) <- paste(unlist(defaultLengthStrata), unlist(defaultLengthStrata)[2:length(defaultLengthStrata)], sep="-")
#' Emulates length stratification by assigning samples to strata and resampling to lower number of fish
#' @param sample, formatted as \code{\link[lotteryEstimator]{NSSH2019}}
#' @param lengthStrata, named list mapping strata names to shortest length that should be included in that strata
#' @param numberPrStrata The number of fish to sample in each strata
#' @noRd
makeLengthStratifiedSubsample <- function(sample, lengthStrata=defaultLengthStrata, numberPrStrata=2){

  namesSample <- names(sample)

  stopifnot(all(!is.na(sample$age)))
  stopifnot(all(!is.na(sample$length)))

  lengthStrata <- lengthStrata[order(unlist(lengthStrata))]
  sample$lengthStrata <- ""
  for (i in 1:length(lengthStrata)){
    sample$lengthStrata[sample$length >= lengthStrata[i]] <- names(lengthStrata)[i]
  }

  nPrStrata <- stats::aggregate(list(N=sample$length), by=list(lengthStrata=sample$lengthStrata, SSUid=sample$SSUid), FUN=length)
  nPrStrata$select <- numberPrStrata
  nPrStrata$select[nPrStrata$N < numberPrStrata] <- nPrStrata$N[nPrStrata$N < numberPrStrata]

  sample <- data.table::merge.data.table(nPrStrata, sample)
  sample$SSUidOriginal <- sample$SSUid
  sample$lengthStrataOriginal <- sample$lengthStrata

  resample <- resample(sample, hierarchy = c("SSUid", "lengthStrata", "FishId"), nSamples = c(NA, NA, "select"), replacement = c(F,F,F), popSize = c(NA, NA, "N"))

  resample$SSUid <- resample$SSUidOriginal
  resample$lengthStrata <- resample$lengthStrataOriginal
  resample$lengthStrataTotal <- resample$N
  resample <- resample[,c(names(resample)[names(resample) %in% namesSample], "lengthStrata", "lengthStrataTotal")]
  return(resample)

}
