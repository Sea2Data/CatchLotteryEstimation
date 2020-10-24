
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
  if (any(selectionProbabilities > 1) | any(selectionProbabilities <= 0)){
    stop("all selectionProbabilities must be in [0,1>")
  }
  if (length(sampleTotals) < 2){
    stop("estimator not defined for less than two samples.")
  }

  sampNames <- names(sampleTotals)
  selNames <- names(selectionProbabilities)

  if (!is.null(sampNames) & !is.null(selNames)){
    if (!(all(sampNames == selNames))){
      stop("Vector names does not match between sampleTotals and selectionProbabilities")
    }
  }

  point <- hansenHurwitz(sampleTotals, selectionProbabilities)

  sumSamples <- outer(sampleTotals[[1]], sampleTotals[[1]])*0

  for (i in 1:length(sampleTotals)){
      dev <- (sampleTotals[[i]] / selectionProbabilities[i]) - point
      sumSamples <- sumSamples + outer(dev, dev)
  }


  n <- length(sampleTotals)
  sumSamples <- sumSamples / (n*(n-1))

  #make matrix if vector
  if (is.vector(sumSamples)){
    nn <- list(names(sumSamples), names(sumSamples))
    dim(sumSamples) <- c(1,1)
    dimnames(sumSamples) <- nn
  }

  return(sumSamples)

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
  if (any(selectionProbabilities > 1) | any(selectionProbabilities <= 0)){
    stop("all selectionProbabilities must be in [0,1>")
  }

  sampNames <- names(sampleCovariances)
  selNames <- names(selectionProbabilities)

  if (!is.null(sampNames) & !is.null(selNames)){
    if (!(all(sampNames == selNames))){
      stop("Vector names does not match between sampleTotals and selectionProbabilities")
    }
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
    stop(paste("inclusionProbabilities does not correspond to the listed sample totals (", length(inclusionProbabilities),",", length(sampleTotals), ")"))
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
#' @param coInclusionProbabilities matrix() matrix of co-inclusion probabilites (joint-inclusionprobabilities) for the samples whoose covariance are listed in 'sampleCovariances'.
#' @return matrix() representing a symmetric matrix with the estimated covariances.
#' @export
horvitzThompsonCovariance <- function(sampleTotals, inclusionProbabilities, coInclusionProbabilities){
  if (length(sampleTotals) != length(inclusionProbabilities)){
    stop(paste("inclusionProbabilities does not correspond to the listed sample covariances", length(sampleTotals), "vs", length(inclusionProbabilities)))
  }
  if (nrow(coInclusionProbabilities) != length(inclusionProbabilities)){
    stop(paste("inclusionProbabilities does not correspond to the listed coInclusionProbabilities", length(inclusionProbabilities), "vs", nrow(coInclusionProbabilities), "x", ncol(coInclusionProbabilities)))
  }
  if (any(inclusionProbabilities > 1) | any(inclusionProbabilities <= 0)){
    stop("all inclusionProbabilities must be in [0,1>")
  }

  sampNames <- names(sampleTotals)
  inclNames <- names(inclusionProbabilities)
  rowNames <- rownames(coInclusionProbabilities)
  colNames <- colnames(coInclusionProbabilities)

  if (!is.null(sampNames) & !is.null(inclNames)){
    if (!(all(sampNames == inclNames))){
      stop("Vector names does not match between sampleTotals and inclusionProbabilities")
    }
  }
  if (!is.null(sampNames) & !is.null(rowNames)){
    if (!(all(sampNames == rowNames))){
      stop("Vector names does not match between sampleTotals and coInclusionProbabilities")
    }
  }
  if (!is.null(sampNames) & !is.null(colNames)){
    if (!(all(sampNames == colNames))){
      stop("Vector names does not match between sampleTotals and coInclusionProbabilities")
    }
  }

  # cache outer products in 4-d array
  ops <- array(dim=c(length(sampleTotals), length(sampleTotals), dim(outer(sampleTotals[[1]], sampleTotals[[1]]))))
  for (i in 1:length(sampleTotals)){
    for (j in i:length(sampleTotals)){
     ops[i,j,,] <- outer(sampleTotals[[i]], sampleTotals[[j]])
     ops[j,i,,] <- ops[i,j,,]
    }
  }

  sumSamples <- (1-inclusionProbabilities[1]) * ops[1,1,,] / inclusionProbabilities[1]**2

  if (length(sampleTotals) > 1){
    for (i in 2:length(sampleTotals)){
      sumSamples <- sumSamples + (1-inclusionProbabilities[i]) * ops[i,i,,] / inclusionProbabilities[i]**2
    }
  }

  for (i in 1:length(sampleTotals)){
    for (j in i:length(sampleTotals)){
      if (i != j){
        sumSamples <- sumSamples + 2*(coInclusionProbabilities[i,j]-inclusionProbabilities[i]*inclusionProbabilities[j]) * ops[i,j,,] / (coInclusionProbabilities[i,j] * inclusionProbabilities[i] * inclusionProbabilities[j])
      }
    }
  }

  #make matrix if vector
  if (is.vector(sumSamples)){
    nn <- list(names(sumSamples), names(sumSamples))
    dim(sumSamples) <- c(1,1)
    dimnames(sumSamples) <- nn
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
#' @export
horvitzThompsonIntra <- function(sampleCovariances, inclusionProbabilities){
  if (length(sampleCovariances) != length(inclusionProbabilities)){
    stop("inclusionProbabilities does not correspond to the listed sample covariances")
  }
  if (any(inclusionProbabilities > 1) | any(inclusionProbabilities <= 0)){
    stop("all inclusionProbabilities must be in [0,1>")
  }

  sampNames <- names(sampleCovariances)
  inclNames <- names(inclusionProbabilities)

  if (!is.null(sampNames) & !is.null(inclNames)){
    if (!(all(sampNames == inclNames))){
      stop("Vector names does not match between sampleTotals and inclusionProbabilities")
    }
  }

  sumSamples <- sampleCovariances[[1]] / inclusionProbabilities[1]

  if (length(sampleCovariances) > 1){
    for (i in 2:length(sampleCovariances)){
      sumSamples <- sumSamples + (sampleCovariances[[i]] / inclusionProbabilities[i])
    }
  }

  return(sumSamples)

}


#
# Hiearchical estimators
#


#' Hierarchical stratified estimate
#' @description
#'  \code{\link[lotteryEstimator]{HierarchicalEstimator}} for estimating totals for a stratified sample
#' @param sample \code{\link[data.table]{data.table}} with sample data
#' @param partitionId character() identifying the columns in 'sample' that identify the strata to estimate from
#' @param subEstimator \code{\link[lotteryEstimator]{ParameterizedEstimator}} for estimating totals for each strata
#' @return vector of parameter estimates
#' @export
hierarchicalStratifiedTotals <- function(sample, partitionId, subEstimator){
  if (!is.character(sample[[partitionId]])){
    stop("partitionId must be a character")
  }

  strataTotals <- list()
  for (id in unique(sample[[partitionId]])){
    sampleUnitData <- sample[sample[[partitionId]] == id,]
    strataTotals[[id]] <- subEstimator(sampleUnitData)
  }
  return(Reduce('+', strataTotals))
}

#' Hierarchical stratified estimate
#' @description
#'  \code{\link[lotteryEstimator]{HierarchicalCovarianceEstimator}} for estimating covariances for a stratified sample
#' @param sample \code{\link[data.table]{data.table}} with sample data
#' @param partitionId character() identifying the columns in 'sample' that identify the strata to estimate from
#' @param subEstimator \code{\link[lotteryEstimator]{ParameterizedCovarianceEstimator}} for estimating covariances for each strata
#' @return matrix of covariance estimates
#' @export
hierarchicalStratifiedCovariance <- function(sample, partitionId, subEstimator){

  #
  # implementation identical to hierarchicalStratifiedTotals,
  # but keep distinct for clarity and adding of format checks
  #

  if (!is.character(sample[[partitionId]])){
    stop("partitionId must be a character")
  }

  strataCovariance <- list()
  for (id in unique(sample[[partitionId]])){
    sampleUnitData <- sample[sample[[partitionId]] == id,]
    strataCovariance[[id]] <- subEstimator(sampleUnitData)
  }
  return(Reduce('+', strataCovariance))

}

#' Hierarchical Hansen-Hurwitz
#' @description
#'  \code{\link[lotteryEstimator]{HierarchicalEstimator}} for applying
#'  \code{\link[lotteryEstimator]{hansenHurwitz}} in hierarchical implementations.
#' @param sample \code{\link[data.table]{data.table}} with sample data
#' @param partitionId character() identifying the columns in 'sample' that identify the sampling units to estimate from
#' @param subEstimator \code{\link[lotteryEstimator]{ParameterizedEstimator}} for estimating totals for each sampled unit
#' @param selectionProbabilities character() indentifying the columns in 'sample' with selection probabilites for the sampling units.
#' @return numeric() estimate of total
#' @examples
#'  numAtAgeSample <- function(sample){countCategorical(sample$age, 2:20)}
#'  numAtAgeHaul <- function(sample){hierarchicalHorvitzThompsonTotals(sample, "SSUid",
#'                                   numAtAgeSample, "SSUinclusionProb")}
#'  exampleSamples <- lotteryEstimator::NSSH2019
#'  exampleSamples$SSUinclusionProb <- exampleSamples$SSUselectionProb * exampleSamples$nSSU
#'  numAtAgeTotal <- hierarchicalHansenHurwitzTotals(exampleSamples, "PSUid",
#'                                             numAtAgeHaul, "PSUselectionProb")
#' @export
hierarchicalHansenHurwitzTotals <- function(sample, partitionId, subEstimator, selectionProbabilities){

  if (!is.character(sample[[partitionId]])){
    stop("partitionId must be a character")
  }

  sampleTotals <- list()
  sProb <- list()
  for (id in unique(sample[[partitionId]])){
    sampleUnitData <- sample[sample[[partitionId]] == id,]
    sProb[[id]] <- sampleUnitData[[selectionProbabilities]][1]
    sampleTotals[[id]] <- subEstimator(sampleUnitData)
  }
  return(hansenHurwitz(sampleTotals, unlist(sProb)))

}

#' Hierarchical Hansen-Hurwitz domain estimate
#' @description
#'  \code{\link[lotteryEstimator]{HierarchicalEstimator}} for estimating domain totals using Hansen-Hurwitz estimators.
#'  The estimator is derived from a ratio estimator and is not unbiased.
#'  The bias is low if the number of sampled units in the domain is correlated with the sample totals for the domain.
#' @param sample \code{\link[data.table]{data.table}} with sample data
#' @param domainIds character() identifying the column in 'sample' that identify a partitioning of the sample into domains
#' @param subEstimator \code{\link[lotteryEstimator]{ParameterizedEstimator}} for estimating totals for each sampled unit
#' @param domain character() name of the domain to estimate for (must occur in the column 'domainIds')
#' @param domainFraction numeric() the fraction of the population that fall in the domain.
#' @examples
#'  numAtAgeSample <- function(sample){countCategorical(sample$age, 2:20)}
#'  numAtAgeHaul <- function(sample){hierarchicalHorvitzThompsonTotals(sample, "SSUid",
#'                                   numAtAgeSample, "SSUinclusionProb")}
#'  exampleSamples <- lotteryEstimator::NSSH2019
#'
#'  exampleSamples$SSUinclusionProb <- exampleSamples$SSUselectionProb * exampleSamples$nSSU
#'  numAtAgeTotal <- function(sample){hierarchicalHansenHurwitzTotals(sample, "PSUid",
#'                                             numAtAgeHaul, "PSUselectionProb")}
#'
#'  #define som arbitrary domains for the sake of the example
#'  exampleSamples$domains <- cut(as.numeric(exampleSamples$PSUid), 3, labels = c("one", "two", "three"))
#'
#'  #obtain domain estimate
#'  hierarchicalHansenHurwitzDomainTotals(exampleSamples, "domains", numAtAgeTotal, "two", 1/3)
#' @export
hierarchicalHansenHurwitzDomainTotals <- function(sample, domainIds, subEstimator, domain, domainFraction){

  if (!(domain %in% sample[[domainIds]])){
    stop(paste("No data for domain", domain))
  }

  sample <- sample[sample[[domainIds]]==domain,]
  domainTotal <- domainFraction*subEstimator(sample)

  return(domainTotal)
}

#' Hierarchical Hansen-Hurwitz domain covariance
#' @description
#'  \code{\link[lotteryEstimator]{HierarchicalCovarianceEstimator}} for estimating the covariance of domain totals using Hansen-Hurwitz estimators
#' @details
#'  The estimator is derived from a ratio estimator and is not unbiased.
#'  In particular, the variation in the selection of domain members is not taken into account.
#'  The bias is low if the number of sampled units in the domain is correlated with the sample totals for the domain.
#'  The 'domainFraction' is assumed known. If this is obtained by estimation from sample,
#'  the variance contribution from that estimate is not taken into account by this function.
#' @param sample \code{\link[data.table]{data.table}} with sample data
#' @param domainIds character() identifying the column in 'sample' that identify a partitioning of the sample into domains
#' @param subCovarianceEstimator \code{\link[lotteryEstimator]{ParameterizedCovarianceEstimator}} for estimating covariances for each sampled unit.
#' @param domain character() name of the domain to estimate for (must occur in the column 'domainIds')
#' @param domainFraction numeric() the fraction of the population that fall in the domain.
#' @examples
#'  numAtAgeSample <- function(sample){countCategorical(sample$age, 2:20)}
#'  numAtAgeHaul <- function(sample){hierarchicalHorvitzThompsonTotals(sample, "SSUid",
#'                                   numAtAgeSample, "SSUinclusionProb")}
#'  exampleSamples <- lotteryEstimator::NSSH2019
#'  exampleSamples$SSUinclusionProb <- exampleSamples$SSUselectionProb * exampleSamples$nSSU
#'  hhCov <- function(samples){hierarchicalHansenHurwitzCovariance(samples, "PSUid",
#'                                                    numAtAgeHaul, function(x){0},
#'                                                    "PSUselectionProb")}
#'
#'  #define som arbitrary domains for the sake of the example
#'  exampleSamples$domains <- cut(as.numeric(exampleSamples$PSUid), 3, labels = c("one", "two", "three"))
#'
#'  #obtain domain estimate
#'  hierarchicalHansenHurwitzDomainTotalCovariance(exampleSamples, "domains", hhCov, "two", 1/3)
#' @export
hierarchicalHansenHurwitzDomainTotalCovariance <- function(sample, domainIds, subCovarianceEstimator, domain, domainFraction){

  if (!(domain %in% sample[[domainIds]])){
    stop(paste("No data for domain", domain))
  }

  sample <- sample[sample[[domainIds]]==domain,]
  domainCov <- (domainFraction**2)*subCovarianceEstimator(sample)

  return(domainCov)
}

#' Hierarchical Hansen-Hurwitz covariance
#' @description
#'  \code{\link[lotteryEstimator]{HierarchicalCovarianceEstimator}} for estimating
#'  covariance of \code{\link[lotteryEstimator]{hansenHurwitz}} in hierarchical implementations.
#' @param sample \code{\link[data.table]{data.table}} with sample data
#' @param partitionId character() identifying the columns in 'sample' that identify the sampling units to estimate from
#' @param subEstimator function, function for estimating totals for each sampled unit
#' @param subCovarianceEstimator \code{\link[lotteryEstimator]{ParameterizedCovarianceEstimator}} for estimating covariances for each sampled unit. If NULL the sampled units will be assumed to have no contribution to the covariance
#' @param selectionProbabilities character() indentifying the columns in 'sample' with selection probabilites for the sampling units.
#' @return numeric() estimate of total
#' @examples
#'  numAtAgeSample <- function(sample){countCategorical(sample$age, 2:20)}
#'  numAtAgeHaul <- function(sample){hierarchicalHorvitzThompsonTotals(sample, "SSUid",
#'                                   numAtAgeSample, "SSUinclusionProb")}
#'  exampleSamples <- lotteryEstimator::NSSH2019
#'  exampleSamples$SSUinclusionProb <- exampleSamples$SSUselectionProb * exampleSamples$nSSU
#'  covariance <- hierarchicalHansenHurwitzCovariance(exampleSamples, "PSUid",
#'                                                    numAtAgeHaul, function(x){0},
#'                                                    "PSUselectionProb")
#' @export
hierarchicalHansenHurwitzCovariance <- function(sample, partitionId, subEstimator, subCovarianceEstimator, selectionProbabilities){
  sampleTotals <- list()
  sampleCovariances <- list()
  sProb <- list()
  for (id in unique(sample[[partitionId]])){
    sampleUnitData <- sample[sample[[partitionId]] == id,]
    sProb[[id]] <- sampleUnitData[[selectionProbabilities]][1]
    sampleTotals[[id]] <- subEstimator(sampleUnitData)
    if (!is.null(subCovarianceEstimator)){
      sampleCovariances[[id]] <- subCovarianceEstimator(sampleUnitData)
    }
  }
  cov <- hansenHurwitzCovariance(sampleTotals, unlist(sProb))

  if (!is.null(subCovarianceEstimator)){
    intra <- hansenHurwitzIntra(sampleCovariances, unlist(sProb))
    cov <- cov + intra
  }

  return(cov)

}

#' Hierarchical Horvitz-Thompson domain estimate
#' @description
#'  \code{\link[lotteryEstimator]{HierarchicalEstimator}} for estimating domain totals using Horvitz-Thompson estimators.
#' @details
#'  Note that no corresponding function for Horvitz-Thompson domain variance estimation is defined.
#' @param sample \code{\link[data.table]{data.table}} with sample data
#' @param domainIds character() identifying the column in 'sample' that identify a partitioning of the sample into domains
#' @param subEstimator \code{\link[lotteryEstimator]{ParameterizedEstimator}} for estimating totals for each sampled unit
#' @param domain character() name of the domain to estimate for (must occur in the column 'domainIds')
#' @export
hierarchicalHorvitzThompsonDomainTotals <- function(sample, domainIds, subEstimator, domain){

  if (!(domain %in% sample[[domainIds]])){
    stop(paste("No data for domain", domain))
  }

  sample <- sample[sample[[domainIds]]==domain,]
  domainTotal <- subEstimator(sample)

  return(domainTotal)
}

#' Hierarchical Horvitz Thompson
#' @description
#'  \code{\link[lotteryEstimator]{HierarchicalEstimator}} for applying
#'  \code{\link[lotteryEstimator]{horvitzThompson}} in hierarchical implementations.
#' @details
#'  Example provided in the documentation for \code{\link[lotteryEstimator]{hierarchicalHansenHurwitzTotals}}
#' @param sample \code{\link[data.table]{data.table}} with sample data
#' @param partitionId character() identifying the columns in 'sample' that identify the sampling units to estimate from
#' @param subEstimator \code{\link[lotteryEstimator]{ParameterizedEstimator}} for estimating totals for each sample
#' @param inclusionProbabilities character() indentifying the columns in 'sample' with inclusion probabilites for the sampling units.
#' @return numeric() estimate of total
#' @export
hierarchicalHorvitzThompsonTotals <- function(sample, partitionId, subEstimator, inclusionProbabilities){

  if (!is.character(sample[[partitionId]])){
    stop("partitionId must be a character")
  }


  sampleTotals <- list()
  iProb <- list()

  data.table::setkeyv(sample, partitionId)

  for (id in unique(sample[[partitionId]])){
    sampleUnitData <- sample[list(id)]
    iProb[[id]] <- sampleUnitData[[inclusionProbabilities]][1]
    sampleTotals[[id]] <- subEstimator(sampleUnitData)
  }
  return(horvitzThompson(sampleTotals, unlist(iProb)))

}


#' Hierarchical Horvitz-Thompson covariance
#' @description
#'  \code{\link[lotteryEstimator]{HierarchicalCovarianceEstimator}} for estimating
#'  covariance of \code{\link[lotteryEstimator]{horvitzThompson}} in hierarchical implementations.
#' @param sample \code{\link[data.table]{data.table}} with sample data
#' @param partitionId character() identifying the columns in 'sample' that identify the sampling units to estimate from
#' @param subEstimator function, function for estimating totals for each sampled unit
#' @param subCovarianceEstimator \code{\link[lotteryEstimator]{ParameterizedCovarianceEstimator}} for estimating covariances for each sampled unit. If NULL the sampled units will be assumed to have no contribution to the covariance
#' @param inclusionProbabilities character() indentifying the columns in 'sample' with inclusion probabilites for the sampling units.
#' @param coInclusionMatrix matrix() named matrix with co-inclusion probabilites, indexed with the values in the column identified by 'partitionId'.
#' @return numeric() estimate of total
#' @export
hierarchicalHorvitzThompsonCovariance <- function(sample, partitionId, subEstimator, subCovarianceEstimator, inclusionProbabilities, coInclusionMatrix){

  sampleTotals <- list()
  sampleCovariances <- list()
  iProb <- list()
  for (id in unique(sample[[partitionId]])){
    sampleUnitData <- sample[sample[[partitionId]] == id,]
    iProb[[id]] <- sampleUnitData[[inclusionProbabilities]][1]
    sampleTotals[[id]] <- subEstimator(sampleUnitData)
    if (!is.null(subCovarianceEstimator)){
      sampleCovariances[[id]] <- subCovarianceEstimator(sampleUnitData)
    }
  }

  cov <- horvitzThompsonCovariance(sampleTotals, unlist(iProb), coInclusionMatrix)

  if (!is.null(subCovarianceEstimator)){
    intra <- horvitzThompsonIntra(sampleCovariances, unlist(iProb))

    if (length(intra) != length(cov)){
      stop("Inter and Intra covariance have different length.")
    }
    cov <- cov + intra
  }

  return(cov)

}

#
# Some datatype conventions
#

#' Parameterized estimator
#'
#' Function contract for parameterized estimators
#' Parameterized estimators are functions that map a data frame to an estimate, with all other parameters fixed.
#'
#' Parameterized estimators take a single argument 'sample' which is a \code{\link[data.table]{data.table}} with sample data
#'
#' Parameterized estimators return a named numeric vector with parameter estimates
#'
#' @examples
#'  # A parameterized estimator can be obtained from
#'  # \code{\link[lotteryEstimator]{HierarchicalEstimator}}
#'  # by fixing parameters.
#'  # For example a parametierized Horwitz-Thomposon estimator:
#'  numAtAgeSample <- function(sample){countCategorical(sample$age, 2:20)}
#'  numAtAgeHaul <- function(sample){hierarchicalHorvitzThompsonTotals(sample, "SSUid",
#'                                   numAtAgeSample, "SSUinclusionProb")}
#'
#'
#' @name ParameterizedEstimator
#'
NULL

#' Parameterized covariance estimator
#'
#' Function contract for parameterized covariance estimators
#' Parameterized covariance estimators are functions that map a data table to a covariance estimate, with all other parameters fixed.
#'
#' Parameterized covariance estimators take a single argument 'sample' which is a \code{\link[data.table]{data.table}} with sample data
#'
#' Parameterized covariance estimators return a named matrix with estimates of covariances between parameters
#'
#' Parameterized covariance estimators can be obtained from \code{\link[lotteryEstimator]{HierarchicalCovarianceEstimator}}
#' by fixing parameters, in the same manner that this is done for \code{\link[lotteryEstimator]{ParameterizedEstimator}}
#'
#' @name ParameterizedCovarianceEstimator
#'
NULL


#' Hierarchical estimator
#'
#' Function contract for hierarchical estimators.
#' Hierarchical estimators are estimators that depend on other estimates being made for partitions of the sample.
#' Such partitions can be lower level sampling units, clusters or strata.
#' A hierachical estimator will divide the sample into subsets with data for each partition
#' and call another function to calculate estimates from each of these subsets.
#'
#' These functions take three arguments:
#' \describe{
#'  \item{sample}{\code{\link{data.table}{data.table}} with the sample to estimate from}
#'  \item{...}{other arguments needed to calculate estimate}
#' }
#'
#' A hierarchical estimators return a vector of estimated parameters.
#'
#' @name HierarchicalEstimator
#'
NULL

#' Hierarchical covariance estimator
#'
#' Function contract for hierarchical covariance estimators.
#' Hierarchical covariance estimators are estimators that depend on other estimates being made for partitions of the sample.
#' A hierachical covariance estimator will divide the sample into subsets with data for each partition
#' and call another function to calculate estimates from each of these subsets.
#'
#' These functions take three arguments:
#' \describe{
#'  \item{sample}{\code{\link{data.table}{data.table}} with the sample to estimate from}
#'  \item{...}{other arguments needed to calculate covariance estimates}
#' }
#'
#' A Hierarchical covariance estimator return a covariance matrix with covariances between estimated parameters
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
