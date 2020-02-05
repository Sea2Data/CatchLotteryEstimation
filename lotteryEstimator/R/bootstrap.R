#' Sample indecies
#' @description
#'  samples integers for indexing data tables in resampling.
#' @details
#'  Sampling without replacement is mimicked by sampling from a pseudopopulation.
#'  Effectively the integers up to an including 'maxindex' are replicated
#'  as many times as necessary to reach the size of the population.
#'  Any remaineder is filled in with number selected without replacement from the range 1:'maxindex'.
#'  This psudopopulation is the sampled wihtout replacement for the final result.
#'
#' @param n integer() the number of indecies to select
#' @param maxindex integer() defines the range of indecies to sample (1:'maxindex')
#' @param replacement logical() whether to sample indecies with replacement
#' @param popsize integer() size of the pseudopopulation to be sampled. May be NA if replacement is TRUE.
#' @return integer() vector with sampled integers in the range 1:'maxindex'
#' @export
sampleIndecies <- function(n, maxindex, replacement=T, popsize=NA){
  if (!replacement & is.na(popsize)){
    stop("Population size must be specified for selection without replacement.")
  }
  if (!is.na(popsize) & popsize < maxindex){
    stop("Population size must be larger than or equal to 'maxindex'.")
  }
  if (replacement){
    return(sample.int(maxindex, size=n, replace = T))
  }
  else{

    # highest integer covered by replicating sample an integral number of times
    div <- maxindex * (popsize %/% maxindex)

    #sample population
    popIndecies <- sample.int(popsize, size=n, replace = F)

    #psudopopulation correspondance of fraction below div
    popIndeciesInDiv <- popIndecies[popIndecies<=div]
    indecies <- maxindex - (popIndeciesInDiv %% maxindex)

    #randomize psudopopulation correspondance for any selection from fraction above div
    popIndeciesInRemainder <- sum(popIndecies > div)
    indecies <- c(indecies, sample.int(maxindex, size=popIndeciesInRemainder, replace = F))
    return(indecies)
  }
}


#' resample data
#' @description
#'  Resamples data from hierarchically clustered data,
#'  such as those ariseing from multi-stage sampling with equal selection probabilities at each stage.
#' @details
#'  The parameter 'samples' should have one row for each set of measurements / observations,
#'  and have have columns identifying the sampling units that observations were sampled from.
#'  The parameter 'hierarchy' enumaretes the columns that identify the sampling units in hierarchical order,
#'  primary sampling unit first.
#'
#'  Equal probability selection is assumed for all stages of sampling.
#' @param samples \code{\link[data.table]{data.table}} with samples
#' @param hierarchy character() vector describing the sampling units in hierarchical order
#' @param nSamples character() vector corresponding to 'hierarchy', indicating columns for that specify how many samples should be included in the resampling at the corresponding stage. NAs signify that the number of samples that is available in 'samples' should be included.
#' @param replacement logical() vector corresponding to 'hierarchy', indicating whether the corresponding sampling units should be resampled with replacement
#' @param popSize character() vector corresponding to 'hierarchy', indicating columns for population sizes, NAs signify that the number of samples that is available in 'samples' should be considered to be the population.
#' @param prefix character() a string to prefix used to rename sampling units in resampled data
#' @return resampled data, with sampling units renamed
#' @examples
#'  # The data set contains SSUs where no replicates where sampled
#'  data(NSSH2019)
#'
#'  # set the number of replicates desired for SSU
#'  NSSH2019$SSUrepl <- 2
#'
#'  # Mimicking replicate SSU sampling as selection without replacement at SSU level and
#'  # selection with replacement as if the fish sampled in an SSU is resampled from the enitre haul
#'  rs <- resample(NSSH2019, c("PSUid", "SSUid","FishId"),
#'                           replacement = c(FALSE,FALSE,TRUE),
#'                           nSamples = c(NA,"SSUrepl",NA),
#'                           popSize = c(NA,"nSSU",NA))
#'
#'  # resampling only the SSUs,
#'  # keeping the PSUs and fish as is, leads to exact duplication of samples
#'  # since there is only one SSU pr PSU in this example
#'  rs <- resample(NSSH2019, c("PSUid", "SSUid"),
#'                           replacement = c(FALSE,TRUE),
#'                           nSamples = c(NA,"SSUrepl"))
#'  nrow(rs) == 2*nrow(NSSH2019)
#'  sum(rs$age) == 2*sum(NSSH2019$age)
#'
#'  # replicate sampling at SSU level can be mimicked
#'  # by resampling within SSU.
#'  # Note that the number of rows are still exactly twice the original samples
#'  # because the SSUs are resampled within their PSU.
#'  rs <- resample(NSSH2019, c("PSUid", "SSUid","FishId"),
#'                           replacement = c(FALSE,FALSE,TRUE),
#'                           nSamples = c(NA,"SSUrepl",NA), popSize = c(NA,"nSSU",NA))
#'  nrow(rs) == 2*nrow(NSSH2019)
#'  sum(rs$age) != 2*sum(NSSH2019$age)
#'
#'  # resampling without replacement at all levels,
#'  # and population sizes equal to sample sizes (default)
#'  # simply reorders the data
#'  rs <- resample(NSSH2019, c("PSUid", "SSUid", "FishId"),
#'                           replacement = c(FALSE,FALSE,FALSE))
#'  nrow(rs) == nrow(NSSH2019)
#'  sum(rs$age) == sum(NSSH2019$age)
#'
#' @export
resample <- function(samples, hierarchy, nSamples=rep(NA, length(hierarchy)), replacement=rep(T, length(hierarchy)), popSize=rep(NA, length(hierarchy)), prefix=""){
  if (!all(hierarchy %in% names(samples))){
    missing <- hierarchy[!(hierarchy %in% names(samples))]
    stop("Not all stages identified in hierarchy are present in sample. Missing: ", paste(missing, collapse=","))
  }
  if (length(hierarchy) != length(replacement)){
    stop("hierarchy and replacement must be of same length.")
  }
  if (length(hierarchy) != length(nSamples)){
    stop("hierarchy and nSamples must be of same length.")
  }
  if (length(hierarchy) != length(popSize)){
    stop("hierarchy and popSize must be of same length.")
  }
  for (h in hierarchy){
    if (!is.character(samples[[h]])){
      stop(paste("Columns identifying sampling units must be character(),",h,"is not."))
    }
  }

  if (length(hierarchy) == 0){
    stop("'hierarchy' must be specified")
  }

  currentStage <- hierarchy[1]
  repl <- replacement[1]
  nSamp <- nSamples[1]
  pop <- popSize[1]
  if (length(hierarchy)>1){
    restStages <- hierarchy[2:length(hierarchy)]
    restRepl <- replacement[2:length(replacement)]
    restN <- nSamples[2:length(replacement)]
    restPop <- popSize[2:length(popSize)]
  }
  else{
    restStages <- NULL
    restRepl <- NULL
    restN <- NULL
    restPop <- NULL
  }

  result <- NULL

  #resample this stage
  availableUnits <- unique(samples[[currentStage]])

  N <- length(availableUnits)
  if (!is.na(pop)){
    N <- samples[[pop]][[1]]
  }

  n <- length(availableUnits)
  if (!is.na(nSamp)){
    n <- samples[[nSamp]][[1]]
  }

  selectedUnits <- availableUnits[sampleIndecies(n, length(availableUnits), repl, N)]

  #this stage is the ultimate sampling stage
  if (is.null(restStages)){

    indecies <- c()
    newnames <- c()
    for (i in 1:length(selectedUnits)){
      u <- selectedUnits[i]
      uSamples <- (1:nrow(samples))[samples[[currentStage]]==u]

      newnames <- c(newnames, rep(paste(prefix,currentStage,":",u,"#",i,sep=""), length(uSamples)))
      indecies <- c(indecies, uSamples)
    }
    result <- samples[indecies,]
    result[[currentStage]] <- newnames
  }
  # resample next stage
  else{
    for (i in 1:length(selectedUnits)){
      u <- selectedUnits[i]
      uSamples <- samples[samples[[currentStage]]==u,]

      newname <- paste(prefix,currentStage,":",u,"#",i,sep="")
      newprefix <- paste(newname,"/", sep="")
      rs <- resample(uSamples[,names(uSamples)[names(uSamples) != currentStage]], restStages, restN, restRepl, restPop, newprefix)
      rs[[currentStage]] <- newname
      result <- rbind(rs, result)

    }
  }

  return(result)
}

#' Bootstrap estimate
#' @description
#'  Computes bootstrap estimate of estimate and covariances.
#'  Bootstrap is based on hierarchical resampling of data.
#'
#'  For non-hierarchical sampling, provide a 'hierarchy' as a vector of length one,
#'  containing a column name that identifies sampled units.
#'
#'  The argument 'popSize' is only of importance for sampling without replacement.
#' @param sample sampe \code{\link[data.table]{data.table}} with sample to bootstrap
#' @param estimator function(), must be a simple function of sample, that is take a single argument, formatted as 'sample', and return a numeric vector.
#' @param iteration integer() the number of bootstrap iterations to run
#' @param hierarchy character() vector describing the sampling units in hierarchical order.
#' @param nSamples character() vector corresponding to 'hierarchy', indicating columns for that specify how many samples should be included in the resampling at the corresponding stage. NAs signify that the number of samples that is available in 'samples' should be included.
#' @param replacement logical() vector corresponding to 'hierarchy', indicating whether the corresponding sampling units should be resampled with replacement
#' @param popSize character() vector corresponding to 'hierarchy', indicating columns for population sizes, NAs signify that the number of samples that is available in 'samples' should be considered to be the population.
#' @return list with three members:
#'  \describe{
#'   \item{iterations}{integer() number of iterations run}
#'   \item{meanEstimate}{numeric() vector with mean estimate over bootstrap iterations.}
#'   \item{covariances}{matrix() with covariance of variables in 'meanEstimate' over bootstrap iterations.}
#'  }
#' @examples
#'  #bootstrap fish in haul, leaving haul-selection intact
#'  data(NSSH2019)
#'  ssu1sample <- NSSH2019[NSSH2019$SSUid == NSSH2019$SSUid[1],]
#'  estimator <- function(sample){countCategorical(sample$age, 1:20)}
#'  bs<-bootstrap(NSSH2019, estimator, 1000,
#'                hierarchy = c("FishId"), replacement=c(T))
#' @export
bootstrap <- function(sample, estimator, iterations, hierarchy, nSamples=rep(NA, length(hierarchy)), replacement=rep(T, length(hierarchy)), popSize=rep(NA, length(hierarchy))){

  estimates <- list()
  for (i in 1:iterations){
    bootSample <- resample(sample, hierarchy, nSamples, replacement, popSize)
    estimates[[i]] <- estimator(bootSample)
  }

  meanEst <- Reduce("+", estimates) / length(estimates)
  estDiff <- lapply(estimates, function(x){x - meanEst})
  estDiffOp <- lapply(estDiff, function(x){outer(x,x)})
  covariances <- Reduce("+", estDiffOp) / (length(estimates) - 1)

  results <- list()
  results$iterations <- iterations
  results$meanEstimate <- meanEst
  results$covariances <- covariances

  return(results)

}
