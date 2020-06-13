#' Makes example population from logbook data
#' @noRd
makeExamplePopulation <- function(logbookfile){
  logb <- RstoxData::readErsFile(logbookfile)

  logb$vesselId <- as.character(as.integer(rank(logb$RC)))
  logb$opid <- as.character(as.integer(rank(paste(logb$RC, logb$STARTTIDSPUNKT))))
  logb$timeid <- as.character(as.integer(rank(substr(logb$STARTTIDSPUNKT,1,10))))
  logb$spaceid <- as.character(as.integer(rank(logb$LOKASJON_START)))
  logb$spacetimeid <- paste(logb$timeid, logb$spaceid, sep="-")

  # keep operations that caught cod
  # but keep records of other catch as well
  codops <- logb$opid[!is.na(logb$FANGSTART_FAO) & logb$FANGSTART_FAO == "COD"]
  codset <- logb[!is.na(logb$FANGSTART_FAO) & logb$opid %in% codops,]

  #keep only one catch pr species
  codset$catchid <- paste(codset$opid, codset$FANGSTART_FAO)
  codset <- codset[!duplicated(codset$catchid),]

  #keep only lines
  codset <- codset[codset$REDSKAP_FAO == "LLS" | codset$REDSKAP_FAO == "LL",]


  codset$speciesFAO <- codset$FANGSTART_FAO
  codset$wholeWeightKg <- codset$RUNDVEKT
  codset$gearFAO <- codset$REDSKAP_FAO
  codset$hookCount <- codset$INNSATS
  codset$soaktimeHours <- codset$VARIGHET/60.0

  export <- data.table::data.table(codset[,c("vesselId", "timeid", "spaceid", "spacetimeid", "opid", "catchid", "gearFAO", "hookCount", "soaktimeHours", "speciesFAO", "wholeWeightKg")])
  return(export)
}

#' Sample example population
#' @details Performs two-stage sampling without replacement from an example population
#' @param nPSU the number of primary sampling units to sample
#' @param nSSU the maximal number of secondary sampling untis to sample
#' @param population \code{\link[data.table]{data.table}} with the entire popoulation to sample. The population must be organised in clusters for two stage sampling, identfied by 'PSU' and 'SSU'.
#' @param PSU character() specifying the column in 'population' that identifies the primary sampling units
#' @param SSU character() specifying the column in 'population' that identifies the secondary sampling units
#' @return \code{\link[data.table]{data.table}} formatted as 'population' containing only the sampled elements, and the following additional columns:
#'  \describe{
#'   \item{Npsu}{total number of PSUs in the population}
#'   \item{Nssu}{total number of SSUs within the given PSU in the population}
#'   \item{nPSU}{total number of PSUs selected}
#'   \item{nSSU}{total number of SSUs selected withing the given PSU.}
#'   \item{ipPSU}{inclusion probability of PSU}
#'   \item{spPSU}{selection probability of PSU}
#'   \item{ipSSU}{inclusion probability of SSU, given that PSU is selected}
#'   \item{spSSU}{selection probability of SSU, given that PSU is selected}
#'  }
#' @examples
#'  sample <- sampleExamplePopulation(10,10)
#' @export
sampleExamplePopulation <- function(nPSU, nSSU, population=lotteryEstimator::examplePopulation, PSU="vesselId", SSU="opid"){

  #annotate total PSUs and SSUs
  psuList <- list(population[[PSU]])
  names(psuList) <- PSU
  ssuTotal <- stats::aggregate(list(Nssu=population[[SSU]]), by=psuList, FUN=function(x){length(unique(x))})

  population <- merge(population, ssuTotal)
  population$Npsu <- length(unique(population[[PSU]]))

  #annotate sample sizes
  population$nPSU <- nPSU
  population$nSSU <- pmin(population$Nssu, nSSU)

  #annotate sampling probabilities
  population$ipPSU <- population$nPSU / population$Npsu
  population$ipSSU <- population$nSSU / population$Nssu

  sample <- resample(population, hierarchy=c(PSU, SSU), nSamples=c("nPSU", "nSSU"), replacement=c(F,F))

  PSUselectionOrder <- unlist(lapply(strsplit(sample[[PSU]], "#"), function(x){as.integer(x[2])}))
  SSUselectionOrder <- unlist(lapply(strsplit(sample[[SSU]], "#"), function(x){as.integer(x[3])}))

  sample$spPSU <- 1 / (sample$Npsu - PSUselectionOrder + 1)
  sample$spSSU <- 1 / (sample$Nssu - SSUselectionOrder + 1)

  return(sample)
}
