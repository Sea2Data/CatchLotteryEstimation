library(RstoxData)
library(lotteryEstimator)

#
# Functions for preparing data for package lotteryEstimator
#

#' prepares an example data set
#' @param datafile biotic 3 annual data file
#' @param species species code to extract
#' @param totalLanded total landed weight of species/stock (Kg), used to calculate normalized selection probabilites
#' @return data.frame() with one row per measured fish
#' @noRd
extract_set <- function(datafile="~/bioticsets/v3/biotic_cruiseNumber_19-2019-1_Kommersielle+fartøy_2020-01-27T23.01.03.301Z.xml", species="161722.G03", totalLanded=430506000, removemissingCatchWeight=T, removeMissingAges=T){
  biotic <- RstoxData::readXmlFile(datafile)

  catches <- biotic$catchsample[biotic$catchsample$catchcategory == species & biotic$catchsample$missiontype == 19,]
  individuals <- biotic$individual[biotic$individual$missiontype == 19,]
  ages <- biotic$agedetermination[biotic$agedetermination$missiontype == 19,]

  if (removemissingCatchWeight){
    catches <- catches[!is.na(catches$catchweight),]
  }

  if (removeMissingAges){
    ages <- ages[!is.na(ages$age),]
  }

  if (nrow(catches) != length(unique(catches$serialnumber))){
    stop("Handle fractioned cacthes (delprøve)")
  }
  if (!all(catches$sampletype %in% c("20", "23"))){
    stop("Check if all sample types are acceptable")
  }

  tab <- merge(individuals, ages) #retain only aged individuals
  tab <- merge(catches, tab, by=names(catches)[names(catches) %in% names(tab)])


  tab$PSUid <- as.character(match(paste(tab$missiontype, tab$startyear, tab$platform, tab$missionnumber, tab$serialnumber), unique(paste(tab$missiontype, tab$startyear, tab$platform, tab$missionnumber, tab$serialnumber))))
  tab$SSUid <- as.character(match(paste(tab$PSUid, tab$catchsampleid), unique(paste(tab$PSUid, tab$catchsampleid))))
  tab$FishId <- as.character(match(paste(tab$SSUid, tab$specimenid), unique(paste(tab$SSUid, tab$specimenid))))
  tab$MeasurmentId <- match(paste(tab$FishId, tab$agedeterminationid), unique(paste(tab$FishId, tab$agedeterminationid)))

  if (length(unique(tab$MeasurmentId)) != length(tab$FishId)){
    stop("Deal with several age determinations for some fish")
  }

  if (any(tab$lengthmeasurement != "E")){
    stop("Some length measurements are not standard length")
  }

  if (any(tab$sampleproducttype != 1)){
    stop("Some samples does not have standard product type")
  }
  if (any(tab$catchproducttype != 1)){
    stop("Some catches does not have standard product type")
  }

  stopifnot(length(unique(tab$FishId)) == nrow(tab))

  tab$PSUselectionProb <- tab$catchweight / totalLanded
  tab$PSUinclusionProb <- lotteryEstimator::poissonInclusionProbability(tab$PSUselectionProb, length(unique(tab$PSUid)))
  tab$SSUselectionProb <- tab$lengthsampleweight / tab$catchweight
  tab$nSSU <- floor(tab$catchweight / tab$lengthsampleweight)
  tab$nFish <- tab$lengthsamplecount
  tab <- tab[,c("PSUid", "SSUid", "FishId", "PSUselectionProb", "PSUinclusionProb", "nSSU", "SSUselectionProb", "nFish", "age", "length")]

  return(tab)
}

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
