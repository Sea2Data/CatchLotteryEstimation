library(RstoxData)

#
# Functions for preparing data for package lotteryEstimator
#

#' prepares an example data set
#' @param datafile biotic 3 annual data file
#' @param species species code to extract
#' @param totalLanded total landed weight of species/stock (Kg), used to calculate normalized selection probabilites
#' @return data.frame() with one row per measured fish
#' @noRd
extract_set <- function(datafile="~/bioticsets/v3/biotic_cruiseNumber_19-2019-1_Kommersielle+fartøy_2020-01-22T23.02.15.698Z.xml", species="161722.G03", totalLanded=430506000, removemissingCatchWeight=T, removeMissingAges=T){
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


  tab$PSUid <- match(paste(tab$missiontype, tab$startyear, tab$platform, tab$missionnumber, tab$serialnumber), unique(paste(tab$missiontype, tab$startyear, tab$platform, tab$missionnumber, tab$serialnumber)))
  tab$SSUid <- match(paste(tab$PSUid, tab$catchsampleid), unique(paste(tab$PSUid, tab$catchsampleid)))
  tab$FishId <- match(paste(tab$SSUid, tab$specimenid), unique(paste(tab$SSUid, tab$specimenid)))
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
  tab$SSUselectionProb <- tab$lengthsampleweight / tab$catchweight
  tab <- tab[,c("PSUid", "SSUid", "PSUselectionProb", "SSUselectionProb", "age", "length")]

  return(tab)
}
