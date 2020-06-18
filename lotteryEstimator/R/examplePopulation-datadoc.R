#' Example data - A population of longline vessels
#'
#' A fictous census of fishing operations with longline.
#' Each row correspond to catch of one species.
#'
#' The data set is compiled for the purpose of testing sampling and estimation routines,
#' variables denoting time, space or vessels are not encoded in any interpretable way.
#'
#' @docType data
#'
#' @usage data(longlinerPopulation)
#'
#' @format \code{\link[data.table]{data.table}} one row for each measured fish and with columns:
#'  \describe{
#'     \item{vesselId}{code for vessel}
#'     \item{opid}{code for fishing operation}
#'     \item{timeid}{code identifying time of operation}
#'     \item{spaceid}{code identifying location of operation}
#'     \item{spacetimeid}{code time and location of operation}
#'     \item{catchid}{code for catch of a specific species}
#'     \item{gearFAO}{FAO code for gear (ISSCFG 1980)}
#'     \item{hookCount}{The total number of hooks on the line}
#'     \item{soaktimeHours}{The total fishing time in hours}
#'     \item{speciesFAO}{FAO code for species (ASFIS 3-apha code)}
#'     \item{wholeWeightKg}{weight of catch in kg}
#'  }
#'
#' @keywords datasets
#'
#' @examples
#' data(longlinerPopulation)
#'
"longlinerPopulation"
