#' Example data - Norwegian Spring Spawning Herring sampled in 2019
#' The example has length stratified age samples taken from each haul-sample.
#' The data was not length stratified when colleted, but has been resampled after the fact for purposes of testing and illustration.
#'
#' Contains the same columns as \code{\link[lotteryEstimator]{NSSH2019}} in addition to the columns 'lengthStrata' and 'lengthStrataTotal'
#'
#' @docType data
#'
#' @usage data(NSSH2019Stratified)
#'
#' @format \code{\link[data.table]{data.table}} one row for each measured fish and with columns:
#'  \describe{
#'     \item{PSUid}{Identifier for the primary sampling unit (PSU) the fish was sampled from}
#'     \item{SSUid}{Identifier for the secondary sampling unit (SSU) the fish was sampled from. Identifier is unique across all PSUs}
#'     \item{Fishid}{Identifier for the fish sampled. Identifier is unique across all PSUs and SSUs}
#'     \item{PSUselectionProb}{Selection probability for the PSU}
#'     \item{nSSU}{The total number of SSUs in the population for this PSU}
#'     \item{SSUselectionProb}{Selection probability for the SSU}
#'     \item{nFish}{The total number of fish in the population for this SSU}
#'     \item{age}{The age determied for the fish}
#'     \item{length}{The length determined for the fish}
#'     \item{lengthStrata}{The length strata for the fish}
#'     \item{lengthStrataTotal}{The total number of fish in the haul-sample (length-measured fish)}
#'  }
#'
#' @keywords datasets
#'
#' @examples
#' data(NSSH2019Stratified)
#'
"NSSH2019Stratified"
