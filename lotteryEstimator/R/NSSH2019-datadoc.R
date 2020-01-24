#' Example data - Norwegian Spring Spawning Herring sampled in 2019
#'
#' @docType data
#'
#' @usage data(NSSH2019)
#'
#' @format \code{\link[data.table]{data.table}} one row for each measured fish and with columns:
#'  \describe{
#'     \item{PSUid}{Identifier for the primary sampling unit (PSU) the fish was sampled from}
#'     \item{SSUid}{Identifier for the secondary sampling unit (SSU) the fish was sampled from. Identifier is unique across all PSUs}
#'     \item{Fishid}{Identifier for the fish sampled. Identifier is unique across all PSUs and SSUs}
#'     \item{PSUselectionProb}{Selection probability for the PSU}
#'     \item{nSSU}{The total number of SSUs in the population for this PSU}
#'     \item{SSUselectionProb}{Selection probability for the SSU}
#'     \item{age}{The age determied for the fish}
#'     \item{length}{The length determined for the fish}
#'  }
#'
#' @keywords datasets
#'
#' @examples
#' data(NSSH2019)
#'
"NSSH2019"
