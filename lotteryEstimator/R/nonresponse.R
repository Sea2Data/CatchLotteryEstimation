
#' Correct total random nonrespons
#' @description
#'  Horvitz-Thompson estimate of a total for random non-response, assuming known response rate
#' @param respTotal Estimate of total based on respondents
#' @param responseRate the response rate.
#' @return Estimate of the total, corrected for random nonresponse
#' @examples
#'  #survey with an estimated total of 55, and a responserate of .7
#'  randomNonResponseCorrectionTotal(55, .7)
#' @export
randomNonResponseCorrectionTotal <- function(respTotal, responseRate){
  nonrespTot <- respTotal * (1-responseRate) / responseRate
  return(respTotal+nonrespTot)
}

#' Correct total random nonrespons
#' @description
#'  Horvitz-Thompson estimate of a standard error for random non-response, assuming known response rate.
#' @param respSE Estimate of standard error of total based on respondents
#' @param responseRate the response rate.
#' @return Estimate of the standard error of the total, corrected for random nonresponse
#' @examples
#'  #survey with an estimated SE of 30, and a responserate of .7
#'  randomNonResponseCorrectionSE(30, .7)
#' @export
randomNonResponseCorrectionSE <- function(respSE, responseRate){

  return(respSE + respSE*((1-responseRate) / responseRate))
}
