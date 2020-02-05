#' Count categorical observations
#' @description
#'  Counts the number of observations of some category in a sample
#' @details
#'  Will count non-observed categories as zero
#' @param sample vector with sample observations for a categorical variable
#' @param categories the categories that was recorded if present in sample.
#' @return named vector with counts in each category
#' @examples
#'  ages <- c(4,5,5,6,10)
#'  categories <- 1:20
#'  countCategorical(ages, categories)
#' @export
countCategorical <- function(sample, categories=NULL){

  if (is.null(categories)){
    categories <- unique(sample)
  }

  a <- tabulate(factor(sample, levels=categories), nbins=length(categories))
  names(a) <- categories
  return(a)
}

#' Covariance of sample proprotions
#' @param proportions vector of size 'n' containging sampled proportions for a categorical variable with 'n' levels
#' @return matrix() 'n'x'n' matrix with proportion covariances.
#' @examples
#'  ages <- c(4,5,5,6,10)
#'  categories <- 1:20
#'  sampleCounts <- countCategorical(ages, categories)
#'  sampleProrportions <- sampleCounts / sum(sampleCounts)
#'  covariance <- calculateSampleProportionCovariance(sampleProrportions)
#' @export
calculateSampleProportionCovariance <- function(proportions){
  return(diag(proportions) + outer(-proportions, proportions))
}

