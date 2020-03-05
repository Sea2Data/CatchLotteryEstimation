

#' Count categorical observations
#' @description
#'  Counts the number of observations of some category in a sample
#' @details
#'  Will count non-observed categories as zero.
#'  if 'categories' is NULL, the present values will be used as the complete set of categories
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

#' Proportions of categorical observations
#' @description
#'  Calculates the proportion of observations of some category in a sample
#' @details
#'  Will count non-observed categories as zero
#'  if 'categories' is NULL, the present values will be used as the complete set of categories
#' @param sample vector with sample observations for a categorical variable
#' @param categories the categories that was recorded if present in sample.
#' @return named vector with proporitons in each category
#' @examples
#'  ages <- c(4,5,5,6,10)
#'  categories <- 1:20
#'  proportionCategorical(ages, categories)
#' @export
proportionCategorical <- function(sample, categories=NULL){

  p <- countCategorical(sample, categories)
  p <- p / sum(p)
  return(p)

}

#' Age-length key
#' @description
#'  Calculates the proportion of ages in each length group.
#' @details
#'  Will count non-observed categories as zero
#'  'ages' and 'lengths' are observations from the same sample, so that observations of the two variables are paired by the index to these vectors.
#'  if 'ageCategories' or 'lengthCategories' is NULL, the present values will be used as the complete set of categories for the corresponding variable.
#' @param ages vector with age observations for a sample
#' @param lenghts vector with length observations for a sample
#' @param ageCategories the categories that was recorded for age, if present in sample.
#' @param lengthCategories the categories that was recorded for length if present in sample.
#' @return named matrix with proprotions of ages (rows) in each length group (columns). All rows sum to 1.
#' @examples
#'  ages <- c(4,5,5,6,10)
#'  lengths <- c("20-30", "30-40","40-50", "30-40","40-50")
#'  categoriesAges <- 1:20
#'  categoriesLengths <- c("0-10","10-20", "20-30","30-40","40-50")
#'  calculateAgeLengthKey(ages, lengths, categoriesAges, categoriesLengths)
#' @export
calculateAgeLengthKey <- function(ages, lenghts, ageCategories=NULL, lengthCategories=NULL){

  if (length(ages) != length(lenghts)){
    stop("'sampleVar1'' and 'sampleVar2' must be paried observations.")
  }

  if (is.null(ageCategories)){
    ageCategories <- unique(ages)
  }
  if (is.null(lengthCategories)){
    lengthCategories <- unique(lenghts)
  }

  m <- as.matrix(table(factor(ages, ageCategories), factor(lenghts, lengthCategories)))
  rownames(m) <- ageCategories
  colnames(m) <- lengthCategories
  m <- t(t(m) / (countCategorical(lenghts, lengthCategories)))
  m[is.nan(m)] <- 0

  return(m)

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

