#' SRS joint inclusion probability
#' @description calculates the pariwise joint inclusion probability for simple random sampling with replacement, based on sample size and population size
#' @details For more convenient incorporation in generic estimators, the result is return as a n X n matrix with the joint inclusion probability repeated
#' @param n sample size
#' @param N population size
#' @return the pairwise joint inclusion probability matrix (n X n)
#' @export
srsworJointInclusionProbabilityMatrix <- function(n,N){

  if (n > N){
    stop("sample size (n) cannot exceed population size (N)")
  }
  if (n < 1){
    stop("sample size (n) must be positive number")
  }
  if (N < 1){
    stop("population size (N) must be positive number")
  }

  coinc <- n * (n-1) / ( (N * (N-1) ))
  coinc <- matrix(coinc, n, n)
  diag(coinc) <- NA

  return(coinc)
}
