library(lotteryEstimator)
data("NSSH2019")

# Calculates the probability of selecting the
probWOreplacement <- function(selectionprob, selections, population){
  if (selections == 1){
    return(selectionprob)
  }
  else{
    estN <- selectionprob * population
    nextselp <- (estN - 1) / (population - 1)
    return(selectionprob * probWreplacement(nextselp, selections-1, population-1))
  }
}

NSSH2019$originalFishId <- NSSH2019$FishId

# Checks that the frequency of making a replicate sample of all equal elements is as expected
# for resampling with replacement
testFrequencyAllEqualWR <- function(size=3, n=3, iter=10, replicates=1000, tolerance=0.05){

  counts <- c()
  for (j in 1:replicates){
    allEq <- 0
    for (i in 1:iter){
      rs<-resample(NSSH2019[1:size,], c("FishId"), nSamples = c(n), replacement = c(T))
      if (length(unique(rs$originalFishId)) == 1){
        allEq <- allEq + 1
      }
    }
    counts <- c(counts, allEq)
  }

  expect <- iter * ((1/size) ** (n-1))

  meanResult <- mean(counts)
  if (abs(expect - meanResult)/expect > tolerance){
    warning(paste("Frequency of re-selection in selection is suspicious for sampling with replacement. Expected", expect, "mean result:", meanResult, "(",replicates,"replications)."))
  }
  else{
    write(paste("Expected", expect, "mean result:", meanResult, "(",replicates,"replications). Within RE:", tolerance),stdout())
  }

}


# Checks that the frequency of making a sampling a certain element repeatedly is as expected
# for resampling without replacement
# This property does not really hold for pseudopopulations when the sample size does not divide the pseudopopulation size
testFrequencyAllEqualWOR <- function(size=3, n=3, popsize=9, iter=1000, replicates=30, tolerance=0.05){

  NSSH2019$popsize <- popsize

  samp <- NSSH2019[1:size,]
  element <- samp$FishId[sample.int(size,1)]

  counts <- c()
  for (j in 1:replicates){
    allEq <- 0
    for (i in 1:iter){
      rs<-resample(samp, c("FishId"), nSamples = c(n), replacement = c(F), popSize = c("popsize"))
      if (all(rs$originalFishId == element)){
        allEq <- allEq + 1
      }
    }
    counts <- c(counts, allEq)
  }

  expect <- iter * probWOreplacement(1/size, n, popsize)

  meanResult <- mean(counts)
  if (abs(expect - meanResult)/expect > tolerance){
    warning(paste("Frequency of re-selection in selection is suspicious for sampling without replacement. Expected", expect, "mean result:", meanResult, "(",replicates,"replications)."))
  }
  else{
    write(paste("Expected", expect, "mean result:", meanResult, "(",replicates,"replications). Within RE:", tolerance),stdout())
  }
}
