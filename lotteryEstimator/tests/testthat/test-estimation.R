makeNamedVector <- function(vec, elementNames){
  names(vec) <- elementNames
  return(vec)
}

context("Hansen-Hurwitz estimator")

parnames <- c("param1", "param2")
samples <- list()
samples$firstsample <- makeNamedVector(c(10,100), parnames)
samples$secondsample <- makeNamedVector(c(5,400), parnames)
samples$thirdsample <- makeNamedVector(c(.1,.6), parnames)
selectionProbabilities <- c(1/10, 2/10, 1/100)

hh <- hansenHurwitz(samples, selectionProbabilities)

#test against manual calculation
manual <- makeNamedVector(c(10*10/3 + 5*5/3 + .1*100/3, 100*10/3 + 400*5/3 + .6*100/3), parnames)
expect_equal(hh[1], manual[1])
expect_equal(hh[2], manual[2])

context("Hansen-Hurwitz covariance estimator")
hhCovar <- hansenHurwitzCovariance(samples, selectionProbabilities)
expect_equal(dim(hhCovar)[1], 2)
expect_equal(dim(hhCovar)[2], 2)
expect_equal(colnames(hhCovar), parnames)
expect_equal(rownames(hhCovar), parnames)

#
# Test against univariate example at: https://online.stat.psu.edu/stat506/node/15/
#
samples <- list()
samples[[1]] <- 420
samples[[2]] <- 1785
samples[[3]] <- 2198
selectionProbabilities <- c(650/15650, 2840/15650, 3200/15650)
hh <- hansenHurwitz(samples, selectionProbabilities)
expect_lte(abs(hh - 10232.75), 1e-2)
hhCovar <- hansenHurwitzCovariance(samples, selectionProbabilities)
expect_lte((hhCovar[1] - 73125.74) / 73125.74, 0.001)



context("Hansen-Hurwitz intra-covariance estimator")
parnames <- c("param1", "param2")
samplesPSU1 <- list()
samplesPSU1$firstsample <- makeNamedVector(c(10,100), parnames)
samplesPSU1$secondsample <- makeNamedVector(c(5,400), parnames)
samplesPSU1$thirdsample <- makeNamedVector(c(.1,.6), parnames)
selectionProbabilitiesPSU1 <- c(1/10, 2/10, 1/100)

samplesPSU2 <- list()
samplesPSU2$firstsample <- makeNamedVector(c(9,800), parnames)
samplesPSU2$secondsample <- makeNamedVector(c(4,600), parnames)
samplesPSU2$thirdsample <- makeNamedVector(c(.2,.4), parnames)
selectionProbabilitiesPSU2 <- c(1/10, 1/10, 1/50)

selectionProbabilities <- c(1/20, 2/20)

covarPSUs <- list()
covarPSUs$PSU1 <- hansenHurwitzCovariance(samplesPSU1, selectionProbabilitiesPSU1)
covarPSUs$PSU2 <- hansenHurwitzCovariance(samplesPSU2, selectionProbabilitiesPSU2)

covarIntra <- hansenHurwitzIntra(covarPSUs, selectionProbabilities)
expect_equal(dim(covarIntra)[1], 2)
expect_equal(dim(covarIntra)[2], 2)
expect_equal(colnames(covarIntra), parnames)
expect_equal(rownames(covarIntra), parnames)

covarCrib <- ( (covarPSUs$PSU1 / selectionProbabilities[1]**2) + (covarPSUs$PSU2 / selectionProbabilities[2]**2)) / 4
expect_equal(covarCrib, covarIntra)


context("Hansen-Hurwitz 2 stage estimator")

#
# Test against univariate example in Lohr (ex 6.6)
#

samplesPSU1 <- list()
samplesPSU1[[1]] <- c(2)
samplesPSU1[[2]] <- c(3)
samplesPSU1[[3]] <- c(2.5)
samplesPSU1[[4]] <- c(3)
samplesPSU1[[5]] <- c(1.5)
selectionProbabilitiesPSU1 <- rep(1/24,5)

samplesPSU2 <- list()
samplesPSU2[[1]] <- c(2.5)
samplesPSU2[[2]] <- c(2)
samplesPSU2[[3]] <- c(3)
samplesPSU2[[4]] <- c(0)
samplesPSU2[[5]] <- c(0.5)
selectionProbabilitiesPSU2 <- rep(1/100,5)

samplesPSU3 <- list()
samplesPSU3[[1]] <- c(3)
samplesPSU3[[2]] <- c(0.5)
samplesPSU3[[3]] <- c(1.5)
samplesPSU3[[4]] <- c(2)
samplesPSU3[[5]] <- c(3)
selectionProbabilitiesPSU3 <- rep(1/100,5)

samplesPSU4 <-list()
samplesPSU4[[1]] <- c(1)
samplesPSU4[[2]] <- c(2.5)
samplesPSU4[[3]] <- c(3)
samplesPSU4[[4]] <- c(5)
samplesPSU4[[5]] <- c(2.5)
selectionProbabilitiesPSU4 <- rep(1/76,5)

samplesPSU5 <- list()
samplesPSU5[[1]] <- c(4)
samplesPSU5[[2]] <- c(4.5)
samplesPSU5[[3]] <- c(3)
samplesPSU5[[4]] <- c(2)
samplesPSU5[[5]] <- c(5)
selectionProbabilitiesPSU5 <- rep(1/44,5)

selectionProbabilities <- c(24,100,100,76,44)/647


totalPSUs <- list()
totalPSUs$PSU1 <- hansenHurwitz(samplesPSU1, selectionProbabilitiesPSU1)
totalPSUs$PSU2 <- hansenHurwitz(samplesPSU2, selectionProbabilitiesPSU2)
totalPSUs$PSU3 <- hansenHurwitz(samplesPSU3, selectionProbabilitiesPSU3)
totalPSUs$PSU4 <- hansenHurwitz(samplesPSU4, selectionProbabilitiesPSU4)
totalPSUs$PSU5 <- hansenHurwitz(samplesPSU5, selectionProbabilitiesPSU5)
total <- hansenHurwitz(totalPSUs, selectionProbabilities)
expect_equal(total, 1617.5) #checked against Lohr

covarPSUs <- list()
covarPSUs$PSU1 <- hansenHurwitzCovariance(samplesPSU1, selectionProbabilitiesPSU1)
covarPSUs$PSU2 <- hansenHurwitzCovariance(samplesPSU2, selectionProbabilitiesPSU2)
covarPSUs$PSU3 <- hansenHurwitzCovariance(samplesPSU3, selectionProbabilitiesPSU3)
covarPSUs$PSU4 <- hansenHurwitzCovariance(samplesPSU4, selectionProbabilitiesPSU4)
covarPSUs$PSU5 <- hansenHurwitzCovariance(samplesPSU5, selectionProbabilitiesPSU5)

semiManual <- sum(unlist(covarPSUs) / selectionProbabilities**2) / 5**2
covarIntra <- hansenHurwitzIntra(covarPSUs, selectionProbabilities) #checked against semi-manual calculation
expect_equal(covarIntra[1], semiManual)

covarInter <- hansenHurwitzCovariance(totalPSUs, selectionProbabilities)
expect_lte((covarInter[1] - 233.28*233.28) / covarInter[1], 0.0001) #checked against Lohr

covar <- covarInter + covarIntra

context("hierarchicalHansenHurwitz and hierarchicalHorwitzThompson")
numAtAgeSample <- function(sample){countCategorical(sample$age, 2:20)}
numAtAgeHaul <- function(sample){hierarchicalHorwitzThompson(sample, "SSUid", numAtAgeSample, "SSUinclusionProb")}
exampleSamples <- lotteryEstimator::NSSH2019
exampleSamples$SSUinclusionProb <- exampleSamples$SSUselectionProb * exampleSamples$nSSU
numAtAgeTotal <- hierarchicalHansenHurwitz(exampleSamples, "PSUid", numAtAgeHaul, "PSUselectionProb")
numAtAgeSimpleNSSH <- simpleNsshEstimatorReference(exampleSamples, minAge = 2, maxAge = 20)$catchAtAge
expect_equivalent(numAtAgeTotal, numAtAgeSimpleNSSH)

context("hierarchicalHansenHurwitzCovariance")
covarianceTotal <- hierarchicalHansenHurwitzCovariance(exampleSamples, "PSUid", numAtAgeHaul, function(x){0}, "PSUselectionProb")
covarianceSimpleNSSH <- simpleNsshEstimatorReference(exampleSamples, minAge = 2, maxAge = 20)$covariance
expect_equivalent(covarianceTotal, covarianceSimpleNSSH)


context("Horvitz-Thompson estimator")

parnames <- c("param1", "param2")
samples <- list()
samples$firstsample <- makeNamedVector(c(10,100), parnames)
samples$secondsample <- makeNamedVector(c(5,400), parnames)
samples$thirdsample <- makeNamedVector(c(.1,.6), parnames)
inclusionProbabilities <- c(1/10, 2/10, 1/10)

ht <- horvitzThompson(samples, inclusionProbabilities)
manual <- makeNamedVector(c(10*10/1 + 5*10/2 + .1*10/1, 100*10/1 + 400*10/2 + .6*10/1), parnames)
expect_equal(ht[1], manual[1])
expect_equal(ht[2], manual[2])


context("Horvitz-Thompson covariance estimator")

#testing against example at http://wiki.awf.forst.uni-goettingen.de/wiki/index.php/Horvitz-Thompson_estimator_example
samples <- list()
samples[[1]] <- 6
samples[[2]] <- 20

inclusionProbabilities <- c(2/3,2/3)
coInclusionProbabilities <- matrix(c(NA,1,
                                     1,NA), nrow=2)

ht <- horvitzThompson(samples, inclusionProbabilities)
expect_equal(ht, 39)

context("strata totals estimator")
mockstrata <- simpleNsshEstimatorReference(NSSH2019)
strata <- list()
strata$strata1 <- mockstrata
strata$strata2 <- mockstrata

result <- estimateFromStrataTotals(strata)
expect_true(all(result$catchAtAge == 2*mockstrata$catchAtAge))
expect_true(all(diag(result$covariance) == 2*diag(mockstrata$covariance)))
