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

context("hierarchicalHansenHurwitz and hierarchicalHorvitzThompson")
numAtAgeSample <- function(sample){countCategorical(sample$age, 2:20)}
numAtAgeHaul <- function(sample){hierarchicalHorvitzThompsonTotals(sample, "SSUid", numAtAgeSample, "SSUinclusionProb")}
exampleSamples <- lotteryEstimator::NSSH2019
exampleSamples$SSUinclusionProb <- exampleSamples$SSUselectionProb #only one SSU sampled pr PSU in example
numAtAgeTotal <- hierarchicalHansenHurwitzTotals(exampleSamples, "PSUid", numAtAgeHaul, "PSUselectionProb")
numAtAgeSimpleNSSH <- simpleNsshEstimatorReference(exampleSamples, minAge = 2, maxAge = 20)$catchAtAge
expect_equivalent(numAtAgeTotal, numAtAgeSimpleNSSH)

context("hierarchicalHansenHurwitzCovariance")
covarianceTotal <- hierarchicalHansenHurwitzCovariance(exampleSamples, "PSUid", numAtAgeHaul, function(x){0}, "PSUselectionProb")
covarianceSimpleNSSH <- simpleNsshEstimatorReference(exampleSamples, minAge = 2, maxAge = 20)$covariance
expect_equivalent(covarianceTotal, covarianceSimpleNSSH)

context("hierarchicalHorvitzThompsonCovariance")
numAtAgeSample <- function(sample){countCategorical(sample$age, 2:20)}
numAtAgeHaul <- function(sample){hierarchicalHorvitzThompsonTotals(sample, "SSUid", numAtAgeSample, "SSUinclusionProb")}
exampleSamples <- lotteryEstimator::NSSH2019
exampleSamples$SSUinclusionProb <- exampleSamples$SSUselectionProb #only one SSU sampled pr PSU in example
sampleNames <-  exampleSamples$PSUid[!duplicated(exampleSamples$PSUid)]
incProb <- exampleSamples$PSUinclusionProb[!duplicated(exampleSamples$PSUid)]
names(incProb) <- sampleNames
coInclusionMatrix <- outer(incProb, incProb)
colnames(coInclusionMatrix) <- sampleNames
rownames(coInclusionMatrix) <- sampleNames
diag(coInclusionMatrix) <- NA
covarianceTotal <- hierarchicalHorvitzThompsonCovariance(exampleSamples, "PSUid", numAtAgeHaul, function(x){outer(numAtAgeHaul(x)*0, numAtAgeHaul(x)*0)}, "PSUinclusionProb", coInclusionMatrix)

expect_true(all(diag(covarianceTotal)>=0))

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


context("Horvitz-Thompson covariance estimator, data cheks")
samples <- list()
samples[["2"]] <- 6
samples[["1"]] <- 20

inclusionProbabilities <- c((1-(1-.1)**2),(1-(1-.3)**2))
names(inclusionProbabilities) <- c("1", "2")
coInclusionProbabilities <- matrix(c(NA,sum(inclusionProbabilities) - (1-(1-.1-.3)**2),
                                     sum(inclusionProbabilities) - (1-(1-.3-.1)**2),NA), nrow=2)
rownames(coInclusionProbabilities) <- c("1", "2")
colnames(coInclusionProbabilities) <- c("1", "2")
expect_error(horvitzThompsonCovariance(samples, inclusionProbabilities, coInclusionProbabilities), "Vector names does not match between sampleTotals and inclusionProbabilities")
names(inclusionProbabilities) <- c("2", "1")
expect_error(horvitzThompsonCovariance(samples, inclusionProbabilities, coInclusionProbabilities), "Vector names does not match between sampleTotals and coInclusionProbabilities")

context("Horvitz-Thompson covariance estimator")

#testing against example at http://wiki.awf.forst.uni-goettingen.de/wiki/index.php/Horvitz-Thompson_estimator_example
samples <- list()
samples[["1"]] <- 6
samples[["2"]] <- 20

inclusionProbabilities <- c((1-(1-.1)**2),(1-(1-.3)**2))
names(inclusionProbabilities) <- c("1", "2")
coInclusionProbabilities <- matrix(c(NA,sum(inclusionProbabilities) - (1-(1-.1-.3)**2),
                                     sum(inclusionProbabilities) - (1-(1-.3-.1)**2),NA), nrow=2)
rownames(coInclusionProbabilities) <- c("1", "2")
colnames(coInclusionProbabilities) <- c("1", "2")

ht <- horvitzThompson(samples, inclusionProbabilities)
htTrue <- 70.79
expect_true(abs(ht-htTrue) < .005)
covMono <- horvitzThompsonCovariance(samples, inclusionProbabilities, coInclusionProbabilities)

covMonoTrue <- (1-inclusionProbabilities[1])*6**2 / inclusionProbabilities[1]**2
covMonoTrue <- covMonoTrue + (1-inclusionProbabilities[2])*20**2 / inclusionProbabilities[2]**2
covMonoTrue <- covMonoTrue + (coInclusionProbabilities[1,2] - inclusionProbabilities[1]*inclusionProbabilities[2]) * 6*20 / (inclusionProbabilities[1]*inclusionProbabilities[2]*coInclusionProbabilities[1,2])
covMonoTrue <- covMonoTrue + (coInclusionProbabilities[2,1] - inclusionProbabilities[2]*inclusionProbabilities[1]) * 20*6 / (inclusionProbabilities[2]*inclusionProbabilities[1]*coInclusionProbabilities[2,1])
expect_equal(covMono[[1]], covMonoTrue[[1]])

context("Horvitz-Thompson covariance estimator, multivariate")
samples <- list()
samples[["1"]] <- c(6,6)
samples[["2"]] <- c(20,20)

inclusionProbabilities <- c((1-(1-.1)**2),(1-(1-.3)**2))
names(inclusionProbabilities) <- c("1", "2")
coInclusionProbabilities <- matrix(c(NA,sum(inclusionProbabilities) - (1-(1-.1-.3)**2),
                                     sum(inclusionProbabilities) - (1-(1-.1-.3)**2),NA), nrow=2)
rownames(coInclusionProbabilities) <- c("1", "2")
colnames(coInclusionProbabilities) <- c("1", "2")

ht <- horvitzThompson(samples, inclusionProbabilities)
expect_true(all(abs(ht-htTrue) < .005))
cov <- horvitzThompsonCovariance(samples, inclusionProbabilities, coInclusionProbabilities)
expect_true(all((diag(cov)-as.vector(covMono)) < 1e-11))

context("Horvitz-Thompson covariance estimator, multivariate covar")
samples <- list()
samples[["1"]] <- c(6,3)
samples[["2"]] <- c(20,10)

inclusionProbabilities <- c((1-(1-.1)**2),(1-(1-.3)**2))
names(inclusionProbabilities) <- c("1", "2")
coInclusionProbabilities <- matrix(c(NA,sum(inclusionProbabilities) - (1-(1-.1-.3)**2),
                                     sum(inclusionProbabilities) - (1-(1-.1-.3)**2),NA), nrow=2)
rownames(coInclusionProbabilities) <- c("1", "2")
colnames(coInclusionProbabilities) <- c("1", "2")

ht <- horvitzThompson(samples, inclusionProbabilities)
expect_true(abs(ht[1]-htTrue) < .005)
cov <- horvitzThompsonCovariance(samples, inclusionProbabilities, coInclusionProbabilities)
expect_equal(cov[1,2], cov[2,1])
browser()
expect_equal(cov[1,1], covMono[1,1])

context("hierarchical strata estimate")
numAtAgeStrata <- function(sample){proportionCategorical(sample$age, 2:20)*sample$lengthStrataTotal[1]}
numAtAgeSample <- function(sample){hierarchicalStratifiedTotals(sample, "lengthStrata", numAtAgeStrata)}
ssuTotal <- numAtAgeSample(NSSH2019Stratified[NSSH2019Stratified$SSUid=="1",])
expect_equal(ssuTotal[["4"]],1)
expect_equal(ssuTotal[["6"]],6+15+1+6)

context("hierarchical strata covariance estimate")
ageCovarianceStrata <- function(sample){calculateSampleProportionCovariance(proportionCategorical(sample$age, 2:20)) * (sample$lengthStrataTotal[1]**2 / ((nrow(sample))*(nrow(sample)-1)))}
ageCovarianceSample <- function(sample){hierarchicalStratifiedCovariance(sample, "lengthStrata", ageCovarianceStrata)}
ssuCov <- ageCovarianceSample(NSSH2019Stratified[NSSH2019Stratified$SSUid=="1",])
expect_equal(ssuCov[["4","4"]],.5)

context("Test integer IDs")
numAtAgeSample <- function(sample){countCategorical(sample$age, 2:20)}
numAtAgeHaul <- function(sample){hierarchicalHorvitzThompsonTotals(sample, "SSUid", numAtAgeSample, "SSUinclusionProb")}

exampleSamples <- lotteryEstimator::NSSH2019
exampleSamples$PSUid <- as.integer(exampleSamples$PSUid)+1
exampleSamples$SSUinclusionProb <- exampleSamples$SSUselectionProb #only one SSU sampled pr PSU in example
expect_error(hierarchicalHansenHurwitzTotals(exampleSamples, "PSUid", numAtAgeHaul, "PSUselectionProb"), "partitionId must be a character")
expect_error(hierarchicalHorvitzThompsonTotals(exampleSamples, "PSUid", numAtAgeHaul, "PSUselectionProb"), "partitionId must be a character")

context("Test domain estimation HH")
numAtAgeSample <- function(sample){countCategorical(sample$age, 2:20)}
numAtAgeHaul <- function(sample){hierarchicalHorvitzThompsonTotals(sample, "SSUid",
                                                                   numAtAgeSample, "SSUinclusionProb")}
exampleSamples <- lotteryEstimator::NSSH2019
exampleSamples$SSUinclusionProb <- exampleSamples$SSUselectionProb * exampleSamples$nSSU
numAtAgeTotal <- function(sample){hierarchicalHansenHurwitzTotals(sample, "PSUid",
                                                                  numAtAgeHaul, "PSUselectionProb")}
exampleSamples$domains <- cut(as.numeric(exampleSamples$PSUid), 3, labels = c("one", "two", "three"))
de1 <- hierarchicalHansenHurwitzDomainTotals(exampleSamples, "domains", numAtAgeTotal, "one", 1/3)
expect_equal(length(de1), 19)

context("Test domain covariance estimation HH")
numAtAgeSample <- function(sample){countCategorical(sample$age, 2:20)}
numAtAgeHaul <- function(sample){hierarchicalHorvitzThompsonTotals(sample, "SSUid",
                                                                   numAtAgeSample, "SSUinclusionProb")}
exampleSamples <- lotteryEstimator::NSSH2019
exampleSamples$SSUinclusionProb <- exampleSamples$SSUselectionProb * exampleSamples$nSSU
hhCov <- function(samples){hierarchicalHansenHurwitzCovariance(samples, "PSUid",
                                                               numAtAgeHaul, function(x){0},
                                                               "PSUselectionProb")}
#define som arbitrary domains for the sake of the example
exampleSamples$domains <- cut(as.numeric(exampleSamples$PSUid), 3, labels = c("one", "two", "three"))
#obtain domain estimate
cov <- hierarchicalHansenHurwitzDomainTotalCovariance(exampleSamples, "domains", hhCov, "two", 1/3)
expect_equal(ncol(cov), 19)
expect_equal(nrow(cov), 19)



context("Test domain estimation HT")
numAtAgeSample <- function(sample){countCategorical(sample$age, 2:20)}
numAtAgeHaul <- function(sample){hierarchicalHorvitzThompsonTotals(sample, "SSUid",
                                                                   numAtAgeSample, "SSUinclusionProb")}
exampleSamples <- lotteryEstimator::NSSH2019
exampleSamples$SSUinclusionProb <- exampleSamples$SSUselectionProb * exampleSamples$nSSU
numAtAgeTotal <- function(sample){hierarchicalHorvitzThompsonTotals(sample, "PSUid",
                                                                  numAtAgeHaul, "PSUinclusionProb")}
exampleSamples$domains <- cut(as.numeric(exampleSamples$PSUid), 3, labels = c("one", "two", "three"))
de1 <- hierarchicalHansenHurwitzDomainTotals(exampleSamples, "domains", numAtAgeTotal, "one", 1/3)
expect_equal(length(de1), 19)

