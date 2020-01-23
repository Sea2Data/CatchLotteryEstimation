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


context("Hansen-Hurwitz 2 stage estimator")
totalPSUs <- list()
totalPSUs$PSU1 <- hansenHurwitz(samplesPSU1, selectionProbabilitiesPSU1)
totalPSUs$PSU2 <- hansenHurwitz(samplesPSU1, selectionProbabilitiesPSU1)
total <- hansenHurwitz(totalPSUs, selectionProbabilities)

covarInter <- hansenHurwitzCovariance(totalPSUs, selectionProbabilities)
covarIntra <- hansenHurwitzIntra(covarPSUs, selectionProbabilities)
covar <- covarInter + covarIntra

expect_equal(dim(covar)[1], 2)
expect_equal(dim(covar)[2], 2)
expect_equal(colnames(covar), parnames)
expect_equal(rownames(covar), parnames)


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
