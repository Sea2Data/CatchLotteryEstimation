context("sample indecies")
ind <- sampleIndecies(5,9,F,100)
expect_equal(length(ind), 5)
expect_lte(max(ind), 9)
expect_error(sampleIndecies(5,10,F,NA))
expect_error(sampleIndecies(5,10,F,9))

ind <- sampleIndecies(5,10,T,NA)
expect_equal(length(ind), 5)
expect_lte(max(ind), 10)

#sampling entire population
ind <- sampleIndecies(10,10,F,10)
expect_equal(sort(ind), 1:10)

context("resample")
#sampling without replacement and population set to sample sizes
rs<-resample(NSSH2019, c("PSUid", "SSUid", "FishId"),replacement = c(F,F,F))
expect_equal(nrow(rs), nrow(NSSH2019))
expect_equal(sum(rs$age), sum(NSSH2019$age))

rs<-resample(NSSH2019, c("PSUid", "SSUid"),replacement = c(F,F))
expect_equal(nrow(rs), nrow(NSSH2019))
expect_equal(sum(rs$age), sum(NSSH2019$age))

NSSH2019$ssize1 <- 2
NSSH2019$ssize2 <- 1
rs<-resample(NSSH2019, c("PSUid", "SSUid", "FishId"),replacement = c(T,T,T), nSamples = c("ssize1","ssize2","ssize1"))
expect_equal(nrow(rs), 4)

context("bootstrap")
estimator <- function(sample){countCategorical(sample$age, 1:19)}
bs<-bootstrap(NSSH2019, estimator, 10, hierarchy = c("FishId"), replacement=c(TRUE))
expect_equal(bs$iterations, 10)
expect_equal(length(bs$meanEstimate), 19)
expect_equal(nrow(bs$covariances), 19)
expect_equal(ncol(bs$covariances), 19)
expect_true(all(diag(bs$covariances)>=0))
expect_true(all(bs$meanEstimate >= 0))

estimator <- function(sample){countCategorical(sample$age, 1:19)}
NSSH2019$PSUrepl <- 2
bs<-bootstrap(NSSH2019, estimator, 10, hierarchy = c("SSUid", "FishId"), nSamples = c("PSUrepl",NA), popSize = c("nSSU", NA), replacement=c(FALSE, TRUE))
expect_equal(bs$iterations, 10)
expect_equal(length(bs$meanEstimate), 19)
expect_equal(nrow(bs$covariances), 19)
expect_equal(ncol(bs$covariances), 19)
expect_true(all(diag(bs$covariances)>=0))
expect_true(all(bs$meanEstimate >= 0))
