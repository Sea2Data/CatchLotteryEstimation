
context("countCategorical")
ages <- c(4,5,5,6,10)
categories <- 1:20
counts <- countCategorical(ages, categories)
expect_equal(length(counts), 20)
expect_equal(counts[["5"]], 2)
expect_equal(names(counts), as.character(1:20))

context("calculateSampleProportionCovariance")
ages <- c(.1,.2,.5,.2)
names(ages) <- 1:4
covariances <- calculateSampleProportionCovariance(ages)
expect_equal(names(ages), rownames(covariances))
expect_equal(names(ages), colnames(covariances))
expect_equal(covariances[1,1], .1 *(1 - .1))
expect_equal(covariances[3,2], -.5*.2)
