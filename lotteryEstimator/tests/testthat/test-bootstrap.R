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
