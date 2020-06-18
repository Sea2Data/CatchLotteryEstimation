context("test sampleExamplePopulation")

sample <- sampleExamplePopulation(10,10)
expect_true(all(sample$nPSU <= sample$Npsu))
expect_true(all(sample$nSSU <= sample$Nssu))


sampleAll <- sampleExamplePopulation(72,365)
expect_equal(nrow(examplePopulation),nrow(sampleAll))
