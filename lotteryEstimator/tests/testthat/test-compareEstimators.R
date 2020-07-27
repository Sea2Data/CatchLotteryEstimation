context("test sampleExamplePopulation")

sample <- sampleExamplePopulation(10,10)
expect_true(all(sample$nPSU <= sample$Npsu))
expect_true(all(sample$nSSU <= sample$Nssu))
expect_true(all(sample$spPSU>=0))
expect_true(all(sample$spSSU>=0))
expect_true(all(sample$spSSU >= 1 / sample$Nssu))

sampleAll <- sampleExamplePopulation(72,365)
expect_equal(nrow(longlinerPopulation),nrow(sampleAll))
