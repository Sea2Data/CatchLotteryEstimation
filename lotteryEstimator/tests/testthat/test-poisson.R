context("Poisson sampling")
catch=200
TAC=400000
sampleSize=200
incSingle <- poissonInclusionProbability(catch/TAC, 1)
expect_equal(incSingle, catch/TAC)
inc <- poissonInclusionProbability(catch/TAC, sampleSize)
expect_gt(inc, incSingle)
sel <- poissonSelectionProbability(inc, sampleSize)
expect_equal(sel, catch/TAC)

vecInc <- poissonInclusionProbability(c(.0005,0.0001),300)
expect_true(all(vecInc>0))
expect_true(all(vecInc<1))
vecSel <- poissonSelectionProbability(vecInc, 300)
expect_true(all(abs(vecSel - c(.0005,0.0001)) < 1e-15))

vecPairInc <- poissonJointInclusionProbabilityMatrix(c(.0005,0.0001))
expect_equal(vecPairInc[1,2], vecPairInc[2,1])
expect_true(all(is.na(diag(vecPairInc))))

corrected <- randomNonResponseCorrection(c(.4,.1), 4)
expect_equal(length(corrected), 2)
expect_lt(abs(corrected[1] - (1-(1-.4)**.5)), 1e-10)

corrected <- randomNonResponseCorrection(.1, 100, 80)
expect_equal(length(corrected), 1)
expect_lt(abs(corrected - (1-(1-.1)**.8)), 1e-10)


sample <- poissionSample(1:1000, 1:1000/1000, 100)
expect_gt(length(sample), 0)

expect_error(poissionSample(1:1000, 1:1000, 100))
