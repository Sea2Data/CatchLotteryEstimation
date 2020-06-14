context("srs functions")

vecPairInc <- srsworJointInclusionProbabilityMatrix(2,10)
expect_equal(vecPairInc[1,2], vecPairInc[2,1])
expect_true(all(is.na(diag(vecPairInc))))
