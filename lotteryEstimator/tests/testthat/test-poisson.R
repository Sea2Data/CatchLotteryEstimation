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
