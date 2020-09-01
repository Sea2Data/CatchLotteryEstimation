context("nonresponse")
inclusionProbabilities <- c(.4,.1,.2,.1)
nonresponse <- c(F,F,T,T)
sampled <- inclusionProbabilities[!nonresponse]
correctionFactor <- randomNonResponseCorrectionFactor(inclusionProbabilities, sampled)
correctedInclusionProbabilites <- sampled * correctionFactor

sampleTotals <- c(3,4)

#estimate
preEst <- horvitzThompson(sampleTotals, correctedInclusionProbabilites)

#correct post-estimate
postEst <- horvitzThompson(sampleTotals, inclusionProbabilities[!nonresponse])/correctionFactor

expect_lt(abs(preEst - postEst), 1e-10)
