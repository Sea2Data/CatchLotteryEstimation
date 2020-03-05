library(lotteryEstimator)

plotResults <- function(estimates, title=""){
  caa <- estimates$catchAtAge / 1e6
  se <- sqrt(diag(estimates$covariance)) / 1e6
  ci <- se*1.96
  plot(1:length(caa), caa, xlab="Age", ylab="catch at age (millions)", main=title, ylim=c(0,max(caa+ci)), xaxt="n")
  segments(1:length(caa), caa+ci, y1 = caa-ci)
  axis(side=1, 1:length(caa), names(caa))
}

data(NSSH2019)
estimates <- simpleNsshEstimator(NSSH2019, 1,20)
plotResults(estimates, "NSSH 2019 simple")

estimates2st <- twoStageNsshEstimator(NSSH2019, 1,20)
plotResults(estimates2st, "NSSH 2019 2stage")

estimates2stBoot <- twoStageNsshEstimatorBootstrapped(NSSH2019, 1,20,10)
plotResults(estimates2stBoot, "NSSH 2019 2stageBs")

str <- NSSH2019
str$lengthStrata <- "all"
estimates2stStratified <- twoStageNsshEstimatorStratified(str, 1,20)
plotResults(estimates2stStratified, "NSSH 2019 2stage single strata")

estimates2stStratifiedAct <- twoStageNsshEstimatorStratified(NSSH2019Stratified, 1,20)
plotResults(estimates2stStratifiedAct, "NSSH 2019 2stage stratified")

hist((diag(estimates$covariance) - diag(estimates2st$covariance)) / sum(diag(estimates$covariance)))
hist((diag(estimates$covariance) - diag(estimates2stBoot$covariance)) / sum(diag(estimates$covariance)))
hist((diag(estimates$covariance) - diag(estimates2stStratified$covariance)) / sum(diag(estimates$covariance)))
