library(lotteryEstimator)

plotResults <- function(estimates, title=""){
  caa <- estimates$catchAtAge / 1e6
  se <- sqrt(diag(estimates$covariance)) / 1e6
  ci <- se*1.96
  plot(names(caa), caa, xlab="Age", ylab="catch at age (millions)", main=title, ylim=c(0,max(caa+ci)))
  segments(1:length(caa), caa+ci, y1 = caa-ci)
}

data(NSSH2019)
estimates <- simpleNsshEstimator(NSSH2019, 1,20)

plotResults(estimates, "NSSH 2019")
