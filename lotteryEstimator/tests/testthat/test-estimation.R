makeNamedVector <- function(vec, elementNames){
  names(vec) <- elementNames
  return(vec)
}

context("Hansen-Hurwitz estimator")

parnames <- c("param1", "param2")
samples <- list()
samples$firstsample <- makeNamedVector(c(10,100), parnames)
samples$secondsample <- makeNamedVector(c(5,400), parnames)
samples$thirdsample <- makeNamedVector(c(.1,.6), parnames)
selectionProbabilites <- c(1/10, 2/10, 1/100)

hh <- hansenHurwitz(samples, selectionProbabilites)

#test against manual calculation
manual <- makeNamedVector(c(10*10/3 + 5*5/3 + .1*100/3, 100*10/3 + 400*5/3 + .6*100/3), parnames)
expect_equal(hh[1], manual[1])
expect_equal(hh[2], manual[2])

hhCovar <- hansenHurwitzInter(samples, selectionProbabilites)

