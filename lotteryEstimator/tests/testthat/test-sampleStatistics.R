
context("countCategorical")
ages <- c(4,5,5,6,10)
categories <- 1:20
counts <- countCategorical(ages, categories)
expect_equal(length(counts), 20)
expect_equal(counts[["5"]], 2)
expect_equal(names(counts), as.character(1:20))
