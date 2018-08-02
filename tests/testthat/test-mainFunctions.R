context("test-mainFunctions.R")

#------------------------------------------------
test_that("human_equilibrium_no_het produces correct results", {
  
  # define parameters
  EIR <- 1
  ft <- 0.4
  age <- default_age()
  p <- load_parameters("parameters_Griffin2014.txt")
  
  # calculate equilibrium solution
  eq1 <- human_equilibrium_no_het(EIR, ft, p, age)
  
  # load the same results from file
  name <- malariaModelFit_file("test_comparisons/human_equilibrium_no_het_test1.RDS")
  saved_solution <- readRDS(name)
  
  # compare solutions
  expect_equal(eq1, saved_solution)
})

#------------------------------------------------
test_that("human_equilibrium produces correct results", {
  
  # define parameters
  EIR <- 1
  ft <- 0.4
  age <- default_age()
  h <- gq_normal(5)
  p <- load_parameters("parameters_Griffin2014.txt")
  
  # calculate equilibrium solution
  eq1 <- human_equilibrium(EIR, ft, p, age, h)
  
  # load the same results from file
  name <- malariaModelFit_file("test_comparisons/human_equilibrium_test1.RDS")
  saved_solution <- readRDS(name)
  
  # compare solutions
  expect_equal(eq1, saved_solution)
})