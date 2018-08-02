context("test-parameters.R")

#------------------------------------------------
test_that("parameter_definitions and all_params contain the same parameter sets, i.e. every parameter is documented", {
  names1 <- sort(names(parameter_definitions()))
  names2 <- sort(all_params())
  expect_equal(names1, names2)
})

#------------------------------------------------
test_that("default parameter set contains all required parameters", {
  names1 <- sort(names(default_parameters()))
  names2 <- sort(all_params())
  expect_equal(names1, names2)
})