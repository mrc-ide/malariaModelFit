# Create scraped likelihood function
likelihood_string <- function_to_string("src/likelihood.cpp")
usethis::use_data(likelihood_string, overwrite = TRUE)