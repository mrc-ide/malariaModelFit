#' Prepare data
#' 
#' Prepare (standardised) data inputs for use with Equilibrium/likelihood master function
#'
#' @param data Standardised data input
#' @param eta eta parameter
#' @param rho rho parameter
#' @param a0 A0 parameter
#' @param output_type Ouput type: "ll" = Log likelihood, "eq" = Equilibrium
#' @param nh Number of heterogeneity groups
#'
#' @return Vector of data in correct format
#' @export
prepare_data <- function(data, eta = 0.0001304631, rho = 0.85, a0 = 2920, output_type = "ll", nh = 5){
  ot <- 0
  if(output_type == "eq"){
    ot <- 1
  }
  # Split data by site
  site_list <- split(data, data$site_index)
  # Study N
  study_n <- length(unqiue(study_index))
  study <- sapply(site_list, function(x){
    x$study_index[1]
  })
  # Site N and size
  site_n <- length(site_list)
  group_n <- sapply(site_list, nrow) + 1
  # Age inputs
  ages <- lapply(site_list, function(x){
    age_prep(c(x$age0[1], x$age1), eta = eta, rho = rho, a0 = a0)
  })
  prop <- unlist(lapply(ages, "[[", "prop"))
  r <- unlist(lapply(ages, "[[", "r"))
  age_days_midpoint <- unlist(lapply(ages, "[[", "age_days_midpoint"))
  psi <- unlist(lapply(ages, "[[", "psi"))
  age20 <- unlist(lapply(ages, "[[", "age20"))
  # Gaussian quandrature nodes and weights
  het_d <- c(nh, unlist(gq_normal(nh)))
  # Numerator
  numer <- data$numer
  # Denominator
  denom <- data$denom
  # Prevalence or incidence
  type = data$type
  # Case detection type
  case_detection <- data$case_detection
  # Age random effect group
  age_bracket <- data$age_bracket
  
  # Return output vector in specific format
  return(list(output_type = ot,
              study_n = study_n,
              study = study, 
              site_n = site_n, group_n = group_n,
              prop = prop, r = r, age_days_midpoint = age_days_midpoint,
              psi = psi, age20 = age20, het_d = het_d, 
              numer = numer, denom = denom, type = type,
              case_detection = case_detection,
              age_bracket = age_bracket))
}

#' Pre-calculation of equilibrium age-related values
#'
#' @param age Vector of age groups
#' @param eta Parameter eta
#' @param rho Parameter rho
#' @param a0 Parameter a0
#'
#' @return A named list of equilibrium age-related values
age_prep <- function(age, eta, rho, a0){
  n_age <- length(age)
  age_days <- age * 365
  prop <- rep(0, n_age)
  r <- rep(0, n_age)
  for (i in 1:n_age) {
    if (i == n_age) {
      r[i] <- 0
    } else {
      age_width <- age_days[i+1] - age_days[i]
      r[i] <- 1/age_width
    }
    if (i == 1) {
      prop[i] <- eta/(r[i] + eta)
    } else {
      prop[i] <- prop[i-1]*r[i-1]/(r[i] + eta)
    }
  }
  age_days_midpoint <- c((age_days[-n_age] + age_days[-1])/2, age_days[n_age])
  age20 <- which.min(abs(age_days_midpoint - (20*365)))
  psi <- 1 - rho*exp(-age_days_midpoint / a0)
  
  out <- list(prop = prop,
              r = r,
              age_days_midpoint = age_days_midpoint,
              psi = psi,
              age20 = age20)
  return(out)
}



