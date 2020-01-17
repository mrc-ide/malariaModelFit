

age_prep <- function(age, eta, rho, a0){
  # get basic properties
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
  
  out <- list(na = n_age,
            age = age_days,
            prop = prop,
            r = r,
            age_days_midpoint = age_days_midpoint,
            psi = psi,
            age20 = age20 - 1)
  return(out)
}

data_prep <- function(output_type, study_n, EIR, ft,
                      age, eta, rho, a0, nh){
  age_d <- unlist(age_prep(age, eta, rho, a0))
  het_d <- c(nh, unlist(gq_normal(nh)))
  
  c(output_type, study_n, EIR, ft, age_d, het_d)
}