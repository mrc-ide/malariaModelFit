

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
            age20 = age20)
  return(out)
}

data_prep <- function(data, eta, rho, a0, output_type = "ll", nh = 5){
  ot <- 0
  if(output_type == "ll"){
    ot <- 1
  }
  # Split data by site
  site_list <- split(data, data$site_index)
  # Site N and size
  site_n <- length(site_list)
  group_n <- sapply(site_list, nrow) + 1
  names(group_n) <- paste0("group_n_", 1:site_n)
  # Age inputs
  ages <- lapply(site_list, function(x){
    age_prep(c(x$age0[1], x$age1), eta = eta, rho = rho, a0 = a0)
  })
  prop <- unlist(lapply(ages, "[[", "prop"))
  names(prop) <- paste0("prop", 1:length(prop))
  r <- unlist(lapply(ages, "[[", "r"))
  names(r) <- paste0("r", 1:length(r))
  age_days_midpoint <- unlist(lapply(ages, "[[", "age_days_midpoint"))
  names(age_days_midpoint) <- paste0("age_days_midpoint", 1:length(r))
  psi <- unlist(lapply(ages, "[[", "psi"))
  names(psi) <- paste0("psi", 1:length(psi))
  age20 <- unlist(lapply(ages, "[[", "age20"))
  names(age20) <- paste0("age20_", 1:length(age20))
  # Gaussian quandrature nodes and weights
  het_d <- c(nh, unlist(gq_normal(nh)))
  
  # Return output vector in specific format
  return(c(output_type = ot, site_n = site_n, group_n, prop, r, age_days_midpoint, psi, age20, het_d))
}

