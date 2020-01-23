#------------------------------------------------
#' @title Simulate data at a single site
#'
#' @description Simulate data in the same format as real world data by drawing
#'   from the likelihood framework. This function simulates data at a single
#'   study and single site (arbitrarily given a value of 1), see
#'   \code{?sim_study()} and \code{?sim_complete()} for alternatives.
#'
#' @param p a parameters object of class "mmfit_params" (see
#'   \code{?load_parameter_set()}).
#' @param age0,age1 lower and upper bounds on observed age ranges.
#' @param n when simulating prevalence, the number of individuals tested in each
#'   age group.
#' @param clin_p_years when simulating incidence, the time at risk of
#'   individuals in each age group.
#' @param EIR the entomological inoculation rate of the site, in units of
#'   infectious bites per person per year.
#' @param case_detection_method code for case detection method: \itemize{ \item
#'   1 Daily active case detection (ACD) \item 2 Weekly ACD \item 3 Passive case
#'   detection (PCD) }
#' @param treatment treatment rate in population.
#' @param sigma_p standard deviation of the normally-distributed study-level
#'   random effect on prevalence.
#' @param theta overdispersion parameter of the beta-binomial distribution on
#'   prevalence. 0 = standard binomial, and theta must be <1.
#' @param sigma_c standard deviation of the normally-distributed study-level
#'   random effect on incidence.
#' @param alpha_c standard deviation of the gamma-distributed site- and
#'   age-level random effect on incidence.
#' @param slide_on whether to report results of slide prevalence.
#' @param PCR_on whether to report results of PCR prevalence.
#' @param clin_on whether to report results of clinical incidence.
#'
#' @importFrom  stats rbinom rpois
#' @export

sim_site <- function(p,
                     age0 = c(0,1,10),
                     age1 = c(1,5,20),
                     n = c(10,10,10),
                     clin_p_years = c(10,10,10),
                     EIR = 5,
                     case_detection_method = 1,
                     treatment = 0.2,
                     sigma_p = 0,
                     theta = 0.1,
                     sigma_c = 0,
                     alpha_c = 1,
                     slide_on = TRUE,
                     PCR_on = TRUE,
                     clin_on = TRUE) {
  
  # check inputs
  assert_custom_class(p, "mmfit_params")
  assert_vector_pos(age0)
  assert_vector_pos(age1)
  assert_same_length(age0, age1)
  assert_gr(age1, age0)
  assert_vector(n)
  assert_pos_int(n, zero_allowed = FALSE)
  assert_same_length(n, age0)
  assert_vector(clin_p_years)
  assert_pos(clin_p_years, zero_allowed = FALSE)
  assert_same_length(clin_p_years, age0)
  assert_single_pos(EIR)
  assert_single_pos_int(case_detection_method)
  assert_in(case_detection_method, 1:3)
  assert_single_bounded(treatment)
  assert_single_pos(sigma_p, zero_allowed = TRUE)
  assert_single_pos(sigma_c, zero_allowed = TRUE)
  assert_single_pos(alpha_c, zero_allowed = TRUE)
  assert_single_bounded(theta, inclusive_right = FALSE)
  assert_single_logical(slide_on)
  assert_single_logical(PCR_on)
  assert_single_logical(clin_on)
  
  # get total number of age categories
  nc <- length(age0)
  
  # create vector of ages that contains all required age intervals as a subset
  age_vec <- sort(unique(c(age0, age1)))
  
  # get equilibrium solution for given parameters
  eq <- human_equilibrium(EIR = EIR, ft = treatment, p = p, age = age_vec)
  
  # subset to defined age intervals
  eq$states <- eq$states[match(age0, age_vec),]
  
  # simulate random effects
  w <- rnorm(1, mean = 0, sd = sigma_p)
  u <- rnorm(1, mean = 0, sd = sigma_c)
  v <- rgamma(nc, shape = 1/alpha_c, rate = 1/alpha_c)
  
  # extract raw prevalence by microscopy and PCR
  pos_M <- eq$states[,"pos_M"]
  pos_PCR <- eq$states[,"pos_PCR"]
  
  # simulate observed prevalence from beta-binomial distribution
  obs_prev <- as.data.frame(mapply(function(x) {
    
    # check for raw prevalence outside possible range
    if (any(x <= 0 | x >= 1)) {
      stop("equilibrium solution produces slide/PCR positivity outside range (0,1)")
    }
    
    # introduce random effects
    if (w == 0) {
      q <- x
    } else {
      q <- logistic(logit(x) + w)
    }
    
    # draw from beta-binomial distribution, which simplifies to binomial
    # distribution when theta == 0
    if (theta == 0) {
      prev_obs <- rbinom(nc, n, prob = x)
    } else {
      alpha <- q/theta - q
      beta <- alpha/q - alpha
      prev_obs <- rbinom(nc, size = n, prob = rbeta(nc, shape1 = alpha, shape2 = beta))
    }
    
    return(prev_obs)
  }, list(pos_M, pos_PCR)))
  names(obs_prev) <- c("slide_k", "pcr_k")
  
  # mask out slide/PCR results if not needed
  if (!slide_on) {
    obs_prev$slide_k <- NA
  }
  if (!PCR_on) {
    obs_prev$pcr_k <- NA
  }
  
  # simulate observed incidence from Poisson distribution
  obs_inc <- rep(NA, nc)
  if (clin_on) {
    inc <- eq$states[,"inc"]
    r <- c(1, 0.5, 0.2)[case_detection_method]
    lambda <- r*exp(u)*v*clin_p_years*inc
    obs_inc <- rpois(nc, lambda)
  }
  
  # create final output dataframe
  df_out <- data.frame(country = NA,
                       study = 1,
                       site = 1,
                       age0 = age0,
                       age1 = age1,
                       case_detection_method = case_detection_method,
                       n = n,
                       slide_k = obs_prev$slide_k,
                       pcr_k = obs_prev$pcr_k,
                       clin_p_years = clin_p_years,
                       clin_k = obs_inc,
                       sev_p_years = NA,
                       sev_k = NA)
  
  return(df_out)
}