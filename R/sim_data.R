#------------------------------------------------
#' @title Simulate data at a single site
#'
#' @description Simulate data in the same format as real world data by drawing
#'   from the likelihood framework. This function simulates data at a single
#'   study and single site (arbitrarily given a value of 1), see
#'   \code{?sim_study()} and \code{?sim_complete()} for alternatives.
#'
#' @param p a parameters object of class "model_params" (see
#'   \code{?load_parameter_set()}).
#' @param EIR the entomological inoculation rate of the site, in units of
#'   infectious bites per person per year.
#' @param treatment treatment rate in population.
#' @param type the data type to be simulated:
#' \itemize{
#'   \item 1 clinical incidence, where \code{denom} is the person-years at risk
#'   and \code{numer} becomes the number of cases in that period.
#'   \item 2 prevalence by microscopy, where \code{denom} is the number of
#'   people tested and \code{numer} becomes the number of positive tests.
#' }
#' @param age0,age1 lower and upper bounds on observed age ranges.
#' @param denom the denominator used in simulation. When simulating prevalence
#'   this represents the number of individuals tested in each age group, when
#'   simulating incidence this represents the time at risk of individuals in
#'   each age group.
#' @param case_detection code for case detection method:
#' \itemize{
#'   \item 1 Daily active case detection (ACD), detection probability 1.0.
#'   \item 2 Weekly ACD, detection probability 0.5.
#'   \item 3 Passive case detection (PCD), detection probability 0.2.
#' }
#' @param sigma_p standard deviation of the normally-distributed study-level
#'   random effect on prevalence.
#' @param theta overdispersion parameter of the beta-binomial distribution on
#'   prevalence. 0 = standard binomial, and theta must be <1.
#' @param sigma_c standard deviation of the normally-distributed study-level
#'   random effect on incidence.
#' @param alpha_c standard deviation of the gamma-distributed site- and
#'   age-level random effect on incidence.
#'
#' @importFrom malariaEquilibrium human_equilibrium
#' @importFrom  stats rbinom rpois
#' @export

sim_site <- function(p,
                     EIR = 5,
                     treatment = 0.2,
                     type = 1,
                     age0 = 0:4,
                     age1 = 1:5,
                     denom = rep(10,5),
                     case_detection = 1,
                     sigma_p = 0,
                     theta = 0.1,
                     sigma_c = 0,
                     alpha_c = 1) {
  
  # check inputs
  assert_custom_class(p, "model_params")
  assert_single_pos(EIR)
  assert_single_bounded(treatment)
  assert_in(type, 1:2)
  assert_vector_pos(age0)
  assert_vector_pos(age1)
  assert_same_length(age0, age1)
  assert_gr(age1, age0)
  assert_vector(denom)
  assert_pos(denom, zero_allowed = FALSE)
  assert_same_length(denom, age0)
  assert_in(case_detection, 1:3)
  assert_single_pos(sigma_p, zero_allowed = TRUE)
  assert_single_bounded(theta, inclusive_right = FALSE)
  assert_single_pos(sigma_c, zero_allowed = TRUE)
  assert_single_pos(alpha_c, zero_allowed = TRUE)
  
  # get total number of age categories
  n_age <- length(age0)
  
  # create vector of ages that contains all required age intervals as a subset
  age_vec <- sort(unique(c(age0, age1)))
  
  # get equilibrium solution for given parameters
  eq <- malariaEquilibrium::human_equilibrium(EIR = EIR, ft = treatment, p = p, age = age_vec)
  
  # subset to defined age intervals
  eq$states <- eq$states[match(age0, age_vec),]
  
  # simulate random effects
  w <- rnorm(1, mean = 0, sd = sigma_p)
  u <- rnorm(1, mean = 0, sd = sigma_c)
  v <- rgamma(n_age, shape = 1/alpha_c, rate = 1/alpha_c)
  
  # simulate incidence (type = 1) or prevalence (type = 2)
  if (type == 1) {
    
    # get daily incidence from equilibrium output
    inc_raw <- eq$states[,"inc"] * 365
    
    # define case detection probability
    r <- c(1, 0.5, 0.2)[case_detection]
    
    # draw from Poisson distribution
    lambda <- r*exp(u)*v*denom*inc_raw
    numer <- rpois(n_age, lambda)
    
  } else {
    
    # extract raw prevalence from equilibrium output
    prev_raw <- eq$states[,"pos_M"]
    
    # introduce random effects
    if (w == 0) {
      q <- prev_raw
    } else {
      q <- logistic(logit(prev_raw) + w)
    }
    
    # draw observed prevalence from beta-binomial distribution, which simplifies
    # to binomial distribution when theta == 0
    if (theta == 0) {
      numer <- rbinom(n_age, denom, prob = q)
    } else {
      alpha <- q/theta - q
      beta <- alpha/q - alpha
      numer <- rbinom(n_age, size = denom, prob = rbeta(n_age, shape1 = alpha, shape2 = beta))
    }
    
  }
  
  # create final output dataframe
  df_out <- data.frame(study_index = 1,
                       site_index = 1,
                       numer = numer,
                       denom = denom,
                       type = type,
                       age0 = age0,
                       age1 = age1,
                       case_detection = case_detection)
  
  return(df_out)
}