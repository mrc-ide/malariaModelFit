#------------------------------------------------
#' @title Laod a Parameter Set from File
#'
#' @description Parameter sets are stored within the package
#'   inst/extdata/parameter_sets folder. Load one of these sets by name.
#'
#' @param file_name the name of a parameter set within the
#'   inst/extdata/parameter_sets folder.
#'
#' @export

load_parameter_set <- function(file_name = "Jamie_parameters.rds") {
  
  # check inputs
  assert_single_string(file_name)
  
  # load parameter set from inst/extdata/parameter_sets folder
  params <- mmfit_file(paste0("parameter_sets/", file_name))
  
  return(params)
}

#------------------------------------------------
#' @title Equilibrium solution without biting heterogeneity
#' 
#' @description Returns the equilibrium states for the model of Griffin et al.
#'   (2014). A derivation of the equilibrium solutions can be found in Griffin
#'   (2016).
#'   
#'   This function does not account for biting heterogeneity - see
#'   \code{human_equilibrium()} for function that takes this into account.
#'
#' @param EIR EIR for adults, in units of infectious bites per person per year
#' @param ft proportion of clinical cases effectively treated
#' @param p vector of model parameters
#' @param age vector of age groups, in units of years
#'
#' @references Griffin et. al. (2014). Estimates of the changing age-burden of
#'   Plasmodium falciparum malaria disease in sub-Saharan Africa.
#'   doi:10.1038/ncomms4136
#'   
#'   Griffin (2016). Is a reproduction number of one a threshold for Plasmodium
#'   falciparum malaria elimination? doi:10.1186/s12936-016-1437-9 (see
#'   supplementary material)
#'
#' @export

human_equilibrium_no_het <- function(EIR, ft, p, age) {
  
  # check inputs
  assert_single_pos(EIR, zero_allowed = FALSE)
  assert_single_bounded(ft)
  assert_custom_class(p, "mmfit_params")
  assert_vector_pos(age)
  assert_noduplicates(age)
  assert_increasing(age)
  
  # get basic properties
  n_age <- length(age)
  age_days <- age*365
  EIR <- EIR/365
  
  # produce population age distribution using eta, which is defined as 1/average
  # age. The population distribution can be envisioned as an exponential
  # distribution, with mass feeding in from the left due to birth at rate eta,
  # people ageing with rates defined based on the width of the age groups, and
  # mass leaking out of all categories with death rate eta. Total birth and
  # death rates are equal, making the distribution stable.
  prop <- r <- rep(0, n_age)
	for (i in 1:n_age) {
	  
    # r[i] can be thought of as the rate of ageing in this age group, i.e.
    # 1/r[i] is the duration of this group
	  if (i == n_age) {
	    r[i] <- 0
	  } else {
	    age_width <- age_days[i+1] - age_days[i]
	    r[i] <- 1/age_width
	  }
    
    # prop is calculated from the relative flows into vs. out of this age group.
    # For the first age group the flow in is equal to the birth rate (eta). For
    # all subsequent age groups the flow in represents ageing from the previous
    # group. The flow out is always equal to the rate of ageing or death.
	  if (i == 1) {
	    prop[i] <- p$eta/(r[i] + p$eta)
	  } else {
	    prop[i] <- prop[i-1]*r[i-1]/(r[i] + p$eta)
	  }
	}
  
  # calculate midpoint of age range. There is no midpoint for the final age
  # group, so use beginning of range instead
  age_days_midpoint <- c((age_days[-n_age] + age_days[-1])/2, age_days[n_age])
  
  # get age category that represents a 20 year old woman
  age20 <- which.min(abs(age_days_midpoint - (20*365)))
  
  # calculate relative biting rate for each age group
  psi <- 1 - p$rho*exp(-age_days_midpoint/p$a0)
  
  # calculate immunity functions and onward infectiousness at equilibrium for
  # all age groups. See doi:10.1186/s12936-016-1437-9 for details of derivation.
  IB <- IC <- ID <- 0
  ICA <- FOI <- q <- cA <- rep(0, n_age)
  for (i in 1:n_age) {
    
    # rate of ageing plus death
    re <- r[i] + p$eta
    
    # update pre-erythrocytic immunity IB
    eps <- EIR*psi[i]
    IB <- (eps/(eps*p$ub + 1) + re*IB)/(1/p$db + re)
    
    # calculate probability of infection from pre-erythrocytic immunity IB via
    # Hill function
    b <- p$b0*(p$b1 + (1-p$b1)/(1+(IB/p$IB0)^p$kb))
    
    # calculate force of infection (lambda)
    FOI[i] <- b*eps
    
    # update clinical immunity IC
    IC <- (FOI[i]/(FOI[i]*p$uc + 1) + re*IC)/(1/p$dc + re)
    ICA[i] <- IC
    
    # update detection immunity ID
    ID <- (FOI[i]/(FOI[i]*p$ud + 1) + re*ID)/(1/p$dd + re)
    
    # calculate probability that an asymptomatic infection (state A) will be
    # detected by microscopy
    fd <- 1 - (1-p$fd0)/(1 + (age_days_midpoint[i]/p$ad0)^p$gd)
    q[i] <- p$d1 + (1-p$d1)/(1 + (ID/p$ID0)^p$kd*fd)
    
    # calculate onward infectiousness to mosquitoes
    cA[i] <- p$cU + (p$cD-p$cU)*q[i]^p$g_inf
  }
  
  # calculate maternal clinical immunity, assumed to be at birth a proportion of
  # the acquired immunity of a 20 year old
  IM0 <- ICA[age20]*p$PM
  ICM <- rep(0, n_age)
  for (i in 1:n_age) {
    
    # rate of ageing plus death
    re <- r[i] + p$eta
    
    # maternal clinical immunity decays from birth
    if (i == 1) {
      ICM_prev <- IM0
    } else {
      ICM_prev <- ICM[i-1]
    }
    ICM[i] <- ICM_prev*re/(1/p$dm + re)
  }
  
  # calculate probability of acquiring clinical disease as a function of 
  # different immunity types
	phi <- p$phi0*(p$phi1 + (1-p$phi1)/(1 + ((ICA+ICM)/p$IC0)^p$kc))
	
	# calculate equilibrium solution of all model states. Again, see
	# doi:10.1186/s12936-016-1437-9 for details
  pos_M <- pos_PCR <- inc <- rep(0, n_age)
	S <- T <- P <- D <- A <- U <- rep(0, n_age)
  for (i in 1:n_age) {
    
    # rate of ageing plus death
    re <- r[i] + p$eta
    
    # calculate beta values
    betaT <- p$rT + re
    betaD <- p$rD + re
    betaA <- FOI[i]*phi[i] + p$rA + re
    betaU <- FOI[i] + p$rU + re
    betaP <- p$rP + re
    
    # calculate a and b values
    aT <- ft*phi[i]*FOI[i]/betaT
    aP <- p$rT*aT/betaP
    aD <- (1-ft)*phi[i]*FOI[i]/betaD
    if (i == 1) {
      bT <- bD <- bP <- 0
    } else {
      bT <- r[i-1]*T[i-1]/betaT
      bD <- r[i-1]*D[i-1]/betaD
      bP <- p$rT*bT + r[i-1]*P[i-1]/betaP
    }
    
    # calculate Y
    Y <- (prop[i] - (bT + bD + bP))/(1 + aT + aD + aP)
    
    # calculate final {T,D,P} solution
    T[i] <- aT*Y + bT
    D[i] <- aD*Y + bD
    P[i] <- aP*Y + bP
    
    # calculate final {A, U, S} solution
    if (i == 1) {
      rA <- rU <- 0
    } else {
      rA <- r[i-1]*A[i-1]
      rU <- r[i-1]*U[i-1]
    }
    A[i] <- (rA + (1-phi[i])*Y*FOI[i] + p$rD*D[i])/(betaA + (1-phi[i])*FOI[i])
    U[i] <- (rU + p$rA*A[i])/betaU
    S[i] <- Y - A[i] - U[i]
    
    # calculate proportion detectable by mocroscopy and PCR
    pos_M[i] <- D[i] + T[i] + A[i]*q[i]
    pos_PCR[i] <- D[i] + T[i] + A[i]*(q[i]^p$aA) + U[i]*(q[i]^p$aU)
    
    # calculate clinical incidence
    inc[i] <- Y*FOI[i]*phi[i]
  }
	
	# calculate incidence of infection
	inf <- p$cD*D + p$cT*T + cA*A + p$cU*U
  
	# return matrix
	ret <- cbind(age = age,
	             S = S,
	             T = T,
	             D = D,
	             A = A,
	             U = U,
	             P = P,
	             inf = inf,
	             prop = prop,
	             psi = psi,
	             pos_M = pos_M,
	             pos_PCR = pos_PCR,
	             inc = inc)
	return(ret)
}

#------------------------------------------------
#' @title Equilibrium solution
#'
#' @description Returns the equilibrium states for the model of Griffin et al.
#'   (2014). A derivation of the equilibrium solutions can be found in Griffin
#'   (2016). Integrates over the distribution of biting heterogeneity using
#'   Gaussian quadrature.
#' 
#' @inheritParams human_equilibrium_no_het
#' @param h a list of Gauss-Hermite nodes and associated weights, used for
#'   integrating over heterogeneity in biting. See \code{?gq_normal} for an
#'   example.
#'
#' @references Griffin et. al. (2014). Estimates of the changing age-burden of
#'   Plasmodium falciparum malaria disease in sub-Saharan Africa.
#'   doi:10.1038/ncomms4136
#'   
#'   Griffin (2016). Is a reproduction number of one a threshold for Plasmodium
#'   falciparum malaria elimination? doi:10.1186/s12936-016-1437-9 (see
#'   supplementary material)
#'
#' @export

human_equilibrium <- function(EIR, ft, p, age, h = gq_normal(10)) {
  
  # check inputs
  assert_single_pos(EIR, zero_allowed = FALSE)
  assert_single_bounded(ft)
  assert_custom_class(p, "mmfit_params")
  assert_vector_pos(age)
  assert_noduplicates(age)
  assert_increasing(age)
  assert_list(h)
  assert_in(c("nodes", "weights"), names(h))
  assert_same_length(h$nodes, h$weights)
  assert_numeric(h$nodes, h$weights)
  
  # get basic properties and initialise
  nh <- length(h$nodes)
  FOIM <- 0 		# overall force of infection on mosquitoes, weighted by onward biting rates
  
  # loop through all Gaussian quadrature nodes
  for (j in 1:nh) {
    zeta <- exp(-p$s2*0.5 + sqrt(p$s2)*h$nodes[j])
    Ej <- human_equilibrium_no_het(EIR = EIR*zeta, ft = ft, p = p, age = age)
    if (j == 1) {
      E <- Ej*h$weights[j]
    } else {
      E <- E + Ej*h$weights[j]
    }
    FOIM <- FOIM + sum(Ej[,"inf"]*Ej[,"psi"])*h$weights[j]*zeta
  }
  
  # calculate overall force of infection on mosquitoes
  omega <- 1 - p$rho*p$eta/(p$eta + 1/p$a0)
  alpha <- p$f*p$Q0
  FOIM <- FOIM*alpha/omega
  
  # return as list
  return(list(states = E,
              FOIM = FOIM))
}

#------------------------------------------------
#' @title Get prevalence from equilibrium output
#'
#' @description Get prevalence by microscopy or PCR in a given age range, using
#'   pre-computed equilibrium solution as input.
#'
#' @param eq equilibrium solution, as returned by \code{human_equilibrium()}.
#' @param age0 lower bound on age range.
#' @param age1 upper bound on age range.
#' @param prev_type whether to calculate prevalence by "microscopy" or "PCR".
#'
#' @export

get_prev <- function(eq, age0, age1, prev_type = "microscopy") {
  
  # check input formats
  assert_in("states", names(eq))
  assert_pos(age0)
  assert_pos(age1)
  assert_length(prev_type, 1)
  assert_in(prev_type, c("microscopy", "PCR"))
  if (prev_type == "microscopy") {
    assert_in(c("age", "pos_M", "prop"), colnames(eq$states))
  } else {
    assert_in(c("age", "pos_PCR", "prop"), colnames(eq$states))
  }
  
  # get weights (proportion of each age bin covered by age range)
  w <-  break_coverage(eq$states[,"age"], age0, age1)
  w <- c(w,0)
  
  # get prevalence in this age range
  type_name <- switch(prev_type,
                      microscopy = "pos_M",
                      PCR = "pos_PCR")
  pos <- eq$states[,type_name]
  prop <- eq$states[,"prop"]
  prev <- sum(w*pos)/sum(w*prop)
  
  return(prev)
}

#------------------------------------------------
#' @title Get incidence from equilibrium output
#' 
#' @description Get incidence in a given age range, using pre-computed
#'   equilibrium solution as input.
#'
#' @param eq equilibrium solution, as returned by \code{human_equilibrium()}.
#' @param age0 lower bound on age range.
#' @param age1 upper bound on age range.
#'
#' @export

get_inc <- function(eq, age0, age1) {
  
  # check input formats
  assert_in("states", names(eq))
  assert_in(c("age", "inc", "prop"), colnames(eq$states))
  assert_pos(age0)
  assert_pos(age1)
  
  # get weights (proportion of each age bin covered by age range)
  w <- break_coverage(eq$states[,"age"], age0, age1)
  w <- c(w,0)
  
  # get clinical incidence in this age range
  pos <- eq$states[,"inc"]
  prop <- eq$states[,"prop"]
  inc <- sum(w*pos)/sum(w*prop)
  
  return(inc)
}
