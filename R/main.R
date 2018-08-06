#------------------------------------------------
# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

#' @useDynLib malariaModelFit
#' @import compiler
#' @import statmod
#' @import ggplot2
#' @importFrom Rcpp evalCpp
#' @import graphics
#' @import stats
#' @import utils
NULL

#------------------------------------------------
#' @title Dummy function
#'
#' @description Simple dummy function to check that R and Rcpp are talking to each other.
#'
#' @export

dummy1 <- function() {
  cat("R function working!\n")
  dummy1_cpp()
}

#------------------------------------------------
#' @title Equilibrium solution without biting heterogeneity
#' 
#' @description Returns the equilibrium states for the model of Griffin et al.
#'   2014 Nature Communications 5. "Estimates of the changing age-burden of
#'   Plasmodium falciparum malaria disease in sub-Saharan Africa". For
#'   derivation of equilibrium solutions see Griffin 2016 Malaria Journal "Is a
#'   reproduction number of one a threshold for Plasmodium falciparum malaria
#'   elimination?" supplementary material.
#'   
#'   This function does not account for biting heterogeneity - see 
#'   human_equilibrium() for function that takes this into account.
#'   
#'   Original code due to Jamie Griffin, later modified by Xiaoyu Liu and Bob 
#'   Verity.
#'
#' @param EIR EIR for adults
#' @param ft proportion of clinical cases effectively treated
#' @param p vector of model parameters
#' @param age vector of age groups
#'
#' @export

human_equilibrium_no_het <- function(EIR, ft, p, age) {
  
  # check inputs
  assert_scalar_pos(EIR)
  assert_scalar_pos(ft)
  assert_bounded(ft)
  check_names(names(p), required_eq_params())
  assert_pos(age)
  
  # get basic properties
  n_age <- length(age)
  age_days <- age*365
  EIR <- EIR/365
  
  # produce population age distribution using eta, which is defined as 1/average
  # age. The population distribution can be envisioned as an exponential 
  # distribution, with mass feeding in from the left due to birth at rate eta, 
  # people ageing with rates defined based on the width of the age groups, and 
  # mass leaking out of all categories with death rate eta, making the 
  # distribution stable.
  prop <- r <- rep(0,n_age)
	for (i in 1:n_age) {
	  
    # r[i] can be thought of as the rate of ageing in this age group, i.e.
    # 1/r[i] is the duration of this group
	  if (i==n_age) {
	    r[i] <- 0
	  } else {
	    age_width <- age_days[i+1] - age_days[i]
	    r[i] <- 1/age_width
	  }
    
    # prop is calculated as the relative flows into vs. out of this age group. 
    # For the first group the flow in is the birth rate (eta) and the flow out 
    # is the age rate. For all subsequent groups the flow in represents ageing
    # from the previous group, and the flow out is the age rate.
	  if (i==1) {
	    prop[i] <- p$eta/(r[i]+p$eta)
	  } else {
	    prop[i] <- prop[i-1]*r[i-1]/(r[i]+p$eta)
	  }
	}
  
  # calculate midpoint of age range
  age_days_midpoint <- c((age_days[-n_age]+age_days[-1])/2, age_days[n_age])
  
  # get age category that represents a 20 year old woman
  age20 <- which.min(abs(age_days_midpoint-(20*365)))
  
  # relative biting rate for each age group
  psi <- 1 - p$rho*exp(-age_days_midpoint/p$a0)
  
  # calculate immunity functions and onward infectiousness at equilibrium for
  # all age groups. See references above for details of derivation.
  IB <- IC <- ID <- 0
  ICA <- FOI <- q <- cA <- rep(0,n_age)
  for (i in 1:n_age) {
    
    # rate of ageing plus death
    re <- r[i] + p$eta
    
    # update pre-erythrocytic immunity IB
    eps <- EIR*psi[i]
    IB <- (eps/(eps*p$ub+1) + re*IB)/(1/p$db + re)
    
    # calculate probability of infection from pre-erythrocytic immunity IB via
    # Hill function
    b <- p$b0*(p$b1 + (1-p$b1)/(1+(IB/p$IB0)^p$kb))
    
    # calculate force of infection (lambda)
    FOI[i] <- b*eps
    
    # update clinical immunity IC
    IC <- (FOI[i]/(FOI[i]*p$uc+1) + re*IC)/(1/p$dc + re)
    ICA[i] <- IC
    
    # update detection immunity ID
    ID <- (FOI[i]/(FOI[i]*p$ud+1) + re*ID)/(1/p$dd + re)
    
    # calculate probability that an asymptomatic infection (state A) will be
    # detected by microscopy
    fd <- 1 - (1-p$fd0)/(1+(age_days_midpoint[i]/p$ad0)^p$gd)
    q[i] <- p$d1 + (1-p$d1)/(1+(ID/p$ID0)^p$kd*fd)
    
    # calculate onward infectiousness to mosquitoes
    cA[i] <- p$cU + (p$cD-p$cU)*q[i]^p$g_inf
  }
  
  # calculate maternal clinical immunity, assumed to be at birth a proportion of
  # the acquired immunity of a 20 year old
  IM0 <- ICA[age20]*p$PM
  ICM <- rep(0,n_age)
  for (i in 1:n_age) {
    
    # rate of ageing plus death
    re <- r[i] + p$eta
    
    # maternal clinical immunity decays from birth
    ICM_prev <- ifelse(i==1, IM0, ICM[i-1])
    ICM[i] <- ICM_prev*re/(1/p$dm + re)
  }
  
  # calculate probability of acquiring clinical disease as a function of 
  # different immunity types
	phi <- p$phi0*(p$phi1 + (1-p$phi1)/(1+((ICA+ICM)/p$IC0)^p$kc))
	
	# calculate equilibrium solution of all model states. Again, see references
	# above for details
  pos_M <- pos_PCR <- inc <- rep(0,n_age)
	S <- T <- P <- D <- A <- U <- rep(0,n_age)
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
    bT <- if(i==1) 0 else r[i-1]*T[i-1]/betaT
    aD <- (1-ft)*phi[i]*FOI[i]/betaD
    bD <- if(i==1) 0 else r[i-1]*D[i-1]/betaD
    aP <- p$rT*aT/betaP
    bP <- (p$rT*bT + (if(i==1) 0 else r[i-1]*P[i-1]))/betaP
    
    # calculate Y
    Y <- (prop[i] - (bT+bD+bP))/(1+aT+aD+aP)
    
    # calculate final {T,D,P} solution
    T[i] <- aT*Y+bT
    D[i] <- aD*Y+bD
    P[i] <- aP*Y+bP
    
    # calculate final {A, U, S} solution
    if (i==1) {
      rA <- rU <- 0
    } else {
      rA <- r[i-1]*A[i-1]
      rU <- r[i-1]*U[i-1]
    }
    A[i] <- (rA + (1-phi[i])*Y*FOI[i] + p$rD*D[i])/(betaA + (1-phi[i])*FOI[i])
    U[i] <- (rU + p$rA*A[i])/betaU
    S[i] <- Y-A[i]-U[i]
    
    # calculate proportion detectable by mocroscopy and PCR
    pos_M[i] <- D[i] + T[i] + A[i]*q[i] # microsopy
    pos_PCR[i] <- D[i] + T[i] + A[i]*(q[i]^p$aA) + U[i]*(q[i]^p$aU) # PCR
    
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
#'   2014 Nature Communications 5. "Estimates of the changing age-burden of
#'   Plasmodium falciparum malaria disease in sub-Saharan Africa". For
#'   derivation of equilibrium solutions see Griffin 2016 Malaria Journal "Is a
#'   reproduction number of one a threshold for Plasmodium falciparum malaria
#'   elimination?" supplementary material.
#'   
#'   Integrates over distribution of biting heterogeneity using Gaussian 
#'   quadrature.
#'   
#'   Original code due to Jamie Griffin, later modified by Xiaoyu Liu and Bob 
#'   Verity.
#'
#' @param EIR EIR for adults
#' @param ft proportion of clinical cases effectively treated
#' @param p vector of model parameters
#' @param age vector of age groups
#' @param h Gauss-Hermite nodes and weights for heterogeneity in biting
#'
#' @export

human_equilibrium <- function(EIR, ft, p, age, h) {
  
  # check inputs
  assert_scalar_pos(EIR)
  assert_scalar_pos(ft)
  assert_bounded(ft)
  check_names(names(p), required_eq_params())
  assert_pos(age)
  assert_in(c("nodes", "weights"), names(h))
  
  # get basic properties and initialise
  nh <- length(h$nodes)
  E <- matrix() 	# states by age, taking the mean over heterogeneity groups
  FOIM <- 0 		# overall force of infection on mosquitoes, weighted by onward biting rates
  
  # loop through all Gaussian quadrature nodes
  for (j in 1:nh) {
    zeta <- exp(-p$s2*0.5 + sqrt(p$s2)*h$nodes[j])
    Ej <- human_equilibrium_no_het(EIR=EIR*zeta, ft=ft, p=p, age=age)
    E <- Ej*h$weights[j] + (if(j==1) 0 else E)
    FOIM <- FOIM + sum(Ej[,"inf"]*Ej[,"psi"])*h$weights[j]*zeta
  }
  
  # complete overall force of infection on mosquitoes
  omega <- 1 - p$rho*p$eta/(p$eta+1/p$a0)
  alpha <- p$f*p$Q0
  FOIM <- FOIM*alpha/omega
  
  # return as list
  return(list(states = E,
              FOIM = FOIM))
}

#------------------------------------------------
#' @title find_prev
#'
#' @description Find prevalence by microscopy or PCR in a given age group, using
#'   pre-computed equilibrium solution as input.
#' 
#' Original code due to Jamie Griffin, later modified by Xiaoyu Liu and Bob 
#' Verity.
#'
#' @param eq equilibrium solution, as returned by \code{human_equilibrium()}
#' @param age0 lower bound on age range over which to calculate prevalence
#' @param age1 upper bound on age range over which to calculate prevalence
#' @param prev_type whether to calculate prevalence by "microscopy" or "PCR"
#'
#' @export

find_prev <- function(eq, age0, age1, prev_type="microscopy") {
  
  # check input formats
  assert_in("states", names(eq))
  assert_pos(age0)
  assert_pos(age1)
  assert_length(prev_type, 1)
  assert_in(prev_type, c("microscopy", "PCR"))
  if (prev_type=="microscopy") {
    assert_in(c("age", "pos_M", "prop"), names(eq$states))
  } else {
    assert_in(c("age", "pos_PCR", "prop"), names(eq$states))
  }
  
  # get weights (proportion of each age bin covered by age range)
  w <-  break_coverage(eq$states[,"age"], age0, age1)
  w <- c(w,0)
  
  # get prevalence in this age range
  type_name <- switch(prev_type, microscopy="pos_M", PCR="pos_PCR")
  pos <- eq$states[,type_name]
  prop <- eq$states[,"prop"]
  prev <- sum(w*pos)/sum(w*prop)
  
  return(prev)
}

#------------------------------------------------
#' @title find_inc
#' 
#' @description Find incidence in a given age group, using pre-computed
#'   equilibrium solution as input.
#' 
#' Original code due to Jamie Griffin, later modified by Xiaoyu Liu and Bob 
#' Verity.
#'
#' @param eq equilibrium solution, as returned by \code{human_equilibrium()}
#' @param age0 lower bound on age range over which to calculate incidence
#' @param age1 upper bound on age range over which to calculate incidence
#'
#' @export

find_inc <- function(eq, age0, age1) {
  
  # check input formats
  assert_in("states", names(eq))
  assert_in(c("age", "inc", "prop"), names(eq$states))
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
