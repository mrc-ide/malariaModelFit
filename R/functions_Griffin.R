
#------------------------------------------------
#' Equilibrium solution without biting heterogeneity
#'
#' Returns the equilibrium states and R0 for the model of Griffin et al. 2014 Nature Communications 5. "Estimates of the changing age-burden of Plasmodium falciparum malaria disease in sub-Saharan Africa". Does not account for biting heterogeneity - see human_equilibrium() for function that takes this into account.
#'
#' Original code due to Jamie Griffin, later modified by Xiaoyu Liu and Bob Verity.
#'
#' @param EIR EIR for adults
#' @param ft proportion of clinical cases effectively treated
#' @param p vector of model parameters
#' @param age vector of age groups
#'
#' @export

human_equilibrium_noHet <- function(EIR, ft, p, age) {
    
	na <- length(age)
	age_days <- age*365
	EIR <- EIR/365
    
	ICA 	<-	vector(length=na, mode="numeric")
	ICM 	<-	vector(length=na, mode="numeric")
	
	S   	<-	vector(length=na, mode="numeric")
	A   	<-	vector(length=na, mode="numeric")
	T   	<-	vector(length=na, mode="numeric")
	D   	<-	vector(length=na, mode="numeric")
	U   	<-	vector(length=na, mode="numeric")
	P		<-	vector(length=na, mode="numeric")
	
	prop	<-  vector(length=na, mode="numeric")
	r		<-  vector(length=na, mode="numeric")
	FOI		<-  vector(length=na, mode="numeric")
	phi		<-  vector(length=na, mode="numeric")
	q		<-  vector(length=na, mode="numeric")
	FOIM	<-  vector(length=na, mode="numeric")
	cA		<-  vector(length=na, mode="numeric")
	pos_M		<-  vector(length=na, mode="numeric")
	pos_PCR	<-  vector(length=na, mode="numeric")
	inc  <-  vector(length=na, mode="numeric")
    
    # produce population age distribution using eta, which is defined as 1/average age. The population distribution can be envisioned as an exponential distribution with mass feeding in from the left due to birth at rate eta, and leaking out of all categories with rate eta (making it stable). Rates of ageing are defined based on the width of the age groups.
	for(i in 1:na){
        r[i] <- if(i==na) 0 else 1/(age_days[i+1]-age_days[i])  # r[i] can be thought of as the rate of ageing in this age group, i.e. 1/r[i] is the duration of this group
        prop[i] <- (if(i==1) p$eta else r[i-1]*prop[i-1])/(r[i]+p$eta)  # prop is calculated as the relative flows into vs. out of the age group. For the first group the flow in is the birth rate (eta) and the flow out is the age rate. For all subsequent groups the flow in the ageing from the previous group, and the flow out is the age rate.
        if(i<na) age_days[i] <- 0.5*(age_days[i]+age_days[i+1]) # convert age_days to midpoint of range
	}
    
    psi <- 1-p$rho*exp(-age_days/p$a0)  # relative biting rate for each age group
	IB <- 0
	IC <- 0
	ID <- 0
    age20 <- which.min(abs(age_days-(20*365)))  # age category that represents 20 year old women
	for(i in 1:na){
        
        re <- r[i] + p$eta  # flow out of this category, i.e. rate of ageing plus death
        
		IB <- (EIR*psi[i]/(EIR*psi[i]*p$ub+1) + re*IB)/(1/p$db+re)
        
        # calculate probability of infection from pre-erythrocytic immunity IB via Hill function
        b <- p$b0*(p$b1+(1-p$b1)/(1+(IB/p$IB0)^p$kb))
        
        # calculate force of infection
		FOI[i] <- EIR*b*psi[i]
		
        IC <- (FOI[i]/(FOI[i]*p$uc+1) + re*IC)/(1/p$dc+re)  # denom = rate immune decay + rate ageing + eta
		ICA[i] <- IC
        
		ID <- (FOI[i]/(FOI[i]*p$ud+1) + re*ID)/(1/p$dd+re)
        q[i] <- p$d1+(1-p$d1)/(1+(ID/p$ID0)^p$kd*(1-(1-p$fd0)/(1+(age_days[i]/p$ad0)^p$gd)))
		cA[i] <- p$cU + (p$cD-p$cU)*q[i]^p$g_inf
	}
    
	IM0 <- ICA[age20]*p$PM
	ICM[1] <- IM0*(r[1] + p$eta)/(1/p$dm+r[1] + p$eta)
	for(i in 2:na) ICM[i] <- (r[i] + p$eta)*ICM[i-1]/(1/p$dm+r[i] + p$eta)
		
	phi <- p$phi0*(p$phi1+(1-p$phi1)/(1+((ICA+ICM)/p$IC0)^p$kc))
    
	for(i in 1:na){
		re <- r[i] + p$eta
		betaS <- FOI[i] + re
		betaT <- p$rT + re
		betaD <- p$rD + re
		betaA <- phi[i]*FOI[i]+ p$rA + re
		betaU <- FOI[i]+ p$rU + re
		betaP <- p$rP + re
		
		aT <- ft*phi[i]*FOI[i]/betaT
		bT <- if(i==1) 0 else r[i-1]*T[i-1]/betaT
		aD <- (1-ft)*phi[i]*FOI[i]/betaD
		bD <- if(i==1) 0 else r[i-1]*D[i-1]/betaD
		aP <- p$rT*aT/betaP
		bP <- (p$rT*bT + (if(i==1) 0 else r[i-1]*P[i-1]))/betaP
	
		Y <- (prop[i]-(bT+bD+bP))/(1+aT+aD+aP)
		
		T[i] <- aT*Y+bT
		D[i] <- aD*Y+bD
		P[i] <- aP*Y+bP
		
		A[i] <- ((if(i==1) 0 else r[i-1]*A[i-1]) + (1-phi[i])*FOI[i]*Y + p$rD*D[i])/(betaA+(1-phi[i])*FOI[i])
		U[i] <- (p$rA*A[i]+(if(i==1) 0 else r[i-1]*U[i-1]))/betaU
		S[i] <- Y-A[i]-U[i]
		
		pos_M[i] <- D[i] + T[i] + A[i]*q[i] # Microsopy
		pos_PCR[i] <- D[i] + T[i] + A[i]*(q[i]^p$aA) + U[i]*(q[i]^p$aU) # PCR
		inc[i] <- phi[i]*FOI[i]*Y
	}
	
	inf <- p$cD*D + p$cT*T + cA*A + p$cU*U

	return(cbind(age=age, S=S, T=T, D=D, A=A, U=U, P=P, inf=inf, prop=prop, psi=psi, pos_M=pos_M, pos_PCR=pos_PCR, inc=inc))
}

#------------------------------------------------
#' Compiled human_equilibrium_no_het
#'
#' Compiled version of human_equilibrium_no_het() function - see original function for details.
#'
#' @export

human_equilibrium_noHet_compiled <- cmpfun(human_equilibrium_noHet, options=list(optimize=3))


#------------------------------------------------
#' Equilibrium solution
#'
#' Returns the equilibrium states and R0 for the model of Griffin et al. 2014 Nature Communications 5. "Estimates of the changing age-burden of Plasmodium falciparum malaria disease in sub-Saharan Africa". Integrates over distribution of biting heterogeneity using Gaussian quadrature.
#'
#' Original code due to Jamie Griffin, later modified by Xiaoyu Liu and Bob Verity.
#'
#' @param EIR EIR for adults
#' @param ft proportion of clinical cases effectively treated
#' @param p vector of model parameters
#' @param age vector of age groups
#' @param h Gauss-Hermite nodes and weights for heterogeneity in biting
#'
#' @export

human_equilibrium <- function(EIR, ft, p, age, h) {
    
    nh <- length(h$nodes)
    E <- matrix() 	# states by age, mean over heterogeneity groups, not weighted by onward biting rate
    FOIM <- 0 		# overall force of infection on mosquitoes, weighted by onward biting rates
    
    # loop through all Gaussian quadrature nodes
    for(j in 1:nh){
        zeta <- exp(-p$s2*0.5 + sqrt(p$s2)*h$nodes[j])
        Ej <- human_equilibrium_noHet_compiled(EIR=EIR*zeta, ft=ft, p=p, age=age)
        E <- Ej*h$weights[j] + (if(j==1) 0 else E)
        FOIM <- FOIM + sum(Ej[,"inf"]*Ej[,"psi"])*h$weights[j]*zeta
    }
    omega <- 1-p$rho*p$eta/(p$eta+1/p$a0)
    alpha <- p$f*p$Q0
    return(list(states=E, FOIM=FOIM*alpha/omega))
}

#------------------------------------------------
#' find_prev
#'
#' Find prevalence by microscopy or PCR in a given age group. Uses pre-computed equilibrium solution as input.
#'
#' Original code due to Jamie Griffin, later modified by Xiaoyu Liu and Bob Verity.
#'
#' @param eq equilibrium solution, as returned by human_equilibrium()
#' @param age0 lower bound on age range over which to calculate prevalence
#' @param age1 upper bound on age range over which to calculate prevalence
#' @param prevType whether to calculate prevalence by "microscopy" or "PCR"
#'
#' @export

find_prev <- function(eq, age0, age1, prevType="microscopy") {
    
    # check input formats
    stopifnot(prevType%in%c("microscopy","PCR"))
    stopifnot(is.list(eq) && !is.data.frame(eq))
    stopifnot(all(names(eq)==c("states","FOIM")))
    
    # get weights (proportion of each age bin covered by age range)
    w <-  breakCoverage(eq$states[,"age"], age0, age1)
    w <- c(w,0)
    
    # get prevalence in this age range
    typeName <- switch(prevType, microscopy="pos_M", PCR="pos_PCR")
    pos <- eq$states[,typeName]
    prop <- eq$states[,"prop"]
    prev <- sum(w*pos)/sum(w*prop)
    
    return(prev)
}

#------------------------------------------------
#' find_inc
#'
#' Find incidence in a given age group. Uses pre-computed equilibrium solution as input.
#'
#' Original code due to Jamie Griffin, later modified by Xiaoyu Liu and Bob Verity.
#'
#' @param eq equilibrium solution, as returned by human_equilibrium()
#' @param age0 lower bound on age range over which to calculate prevalence
#' @param age1 upper bound on age range over which to calculate prevalence
#'
#' @export

find_inc <- function(eq, age0, age1) {
    
    # check input formats
    stopifnot(is.list(eq) && !is.data.frame(eq))
    stopifnot(all(names(eq)==c("states","FOIM")))
    
    # get weights (proportion of each age bin covered by age range)
    w <- breakCoverage(eq$states[,"age"], age0, age1)
    w <- c(w,0)
    
    # get clinical incidence in this age range
    pos <- eq$states[,"inc"]
    prop <- eq$states[,"prop"]
    inc <- sum(w*pos)/sum(w*prop)
    
    return(inc)
}

#------------------------------------------------

find_R0_simple <- function(eq, EIR, ft, p, age, h){
    
    # R0 that doesn't account for ageing during an infection
    omega <- 1-p$rho*p$eta/(p$eta+1/p$a0)
    alpha <- p$f*p$Q0
    C <- p$phi0*ft*p$cT/(p$rT+p$eta) + p$phi0*(1-ft)*p$cD/(p$rD+p$eta) + (1-p$phi0 + p$phi0*(1-ft)*p$rD/(p$rD+p$eta))*(p$cD/(p$rA+p$eta)+p$rA/(p$rA+p$eta)*p$cU/(p$rU+p$eta))
    da <- 1/omega^2*(1-2*p$rho/(1+1/(p$a0*p$eta))+p$rho^2/(1+2/(p$a0*p$eta)))
    dh <- exp(p$s2)
    
    # multiply by omega to go from mean EIR for adults to mean EIR in population
    EIRd <- EIR/365*omega
    gamma <- exp(-p$tau*p$mu)/p$mu
    
    m <- EIRd*(eq$FOIM+p$mu)/(alpha*eq$FOIM*exp(-p$tau*p$mu))
    R0 <- m*gamma*alpha^2*da*dh*p$b0*C
    return(R0)
}
