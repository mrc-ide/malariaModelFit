
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
    
	for(i in 1:na){
		r[i] <- if(i==na) 0 else 1/(age_days[i+1]-age_days[i])
		prop[i] <- (if(i==1) p$eta else r[i-1]*prop[i-1])/(r[i]+p$eta)
		if(i<na) age_days[i] <- 0.5*(age_days[i]+age_days[i+1])
	}
	psi <- 1-p$rho*exp(-age_days/p$a0)
	IB <- 0
	IC <- 0
	ID <- 0
	age20 <- 0
	for(i in 1:na){
		if(i<na & age_days[i]<=20*365 & age_days[i+1]>20*365) age20 <- i
		re <- r[i] + p$eta
		IB <- (EIR*psi[i]/(EIR*psi[i]*p$ub+1) + re*IB)/(1/p$db+re)
		b <- p$b0*(p$b1+(1-p$b1)/(1+(IB/p$IB0)^p$kb))
		FOI[i] <- EIR*b*psi[i]
		
		IC <- (FOI[i]/(FOI[i]*p$uc+1) + re*IC)/(1/p$dc+re)
		ICA[i] <- IC
		ID <- (FOI[i]/(FOI[i]*p$ud+1) + re*ID)/(1/p$dd+re)
		q[i] <- p$d1+(1-p$d1)/(1+(ID/p$ID0)^p$kd*(1-(1-p$fd0)/(1+(age[i]/p$ad0)^p$gd)))
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
#' Find prevalence by microscopy or PCR in a given age group. Reads in equilibrium solution as input.
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
    w <-  breakCoverageTest(eq$states[,"age"], age0, age1)
    w <- c(w,0)
    
    # get prevalence in this age range
    typeName <- switch(prevType, microscopy="pos_M", PCR="pos_PCR")
    pos <- eq$states[,typeName]
    prev <- sum(w*pos)
    
    return(prev)
}
