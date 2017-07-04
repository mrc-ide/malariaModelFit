
#------------------------------------------------
#' Default parameters
#'
#' Returns list of default parameters, taken from Griffin et al. 2014 Nature Communications 5, "Estimates of the changing age-burden of Plasmodium falciparum malaria disease in sub-Saharan Africa". All time units are in days.
#'
#' @export

default_parameters <- function(){
	
	p <- list()
	
	# age, heterogeneity in exposure
	p$eta	<-	0.0001305
	p$rho	<-	0.85
	p$a0	<-	2920
	p$s2	<-	1.67

	# rate of leaving infection states
	p$rA	<-	0.00512821
	p$rT	<-	0.2
	p$rD	<-	0.2
	p$rU	<-	0.00906627
	p$rP	<-	0.2
	
	# human latent period and time lag from asexual parasites to infectiousness
	p$dE 	<- 12
	p$tl 	<- 12.5

	# infectiousness to mosquitoes
	p$cD	<-	0.0676909
	p$cT	<-	0.0034482
	p$cU	<-	0.006203
	p$g_inf	<-	1.82425

	# anti-parasite immunity
	p$d1	<-	0.160527
	p$dd	<-	3650
	p$ID0	<-	1.577533
	p$kd	<-	0.476614
	p$ud	<-	9.44512
	p$ad0	<-	8001.99
	p$fd0	<-	0.007055
	p$gd	<-	4.8183

	# anti-infection immunity
	p$b0	<-	0.590076
	p$b1	<-	0.5
	p$db	<-	3650
	p$IB0	<-	43.8787
	p$kb	<-	2.15506
	p$ub	<-	7.19919

	# clinical immunity
	p$phi0	<-	0.791666
	p$phi1	<-	0.000737
	p$dc	<-	10950
	p$IC0	<-	18.02366
	p$kc	<-	2.36949
	p$uc	<-	6.06349
	p$PM	<-	0.774368
	p$dm	<-	67.6952

	# mosquito parameters
	p$tau	<-	10
	p$mu	<-	0.132
	p$f		<-	0.33333333
	p$Q0	<-	0.92
	
	return(p)
}

#------------------------------------------------
#' Set remaining parameters
#'
#' Replace some default parameters with values from a parameter vector v.
#'
#' @param v vector of parameters, with names corresponding to those returned by default_parameters().
#'
#' @export

set_parameters <- function(v){
	
    # start with default parameters
	p <- default_parameters()
    
    # replace values based on parameter name
	s <- names(v)
	for(i in 1:length(v)){
		p[s[i]] <- v[i]
	}
    
	return(p)
}

#------------------------------------------------
#' Draw random parameters
#'
#' Set some parameters as default, and draw others from distributions. TODO - find out where these distributions come from.
#'
#' @param n nuber of random parameter sets to draw.
#'
#' @export

random_parameters <- function(n=1)
{
    ret <-list()
    for(i in 1:n)
    {
        p <- list()
        
        # age, heterogeneity in exposure
        p$eta	<-	0.0001305
        p$rho	<-	0.85
        p$a0	<-	2920
        p$s2	<-	1.67
        
        # rate of leaving infection states
        p$rA	<-	0.00512821
        p$rT	<-	0.2
        p$rD	<-	0.2
        p$rU	<-	1/(rlnorm(1,4.0,0.31)*365)
        p$rP	<-	0.2
        
        # human latent period and time lag from asexual parasites to infectiousness
        p$dE 	<- 12
        p$tl 	<- 12.5
        
        # infectiousness to mosquitoes
        p$cD	<-	0.0676909
        p$cT	<-	0.0034482
        p$cU	<-	0.006203
        p$g_inf	<-	1.82425
        p$aA  <-  rbeta(1,1,1)
        p$aU  <-  rbeta(1,4,1)
        
        # anti-parasite immunity
        p$d1	<-	rbeta(1,10,30)
        p$dd	<-	3650
        p$ID0	<-	rlnorm(1,2.38,1.08)
        p$kd	<-	rlnorm(1,0.50,0.32)
        p$ud	<-	rlnorm(1,1.40,0.80)
        p$ad0	<-	rlnorm(1,2.02,0.43)*365
        p$fd0	<-	rbeta(1,1,1)
        p$gd	<-	rlnorm(1,0.22,0.64)
        
        # anti-infection immunity
        p$b0	<-	rbeta(1,1.3,1.3)
        p$b1	<-	0.5
        p$db	<-	3650
        p$IB0	<-	rlnorm(1,3.07,0.93)
        p$kb	<-	rlnorm(1,0.50,0.32)
        p$ub	<-	rlnorm(1,1.40,0.80)
        
        # clinical immunity
        p$phi0	<-	rbeta(1,8.3,2.1)
        p$phi1	<-	rbeta(1,1.1,2)
        p$dc	<-	10950
        p$IC0	<-	rlnorm(1,3.76,1.02)
        p$kc	<-	rlnorm(1,0.50,0.32)
        p$uc	<-	rlnorm(1,1.40,0.80)
        p$PM	<-	rbeta(1,1,1)
        p$dm	<-	rlnorm(1,5.26,0.33)
        
        # mosquito parameters
        p$tau	<-	10
        p$mu	<-	0.132
        p$f		<-	0.33333333
        p$Q0	<-	0.92
        
        ret[[i]] <- p
    }
    
    return(ret)
}

#------------------------------------------------
#' Default age range
#'
#' Returns vector of ages used by default in model fitting.
#'
#' @export

default_age <- function() {
    data(age)
    return(as.vector(unlist(age)))
}
