
#------------------------------------------------
#' @title Return names of required equilibrium parameters
#'   
#' @description Return names of required equilibrium parameters. Note that this 
#'   is a subset of all parameters, as some parameters are only needed in a 
#'   dynamic setting.
#'   
#' @export

required_eq_params <- function() {
  c("eta", "rho", "a0", "ub", "db", "b0", "b1", "IB0", "kb", "uc", 
    "dc", "IC0", "kc", "ud", "dd", "fd0", "ad0", "gd", "d1", "ID0", 
    "kd", "cU", "cD", "g_inf", "PM", "dm", "phi0", "phi1", "rT", 
    "rD", "rA", "rU", "rP", "aA", "aU", "cT", "s2", "f", "Q0")
}

#------------------------------------------------
#' @title Return names of all model parameters
#' 
#' @description Return names of all model parameters.
#'   
#' @export

all_params <- function() {
  c(required_eq_params(), "dE", "tl", "tau", "mu")
}

#------------------------------------------------
#' @title Default age range
#' 
#' @description Returns vector of ages used by default in model fitting.
#'   
#' @export

default_age <- function() {
    age <- read.table(malariaModelFit_file("age_distributions/age.txt"))
    ret <- as.vector(unlist(age))
    return(ret)
}

#------------------------------------------------
#' @title Load parameter set
#'   
#' @description Load a named parameter set from the inst/parameter_sets folder
#'   of this package. Input is filtered to relevant parameters only, and
#'   returned as a named list.
#'   
#' @param name name of parameter set to load
#'   
#' @export

load_parameters <- function(name = "parameters_Griffin2014.txt") {
  
  # read in parameters file
  name <- paste0("parameter_sets/", name)
  params <- read.table(malariaModelFit_file(name), sep="\t")
  
  # check format
  assert_dataframe(params)
  assert_eq(ncol(params), 2)
  names(params) <- c("name", "value")
  check_names(params$name, all_params())
  
  # subset to relevant parameters only
  params <- subset(params, name %in% all_params())
  
  # coerce to list
  ret <- as.list(as.numeric(as.character(params$value)))
  names(ret) <- as.character(params$name)
  return(ret)
}

#------------------------------------------------
#' @title Parameter definitions
#'   
#' @description Returns list of definitions for each parameter.
#'   
#' @export

parameter_definitions <- function() {
  
  p <- list()
  
  # age, heterogeneity in exposure
  p$eta	<-	"inverse of average age, assumed to be exponentially distributed (1/days)"
  p$rho	<-	"age-dependent biting parameter"
  p$a0	<-	"age-dependent biting parameter (days)"
  p$s2	<-	"variance of log-heterogeneity in biting rates"
  
  # rate of leaving infection states
  p$rA	<-	"inverse duration patent infection (1/days)"
  p$rT	<-	"inverse duration treated clinical disease (1/days)"
  p$rD	<-	"inverse duration untreated clinical disease (1/days)"
  p$rU	<-	"inverse duration sub-patent infection (1/days)"
  p$rP	<-	"inverse duration prophylactic period (1/days)"
  
  # human latent period and time lag from asexual parasites to infectiousness
  p$dE 	<- "human latent period (days)"
  p$tl 	<- "delay from emergence of parasites to onward infectivity"
  
  # infectiousness to mosquitoes
  p$cD	<-	"infectiousness to mosquitoes with no immunity, or with untreated disease"
  p$cT	<-	"infectiousness to mosquitoes after treatment"
  p$cU	<-	"infectiousness to mosquitoes in sub-patent infection"
  p$g_inf	<-	"relates infectiousness to mosquitoes to probability of detection"
  
  # anti-parasite immunity
  p$d1	<-	"anti-parasite immunity: probability with maximum immunity"
  p$dd	<-	"anti-parasite immunity: inverse decay rate, or duration of immunity (days)"
  p$ID0	<-	"anti-parasite immunity: scale parameter"
  p$kd	<-	"anti-parasite immunity: shape parameter"
  p$ud	<-	"anti-parasite immunity: duration in which immunity is not boosted (days)"
  p$ad0	<-	"anti-parasite immunity: scale parameter relating age to immunity"
  p$fd0	<-	"anti-parasite immunity: parameter relating age to immunity"
  p$gd	<-	"anti-parasite immunity: shape parameter relating age to immunity"
  
  # detectibility by PCR
  p$aA <- "PCR detectibility parameter for asymptomatic infection"
  p$aU <- "PCR detectibility parameter for sub-patent infection"
  
  # anti-infection immunity
  p$b0	<-	"anti-infection immunity: probability with no immunity"
  p$b1	<-	"anti-infection immunity: maximum relative reduction"
  p$db	<-	"anti-infection immunity: inverse decay rate, or duration of immunity (days)"
  p$IB0	<-	"anti-infection immunity: scale parameter"
  p$kb	<-	"anti-infection immunity: shape parameter"
  p$ub	<-	"anti-infection immunity: duration in which immunity is not boosted (days)"
  
  # clinical immunity
  p$phi0	<-	"probability with no immunity"
  p$phi1	<-	"maximum relative reduction"
  p$dc	<-	"inverse decay rate, or duration of immunity (days)"
  p$IC0	<-	"scale parameter"
  p$kc	<-	"shape parameter"
  p$uc	<-	"duration in which immunity is not boosted"
  p$PM	<-	"new-born immunity relative to mother's"
  p$dm	<-	"inverse decay rate of maternal immunity (days)"
  
  # mosquito parameters
  p$tau	<-	"extrinsic incubation period (days)"
  p$mu	<-	"mosquito death rate (1/life expectancy in days)"
  p$f		<-	"feeding rate (1/gonotrophic cycle length in days)"
  p$Q0	<-	"Human Blood Index; proportion bloodmeals taken on humans"
  
  return(p)
}

#------------------------------------------------
#' @title Default parameters
#'   
#' @description Load the current default parameter set.
#'   
#' @export

default_parameters <- function() {
  load_parameters("parameters_Xiaoyu.txt")
}

#------------------------------------------------
#' @title Set remaining parameters
#'   
#' @description Replace some default parameters with values from a parameter
#'   vector v.
#'   
#' @param v vector of parameters, with names corresponding to those returned by
#'   \code{default_parameters()}.
#'   
#' @export

set_parameters <- function(v) {
  
  # check inputs
  s <- names(v)
  assert_in(s, all_params())
  
  # start with default parameters
  p <- default_parameters()
  
  # replace values based on parameter name
  for (i in 1:length(v)) {
    p[s[i]] <- v[i]
  }
  
  return(p)
}

#------------------------------------------------
#' @title Draw random parameters
#'   
#' @description Set some parameters as default, and draw others from
#'   distributions. TODO - find out where these distributions and values came
#'   from.
#'   
#' @param n nuber of random parameter sets to draw.
#'   
#' @export

random_parameters <- function(n = 1) {
  
  ret <-list()
  
  for(i in 1:n) {
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
