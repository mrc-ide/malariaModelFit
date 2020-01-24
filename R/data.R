#' Details of data from Griffin et al 2014 refit
#' 
#' \describe{
#'   \item{country_name}{Country name (metadata not used in fitting)}
#'   \item{site_name}{Site name (metadata not used in fitting)}
#'   \item{reference}{Citation (metadata not used in fitting)}
#'   \item{study_index}{Unique numerical index for each study, used for study-level random effects}
#'   \item{site_index}{Unique numerical index for each site, used for site-level random effects}
#'   \item{numer}{Data numerator. This will either be number +ve by microscopy or number of cases}
#'   \item{denom}{Data denominator. This will either be number tested or person-years}
#'   \item{type}{Data type:
#'     \itemize{
#'       \item{1}{ Incidence data}
#'       \item{2}{ Prevalence data}
#'     }}   
#'   \item{age0}{Lower end of age group}
#'   \item{age1}{Upper end of age group}
#'   \item{case_detection}{Code for case detection method: 
#'     \itemize{
#'       \item{1}{ Daily active case detection (ACD)}
#'       \item{2}{ Weekly ACD}
#'       \item{3}{ Passive case detection (PCD)}
#'     }
#'   }
#'   \item{age_bracket}{Age bracket:  used to calculated age-group random effects. Where the index
#'   indicated that the mid point of the study age group falls in
#'     \itemize{
#'       \item{1}{ 0 - 2 years old}
#'       \item{2}{ 2 - 5 years old}
#'       \item{3}{ 5 - 15 years old}
#'       \item{4}{ 15+  years old}
#'     }
#'   }
#' }
"study_data_2014"

#' Extracted string version of the c++ loglikelihood functon for use in drjacoby
"likelihood_string"