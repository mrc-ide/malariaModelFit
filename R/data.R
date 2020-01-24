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
#' }
"study_data_2014"