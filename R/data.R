#' Details of data from Griffin et al 2014 refit
#' 
#' \describe{
#'   \item{study}{Dataset name (study/group)}
#'   \item{slide_k}{Number positive by slide microscopy}
#'   \item{n}{Number sampled for slide and PCR}
#'   \item{pcr_k}{Number positive by PCR}
#'   \item{clin_k}{Number of clinical events}
#'   \item{clin_p_years}{Number of person-years for clinical events}
#'   \item{sev_k}{Severe events}
#'   \item{sev_p_years}{Number of person-years for severe events}
#'   \item{age0}{Lower end of age group}
#'   \item{age1}{Upper end of age group}
#'   \item{clin_group}{Group for clinical data random effects by age}
#'   \item{country}{Study location}
#'   \item{eir_mean}{Mean of prior for log(EIR)}
#'   \item{eir_sd}{SD of prior for log(EIR)}
#'   \item{place_index}{Place index}
#'   \item{treated_alpha}{Beta distribution prior parameter for proportion treated}
#'   \item{treated_beta}{Beta distribution prior parameter for proportion treated}
#'   \item{clinical_case_detection}{Code for case detection method: 
#'     \itemize{
#'       \item{1}{ Daily active case detection (ACD)}
#'       \item{2}{ Weekly ACD}
#'       \item{3}{ Passive case detection (PCD)}
#'     }
#'   }
#'   \item{study_index}{Study index (random effects for clinical and prevelence detection operate at level of study)}
#' }
"study_data_2014"