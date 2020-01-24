
#' Run model fitting
#'
#' @param project A mmfit project containing description, model priors, fitting priors and data
#' @param ... Additional named parameters passed to the \code{\link[drjacoby]{run_mcmc}} function
#'
#' @return Project with fitting output
#' @export
run_fit <- function(project, ...){
  
  # check inputs
  assert_custom_class(project, "mmfit_project")
  assert_non_null(project$model_priors, message = "Project does not contain model priors")
  assert_non_null(project$fitting_priors, message = "Project does not contain fitting priors")
  assert_non_null(project$data, message = "Project does not contain data")
  
  # Prepare vector of data for use in DrJacoby
  pd <- prepare_data(project$data)
  x <- unlist(pd)
  
  # logLikelihood string for use in DrJacoby
  lL <- likelihood_string
  
  # logPrior string for use in DrJacoby
  lP <- create_prior_string(project)

  # Parameter dataframe
  df_params <- create_df_params(project)
  
  # Hacks !!
  df_params[df_params$name == "rD", "max"] <- 50 
  additional <- data.frame(name = c("EIR", "ft"), min = c(10, 0.5), max = c(10, 0.5))
  df_params <- dplyr::bind_rows(df_params, additional)

  # Run DrJacoby
  project$output_raw <- drjacoby::run_mcmc(data = x,
                                 df_params = df_params,
                                 loglike = lL,
                                 logprior = lP,
                                 ...)
  
  return(project)
}