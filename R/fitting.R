#------------------------------------------------
#' @title Load model priors into a project
#'
#' @description Model priors are stored within the package
#'   inst/extdata/model_priors folder. Load one of these objects by name, and
#'   attach to an existing project.
#'
#' @param project an object of class "mmfit_project" (see
#'   \code{?mmfit_project()}).
#' @param file_name the name of a file within the inst/extdata/model_priors
#'   folder.
#' 
#' @export

load_model_priors <- function(project, file_name = "refit2020_model_priors.rds") {
  
  # check inputs
  assert_custom_class(project, "mmfit_project")
  assert_single_string(file_name)
  
  # load model parameters from inst/extdata/model_priors folder
  params_df <- mmfit_file(paste0("model_priors/", file_name))
  
  # check parameters
  check_priors(params_df)
  
  # add to project and return
  project$model_priors <- params_df
  return(project)
}

#------------------------------------------------
#' @title Load fitting priors into a project
#'
#' @description Priors involved in fitting to the data are stored within the
#'   package inst/extdata/fitting_priors folder. Load one of these objects by
#'   name, and attach to an existing project.
#'
#' @param project an object of class "mmfit_project" (see
#'   \code{?mmfit_project()}).
#' @param file_name the name of a file within the inst/extdata/fitting_priors
#'   folder.
#' 
#' @export

load_fitting_priors <- function(project, file_name = "refit2020_fitting_priors.rds") {
  
  # check inputs
  assert_custom_class(project, "mmfit_project")
  assert_single_string(file_name)
  
  # load model parameters from inst/extdata/fitting_priors folder
  params_df <- mmfit_file(paste0("fitting_priors/", file_name))
  
  # check parameters
  check_priors(params_df)
  
  # add to project and return
  project$fitting_priors <- params_df
  return(project)
}

#------------------------------------------------
#' @title Define a prior distribution over a parameter within a project 
#'
#' @description Given a model fitting project with model parameters and/or
#'   fitting parameters loaded, redefine the prior distribution of a named
#'   parameter. This function will try to find the given parameter name within
#'   both \code{project$model_params} and \code{project$fitting_params}, hence
#'   this function can be used to specify priors over both objects.
#'
#' @param project an object of class "mmfit_project" (see
#'   \code{?mmfit_project()}).
#' @param name the name of a parameter within the project model parameters or
#'   fitting parameters.
#' @param prior_dist the type of prior distribution to define. One of {"fixed",
#'   "dbeta", "dnorm", "dlnorm", "dgamma"}.
#' @param prior_params parameters associated with the prior distribution.
#' @param range_min the minimum allowed parameter value.
#' @param range_max the maximum allowed parameter value.
#'
#' @export

define_prior <- function(project, name, prior_dist, prior_params, range_min = NULL, range_max = NULL) {
  
  # check inputs
  assert_custom_class(project, "mmfit_project")
  model_priors_loaded <- !is.null(project$model_priors)
  fitting_priors_loaded <- !is.null(project$fitting_priors)
  if (!model_priors_loaded && !fitting_priors_loaded) {
    stop("no model priors loaded")
  }
  assert_single_string(name)
  name_in_model_priors <- FALSE
  if (model_priors_loaded) {
    name_in_model_priors <- (name %in% project$model_priors$name)
  }
  name_in_fitting_priors <- FALSE
  if (fitting_priors_loaded) {
    name_in_fitting_priors <- (name %in% project$fitting_priors$name)
  }
  if (!name_in_model_priors && !name_in_fitting_priors) {
    stop(sprintf("could not find parameter %s in model priors or fitting priors", name))
  }
  if (name_in_model_priors && name_in_fitting_priors) {
    stop(sprintf("parameter %s found in both model priors and fitting priors; cannot continue", name))
  }
  assert_single_string(prior_dist)
  assert_in(prior_dist, c("fixed", "dbeta", "dnorm", "dlnorm", "dgamma"))
  
  # define default ranges and perform checks on ranges
  if (is.null(range_min)) {
    if (prior_dist == "fixed") {
      stop("range_min and range_max must be specified for fixed parameters")
    } else if (prior_dist == "dbeta") {
      range_min <- 0
      range_max <- 1
    } else if (prior_dist == "dnorm") {
      range_min <- -Inf
      range_max <- Inf
    } else if (prior_dist %in% c("dlnorm", "dgamma")) {
      range_min <- 0
      range_max <- Inf
    }
  }
  assert_single_numeric(range_min)
  assert_single_numeric(range_max)
  
  # create copy of priors dataframe
  if (name_in_model_priors) {
    new_df <- project$model_priors
  } else {
    new_df <- project$fitting_priors
  }
  
  # update specified parameter
  w <- which(new_df$name == name)
  new_df$range_min[w] <- range_min
  new_df$range_max[w] <- range_max
  new_df$prior_dist[w] <- prior_dist
  new_df$prior_params[w] <- list(prior_params)
  
  # perform checks on new priors dataframe
  check_priors(new_df)
  
  # load new priors dataframe into project
  if (name_in_model_priors) {
    project$model_priors <- new_df
  } else {
    project$fitting_priors <- new_df
  }
  
  # return project
  return(project)
}

#------------------------------------------------
# perform checks on prior distributions
#' @noRd

check_priors <- function(parameters) {
  
  # check that dataframe with the correct columns
  assert_dataframe(parameters)
  assert_in(c("name", "definition", "range_min", "range_max", "prior_dist", "prior_params"), names(parameters),
            message = "parameters dataframe must contain the following column names: {name, definition, range_min, range_max, prior_dist, prior_params}")
  
  # check format of every column
  assert_string(parameters$name)
  assert_string(parameters$definition)
  assert_numeric(parameters$range_min)
  assert_numeric(parameters$range_max)
  assert_in(parameters$prior_dist, c("fixed", "dbeta", "dnorm", "dlnorm", "dgamma"))
  assert_list(parameters$prior_params)
  
  # check that range_max is greater than or equal to range_min
  assert_greq(parameters$range_max, parameters$range_min,
              message = "parameters$range_max must be greater than or equal to parameters$range_min")
  
  # check that all parameters have the correct range and correct set of prior
  # parameters
  if (any(parameters$prior_dist == "fixed")) {
    p_sub <- subset(parameters, parameters$prior_dist == "fixed")
    apply(p_sub, 1, function(x) {
      assert_eq(x$range_min, x$range_max,
                message = "fixed parameters must have the same value for range_min and range_max")
      assert_null(x$prior_params, message = "prior_params must be NULL for fixed parameters")
    })
  }
  if (any(parameters$prior_dist == "dbeta")) {
    p_sub <- subset(parameters, parameters$prior_dist == "dbeta")
    dnorm_range_message <- "dbeta parameters must have range 0 to 1"
    dnorm_params_message <- "dbeta parameters must have two prior_params values specifying the shape of the distribution"
    apply(p_sub, 1, function(x) {
      assert_eq(x$range_min, 0, message = dnorm_range_message)
      assert_eq(x$range_max, 1, message = dnorm_range_message)
      assert_vector(x$prior_params, message = dnorm_params_message)
      assert_length(x$prior_params, 2, message = dnorm_params_message)
      assert_pos(x$prior_params, zero_allowed = FALSE, message2 = dnorm_params_message)
    })
  }
  if (any(parameters$prior_dist == "dnorm")) {
    p_sub <- subset(parameters, parameters$prior_dist == "dnorm")
    dnorm_range_message <- "dnorm parameters must have range -Inf to Inf"
    dnorm_params_message <- "dnorm parameters must have two prior_params values specifying the mean and standard deviation of the distribution"
    apply(p_sub, 1, function(x) {
      assert_eq(x$range_min, -Inf, message = dnorm_range_message)
      assert_eq(x$range_max, Inf, message = dnorm_range_message)
      assert_vector(x$prior_params, message = dnorm_params_message)
      assert_length(x$prior_params, 2, message = dnorm_params_message)
      assert_numeric(x$prior_params, message = dnorm_params_message)
      assert_pos(x$prior_params[2], zero_allowed = FALSE, message2 = dnorm_params_message)
    })
  }
  if (any(parameters$prior_dist == "dlnorm")) {
    p_sub <- subset(parameters, parameters$prior_dist == "dlnorm")
    dnorm_range_message <- "dlnorm parameters must have range 0 to Inf"
    dnorm_params_message <- "dlnorm parameters must have two prior_params values specifying the mean and standard deviation of the distribution prior to exponentiating"
    apply(p_sub, 1, function(x) {
      assert_eq(x$range_min, 0, message = dnorm_range_message)
      assert_eq(x$range_max, Inf, message = dnorm_range_message)
      assert_vector(x$prior_params, message = dnorm_params_message)
      assert_length(x$prior_params, 2, message = dnorm_params_message)
      assert_numeric(x$prior_params, message = dnorm_params_message)
      assert_pos(x$prior_params[2], zero_allowed = FALSE, message2 = dnorm_params_message)
    })
  }
  if (any(parameters$prior_dist == "dgamma")) {
    p_sub <- subset(parameters, parameters$prior_dist == "dgamma")
    dnorm_range_message <- "dgamma parameters must have range 0 to Inf"
    dnorm_params_message <- "dgamma parameters must have two prior_params values specifying the shape and rate of the distribution"
    apply(p_sub, 1, function(x) {
      assert_eq(x$range_min, 0, message = dnorm_range_message)
      assert_eq(x$range_max, Inf, message = dnorm_range_message)
      assert_vector(x$prior_params, message = dnorm_params_message)
      assert_length(x$prior_params, 2, message = dnorm_params_message)
      assert_pos(x$prior_params, zero_allowed = FALSE, message2 = dnorm_params_message)
    })
  }
  
}

#------------------------------------------------
#' @title Draw a set of model parameters from the prior
#'
#' @description Takes a fitting project with model priors already loaded, and
#'   draws a new parameter set from the prior.
#'
#' @param project an object of class "mmfit_project" (see
#'   \code{?mmfit_project()}), with model parameters already loaded.
#'
#' @export

draw_model_prior <- function(project) {
  
  # check inputs
  assert_custom_class(project, "mmfit_project")
  assert_non_null(project$model_priors)
  
  # draw from priors
  draw_prior(project$model_priors)
}

#------------------------------------------------
#' @title Draw a set of fitting parameters from the prior
#'
#' @description Takes a fitting project with fitting priors already loaded, and
#'   draws a new parameter set from the prior.
#'
#' @param project an object of class "mmfit_project" (see
#'   \code{?mmfit_project()}), with fitting parameters already loaded.
#'
#' @export

draw_fitting_prior <- function(project) {
  
  # check inputs
  assert_custom_class(project, "mmfit_project")
  assert_non_null(project$fitting_priors)
  
  # draw from priors
  draw_prior(project$fitting_priors)
}

#------------------------------------------------
# general function for drawing from any priors dataframe
#' @importFrom stats rbeta rnorm rlnorm rgamma
#' @noRd

draw_prior <- function(x) {
  
  # draw from specified distributions
  ret <- apply(x, 1, function(x) {
    switch(x$prior_dist,
           "fixed" = x$range_min,
           "dbeta" = rbeta(1, x$prior_params[1], x$prior_params[2]),
           "dnorm" = rnorm(1, x$prior_params[1], x$prior_params[2]),
           "dlnorm" = rlnorm(1, x$prior_params[1], x$prior_params[2]),
           "dgamma" = rgamma(1, x$prior_params[1], x$prior_params[2])
           )
  })
  
  # return as mmfit_params class
  ret <- as.list(ret)
  names(ret) <- x$name
  class(ret) <- "mmfit_params"
  
  return(ret)
}

#------------------------------------------------
#' @title Load data into a project
#'
#' @description Data for use in model fitting is stored within the package
#'   inst/extdata/data folder. Load a data object from this folder by name, and
#'   attach to an existing project.
#'
#' @param project an object of class "mmfit_project" (see
#'   \code{?mmfit_project()}).
#' @param file_name the name of a file within the inst/extdata/data folder.
#'
#' @export

load_data <- function(project, file_name = "refit2020_data.rds") {
  
  # check inputs
  assert_custom_class(project, "mmfit_project")
  assert_single_string(file_name)
  
  # load data from inst/extdata/data folder
  data_df <- mmfit_file(paste0("data/", file_name))
  
  # check data
  check_data(data_df)
  
  # add to project and return
  project$data <- data_df
  return(project)
}

#------------------------------------------------
# perform checks on data format
#' @noRd

check_data <- function(dat) {
  # TODO - some checks on data
}

