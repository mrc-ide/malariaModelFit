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
#' @param prior_dist the type of prior distribution to define.
#' @param prior_params parameters associated with the prior distribution.
#'
#' @export

define_prior <- function(project, name, prior_dist, prior_params) {
  
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
  
  # create copy of priors dataframe
  if (name_in_model_priors) {
    new_df <- project$model_priors
  } else {
    new_df <- project$fitting_priors
  }
  
  # update specified parameter
  w <- which(new_df$name == name)
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
  assert_in(c("name", "definition", "prior_dist", "prior_params"), names(parameters),
            message = "parameters dataframe must contain the following column names: {name, definition, prior_dist, prior_params}")
  
  # check format of every column
  assert_string(parameters$name)
  assert_string(parameters$definition)
  assert_in(parameters$prior_dist, c("fixed", "beta", "norm", "lnorm", "gamma"))
  assert_list(parameters$prior_params)
  
  # check that all parameters have the correct correct set of prior parameters
  if (any(parameters$prior_dist == "fixed")) {
    p_sub <- subset(parameters, parameters$prior_dist == "fixed")
    params_message <- "for fixed parameters, prior_params must contain the value of the parameter"
    apply(p_sub, 1, function(x) {
      assert_length(x$prior_params, 1, message = params_message)
      assert_numeric(x$prior_params, message = params_message)
    })
  }
  if (any(parameters$prior_dist == "beta")) {
    p_sub <- subset(parameters, parameters$prior_dist == "beta")
    params_message <- "beta distribution must have two prior_params values specifying the shape of the distribution. Both parameters must be in the interval (0,infinity)"
    apply(p_sub, 1, function(x) {
      assert_vector(x$prior_params, message = params_message)
      assert_length(x$prior_params, 2, message = params_message)
      assert_pos(x$prior_params, zero_allowed = FALSE, message2 = params_message)
    })
  }
  if (any(parameters$prior_dist == "norm")) {
    p_sub <- subset(parameters, parameters$prior_dist == "norm")
    params_message <- "normal distribution must have two prior_params values specifying the mean and standard deviation of the distribution. The mean must be in the interval (-infinity,infinity), and the standard deviation must be in the interval (0,infinity)"
    apply(p_sub, 1, function(x) {
      assert_vector_numeric(x$prior_params, message = params_message)
      assert_length(x$prior_params, 2, message = params_message)
      assert_pos(x$prior_params[2], zero_allowed = FALSE, message2 = params_message)
    })
  }
  if (any(parameters$prior_dist == "lnorm")) {
    p_sub <- subset(parameters, parameters$prior_dist == "lnorm")
    params_message <- "lognormal distribution must have two prior_params values specifying the mean and standard deviation of the normal distribution that is exponentiated to produce the lognormal distribution. The mean must be in the interval (-infinity,infinity), and the standard deviation must be in the interval (0,infinity)"
    apply(p_sub, 1, function(x) {
      assert_vector_numeric(x$prior_params, message = params_message)
      assert_length(x$prior_params, 2, message = params_message)
      assert_pos(x$prior_params[2], zero_allowed = FALSE, message2 = params_message)
    })
  }
  if (any(parameters$prior_dist == "gamma")) {
    p_sub <- subset(parameters, parameters$prior_dist == "gamma")
    params_message <- "gamma distribution must have two prior_params values specifying the shape and scale of the distribution. Both parameters must be in the interval (0,infinity)"
    apply(p_sub, 1, function(x) {
      assert_vector(x$prior_params, message = params_message)
      assert_length(x$prior_params, 2, message = params_message)
      assert_pos(x$prior_params, zero_allowed = FALSE, message2 = params_message)
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
  ret <- apply(x, 1, function(y) {
    if (y$prior_dist == "fixed") {
      return(y$prior_params)
    } else {
      
      # gamma distribution is defined in terms of scale, but rgamma needs in terms of rate
      if (y$prior_dist == "gamma") {
        y$prior_params[2] <- 1/y$prior_params[2]
      }
      
      # create expression string and evaluate
      rand_expression <- sprintf("r%s(1, %s)", y$prior_dist, paste(unlist(y$prior_params), collapse = ", "))
      ret <- eval(parse(text = rand_expression))
      return(ret)
    }
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

#------------------------------------------------
#' @title Create a C++ prior function string
#'
#' @description Given a project with model and fitting priors already defined,
#'   creates a C++ function in the form of a string that will return the joint
#'   density of these priors in log space. This string can then be passed to
#'   \code{drjacoby} and used in MCMC.
#'
#' @param project an object of class "mmfit_project" (see
#'   \code{?mmfit_project()}).
#'
#' @export

create_prior_string <- function(project) {
  
  # check inputs
  assert_custom_class(project, "mmfit_project")
  assert_non_null(project$model_priors, message = "no model priors defined")
  
  # extract model and fitting priors and combine
  prior_df <- rbind(project$model_priors, project$fitting_priors)
  
  # initialise first lines of function string
  s_vec <- "SEXP logprior(std::vector<double> params) {
double ret = 0;"
  
  # add all prior distributions into function string
  for (i in 1:nrow(prior_df)) {
    if (prior_df$prior_dist[i] != "fixed") {
      prior_expression <- sprintf("ret += R::d%s(params[%s], %s, true);",
                                  prior_df$prior_dist[i],
                                  i-1,
                                  paste(prior_df$prior_params[[i]], collapse = ", "))
      s_vec <- c(s_vec, prior_expression)
    }
  }
  
  # finalise function string
  s_vec <- c(s_vec, "return Rcpp::wrap(ret);\n}")
  s <- paste(s_vec, collapse = "\n")
  
  # return string
  return(s)
}

