#' Markov Chain Monte Carlo (MCMC) for time-varying parameter survival models with shrinkage
#'
#' \code{shrinkDSM} samples from the joint posterior distribution of the parameters of a time-varying
#' parameter survival model with shrinkage and returns the MCMC draws.
#' See also \code{\link[shrinkTVP]{shrinkTVP}} to see more examples of how to modify the prior setup of the time-varying
#' component of the model.
#'
#' @param formula an object of class "formula": a symbolic representation of the model, as in the
#' function \code{lm}. For details, see \code{\link{formula}}. The response may by a survival object as returned by the \code{Surv} function.
#' @param data \emph{optional} data frame containing the response variable and the covariates. If not found in \code{data},
#' the variables are taken from \code{environment(formula)}, typically the environment from which \code{shrinkDSM}
#' is called. No \code{NA}s are allowed in the response variable and the covariates.
#' @param mod_type character string that reads either \code{"triple"}, \code{"double"} or \code{"ridge"}.
#' Determines whether the triple gamma, double gamma or ridge prior are used for \code{theta_sr} and \code{beta_mean}.
#' The default is "double".
#' @param delta The status indicator of the last period, 0 = censored or 1 = event observed. If the \code{formula} is a survival object this parameter is ignored.
#' @param S integer vector of time points that start a new interval.
#' Parameters are fixed within an interval and vary across intervals.
#' @param group \emph{optional} grouping indicator for latent factor.
#' @param subset \emph{optional} vector specifying a subset of observations to be used in the fitting process.
#' @param niter positive integer, indicating the number of MCMC iterations
#' to perform, including the burn-in. Has to be larger than or equal to \code{nburn} + 2. The default value is 10000.
#' @param nburn non-negative integer, indicating the number of iterations discarded
#' as burn-in. Has to be smaller than or equal to \code{niter} - 2. The default value is \code{round(niter / 2)}.
#' @param nthin positive integer, indicating the degree of thinning to be performed. Every \code{nthin} draw is kept and returned.
#' The default value is 1, implying that every draw is kept.
#' @param learn_a_xi  logical value indicating whether to learn a_xi, the spike parameter of the state variances.
#' The default value is \code{TRUE}.
#' @param learn_a_tau logical value indicating whether to learn a_tau, the spike parameter of the mean of the
#' initial values of the states. The default value is \code{TRUE}.
#' @param a_xi positive, real number, indicating the (fixed) value for a_xi. Ignored if
#' \code{learn_a_xi} is \code{TRUE} or \code{mod_type} is set to \code{"ridge"}. The default value is 0.1.
#' @param a_tau positive, real number, indicating the (fixed) value for a_tau. Ignored if
#' \code{learn_a_tau} is \code{TRUE} or \code{mod_type} is set to \code{"ridge"}. The default value is 0.1.
#' @param learn_c_xi logical value indicating whether to learn c_xi, the tail parameter of the state variances.
#' Ignored if \code{mod_type} is not set to \code{"triple"} or \code{a_eq_c_xi} is set to \code{TRUE}.
#' The default value is \code{TRUE}.
#' @param learn_c_tau logical value indicating whether to learn c_tau, the tail parameter of the mean of the
#' initial values of the states. Ignored if \code{mod_type} is not set to \code{"triple"} or \code{a_eq_c_tau} is set to \code{TRUE}.
#' The default value is \code{TRUE}.
#' @param c_xi positive, real number, indicating the (fixed) value for c_xi. Ignored if
#' \code{learn_c_xi} is \code{TRUE}, \code{mod_type} is not set to \code{"triple"} or \code{a_eq_c_xi} is set to \code{TRUE}.
#' The default value is 0.1.
#' @param c_tau positive, real number, indicating the (fixed) value for c_tau. Ignored if
#' \code{learn_c_xi} is \code{TRUE}, \code{mod_type} is not set to \code{"triple"}  or \code{a_eq_c_tau} is set to \code{TRUE}.
#' The default value is 0.1.
#' @param a_eq_c_xi logical value indicating whether to force \code{a_xi} and \code{c_xi} to be equal.
#' If set to \code{TRUE}, \code{beta_a_xi} and \code{alpha_a_xi} are used as the hyperparameters and \code{beta_c_xi} and \code{alpha_c_xi} are ignored.
#' Ignored if \code{mod_type} is not set to \code{"triple"}. The default value is \code{FALSE}.
#' @param a_eq_c_tau logical value indicating whether to force \code{a_tau} and \code{c_tau} to be equal.
#' If set to \code{TRUE}, \code{beta_a_tau} and \code{alpha_a_tau} are used as the hyperparameters and \code{beta_c_tau} and \code{alpha_c_tau} are ignored.
#' Ignored if \code{mod_type} is not set to \code{"triple"}. The default value is \code{FALSE}.
#' @param learn_kappa2_B logical value indicating whether to learn kappa2_B, the global level of shrinkage for
#' the state variances. The default value is \code{TRUE}.
#' @param learn_lambda2_B logical value indicating whether to learn the lambda2_B parameter,
#' the global level of shrinkage for the mean of the initial values of the states. The default value is \code{TRUE}.
#' @param kappa2_B positive, real number. Initial value of kappa2_B. The default value is 20.
#' @param lambda2_B positive, real number. Initial value of lambda2_B. The default value is
#' @param hyperprior_param \emph{optional} named list containing hyperparameter values.
#' Not all have to be supplied, with those missing being replaced by the default values.
#' Any list elements that are misnamed will be ignored and a warning will be thrown.
#' All hyperparameter values have to be positive, real numbers. The following hyperparameters can be
#' supplied:
#' \itemize{
#' \item \code{e1}: The default value is 0.001.
#' \item \code{e2}: The default value is 0.001.
#' \item \code{d1}: The default value is 0.001.
#' \item \code{d2}: The default value is 0.001.
#' \item \code{beta_a_xi}: The default value is 10.
#' \item \code{beta_a_tau}: The default value is 10.
#' \item \code{alpha_a_xi}: The default value is 5.
#' \item \code{alpha_a_tau}: The default value is 5.
#' }
#' @param sv_param \emph{optional} named list containing hyperparameter values for the stochastic volatility
#' parameters. Not all have to be supplied, with those missing being replaced by the default values.
#' Any list elements that are misnamed will be ignored and a warning will be thrown. Ignored if
#' \code{group} is missing. The following elements can be supplied:
#' \itemize{
#' \item \code{Bsigma_sv}: positive, real number. The default value is 1.
#' \item \code{a0_sv}: positive, real number. The default value is 5.
#' \item \code{b0_sv}: positive, real number. The default value is 1.5.
#' }
#' @param MH_tuning \emph{optional} named list containing values used to tune the MH steps for \code{a_xi} and \code{a_tau}. Not all have to be supplied, with those missing being replaced by the default values.
#' Any list elements that are misnamed will be ignored and a warning will be thrown.
#' The arguments for \code{a_xi}(\code{a_tau}) are only used if \code{learn_a_xi}(\code{learn_a_tau})
#' is set to \code{TRUE}. Arguments ending in "adaptive" are
#' logical values indicating whether or not to make the MH step for the respective parameter adaptive. Arguments ending in "tuning_par" serve two different purposes.
#' If the respective MH step is not set to be adaptive, it acts as the standard deviation of the proposal distribution. If the respective MH step
#' is set to be adaptive, it acts as the initial standard deviation. Arguments ending in "target_rate" define the acceptance rate the algorithm aims to achieve.
#' Arguments ending in "max_adapt" set the maximum value by which the logarithm of the standard deviation of the proposal distribution is adjusted. Finally,
#' arguments ending in "batch_size" set the batch size after which the standard deviation of the proposal distribution is adjusted.
#' The following elements can be supplied:
#' \itemize{
#' \item \code{a_xi_adaptive}: logical value. The default is \code{TRUE}.
#' \item \code{a_xi_tuning_par}: positive, real number. The default value is 1.
#' \item \code{a_xi_target_rate}: positive, real number, between 0 and 1. The default value is 0.44.
#' \item \code{a_xi_max_adapt}: positive, real number. The default value is 0.01.
#' \item \code{a_xi_batch_size}: positive integer. The default value is 50.
#' \item \code{a_tau_adaptive}: logical value. The default is \code{TRUE}.
#' \item \code{a_tau_tuning_par}: positive, real number. The default value is 1.
#' \item \code{a_tau_target_rate}: positive, real number, between 0 and 1. The default value is 0.44.
#' \item \code{a_tau_max_adapt}: positive, real number. The default value is 0.01.
#' \item \code{a_tau_batch_size}: positive integer. The default value is 50.
#' }
#' @param phi_param \emph{optional} named list containing hyperparameter values for the grouped factor
#' and values to tune the MH steps for \code{a_phi} and \code{c_phi}. Not all have to be supplied, with
#' those missing being replaced by the default values. Any list elements that are misnamed will be ignored
#' and a warning will be thrown. Ignored if \code{group} is missing. The following elements can be supplied:
#' \itemize{
#' \item \code{mod_type_phi} character string that reads either \code{"triple"}, \code{"double"} or \code{"ridge"}. Determines whether the triple gamma, double gamma or ridge prior are used for \code{phi}. The default is "double".
#' \item \code{learn_a_phi}: logical value. The default is \code{TRUE}.
#' \item \code{a_phi}: positive, real number. The default value is 0.1.
#' \item \code{learn_c_phi}: logical value. The default is \code{TRUE}.
#' \item \code{c_phi}: positive, real number. The default value is  0.1,
#' \item \code{a_phi_eq_c_phi}: logical value. The default is \code{FALSE}.
#' \item \code{learn_lambda2_B_phi}: logical value. The default is \code{TRUE}.
#' \item \code{lambda2_B_phi}: positive, real number. The default value is 20.
#' \item \code{e1_phi}: positive, real number. The default value is 0.001.
#' \item \code{e2_phi}: positive, real number. The default value is 0.001.
#' \item \code{beta_a_phi}: positive, real number. The default value is 10.
#' \item \code{alpha_a_phi}: positive, real number. The default value is 5.
#' \item \code{beta_c_phi}: positive, real number. The default value is 10.
#' \item \code{alpha_c_phi}: positive, real number. The default value is 5.
#' \item \code{a_phi_adaptive}: logical value. The default is \code{TRUE}.
#' \item \code{a_phi_tuning_par}: positive, real number. The default value is 1.
#' \item \code{a_phi_target_rate}: positive, real number, between 0 and 1. The default value is 0.44.
#' \item \code{a_phi_max_adapt}: positive, real number. The default value is 0.01.
#' \item \code{a_phi_batch_size}: positive integer. The default value is 50.
#' \item \code{c_phi_adaptive}: logical value. The default is \code{TRUE}.
#' \item \code{c_phi_tuning_par}: positive, real number. The default value is 1.
#' \item \code{c_phi_target_rate}:  positive, real number, between 0 and 1. The default value is 0.44.
#' \item \code{c_phi_max_adapt}: positive, real number. The default value is 0.01.
#' \item \code{c_phi_batch_size}: positive integer. The default value is 50.
#' }
#'
#' @param display_progress logical value indicating whether the progress bar and other informative output should be
#' displayed. The default value is \code{TRUE}.
#'
#' @return  The value returned is a list object of class \code{shrinkDSM} containing
#' \item{\code{beta}}{\code{list} object containing an \code{mcmc.dsm.tvp} object for the parameter draws from the posterior distribution of the centered
#' states, one for each covariate. In the case that there is only one covariate, this becomes just a single \code{mcmc.dsm.tvp} object.}
#' \item{\code{beta_mean}}{\code{mcmc} object containing the parameter draws from the posterior distribution of beta_mean.}
#' \item{\code{theta_sr}}{\code{mcmc} object containing the parameter draws from the posterior distribution of the square root of theta.}
#' \item{\code{tau2}}{\code{mcmc} object containing the parameter draws from the posterior distribution of tau2.}
#' \item{\code{xi2}}{\code{mcmc} object containing the parameter draws from the posterior distribution of xi2.}
#' \item{\code{lambda2}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of lambda2.
#' Not returned if \code{mod_type} is not \code{"triple"}.}
#' \item{\code{kappa2}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of kappa2.
#' Not returned if \code{mod_type} is not \code{"triple"}.}
#' \item{\code{a_xi}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of a_xi.
#' Not returned if \code{learn_a_xi} is \code{FALSE} or \code{mod_type} is \code{"ridge"}.}
#' \item{\code{a_tau}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of a_tau.
#' Not returned if \code{learn_a_tau} is \code{FALSE} or \code{mod_type} is \code{"ridge"}.}
#' \item{\code{c_xi}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of c_xi.
#' Not returned if \code{learn_c_xi} is \code{FALSE} or \code{mod_type} is not \code{"triple"}.}
#' \item{\code{c_tau}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of c_tau.
#' Not returned if \code{learn_c_tau} is \code{FALSE} or \code{mod_type} is not \code{"triple"}.}
#' \item{\code{lambda2_B}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of lambda2_B.
#' Not returned if \code{learn_lambda2_B} is \code{FALSE} or \code{mod_type} is \code{"ridge"}.}
#' \item{\code{kappa2_B}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of kappa2_B.
#' Not returned if \code{learn_kappa2_B} is \code{FALSE} or \code{mod_type} is \code{"ridge"}.}
#' \item{\code{MH_diag}}{\emph{(optional)} named list containing statistics for assessing MH performance. Not returned if no MH steps are required
#' or none of them are specified to be adaptive.}
#' \item{\code{priorvals}}{\code{list} object containing hyperparameter values of the prior distributions, as specified by the user.}
#' \item{\code{model}}{\code{list} object containing the model matrix, model response and formula used.}
#' \item{\code{summaries}}{\code{list} object containing a collection of summary statistics of the posterior draws.}
#'
#' To display the output, use \code{plot} and \code{summary}. The \code{summary} method displays the specified prior values stored in
#' \code{priorvals} and the posterior summaries stored in \code{summaries}, while the \code{plot} method calls \code{coda}'s \code{plot.mcmc}
#' or the \code{plot.mcmc.dsm.tvp} method. Furthermore, all functions that can be applied to \code{coda::mcmc} objects
#' (e.g. \code{coda::acfplot}) can be applied to all output elements that are \code{coda} compatible.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' data("gastric")
#'
#' # Create intervals for piecewise exponential model
#' intervals <- divisionpoints(gastric$time, gastric$status, 2)
#'
#' # Estimate baseline model
#' mod <- shrinkDSM(time ~ radiation, gastric,
#'   delta = gastric$status, S = intervals
#' )
#'
#' # Alternative formula interface
#' mod_surv <- shrinkDSM(Surv(time, status) ~ radiation, gastric,
#'   S = intervals
#' )
#'
#' # Estimate model with different prior setup
#' mod2 <- shrinkDSM(time ~ radiation, gastric,
#'   delta = gastric$status, S = intervals,
#'   mod_type = "triple"
#' )
#'
#' # Change some of the hyperparameters
#' mod3 <- shrinkDSM(time ~ radiation, gastric,
#'   delta = gastric$status, S = intervals,
#'   mod_type = "triple",
#'   hyperprior_param = list(
#'     beta_a_xi = 5,
#'     alpha_a_xi = 10
#'   )
#' )
#' }
#' @author Daniel Winkler \email{daniel.winkler@@wu.ac.at}
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @export
shrinkDSM <- function(formula,
                      data,
                      mod_type = "double",
                      delta,
                      S,
                      group,
                      subset,
                      niter = 10000,
                      nburn = round(niter / 2),
                      nthin = 1,
                      learn_a_xi = TRUE,
                      learn_a_tau = TRUE,
                      a_xi = 0.1,
                      a_tau = 0.1,
                      learn_c_xi = TRUE,
                      learn_c_tau = TRUE,
                      c_xi = 0.1,
                      c_tau = 0.1,
                      a_eq_c_xi = FALSE,
                      a_eq_c_tau = FALSE,
                      learn_kappa2_B = TRUE,
                      learn_lambda2_B = TRUE,
                      kappa2_B = 20,
                      lambda2_B = 20,
                      hyperprior_param,
                      sv_param,
                      MH_tuning,
                      phi_param,
                      display_progress = TRUE) {
  # Check if time-varying inputs are present
  tv_inputs <- "tvsurv" %in% class(data)

  assert(missing(group) || length(group) == nrow(data), "Grouping indicator, group, has to be omitted or the same length as data")


  ## Implementation of formula interface:
  # Check if formula is a formula
  if (inherits(formula, "formula") == FALSE) {
    stop("formula is not of class formula")
  }

  mf <- match.call(expand.dots = FALSE)
  m <- match(
    x = c("formula", "data", "subset"),
    table = names(mf), nomatch = 0L
  )
  mf <- mf[c(1L, m)]

  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())
  for (name in names(mf)) {
    assert(!any(class(mf[[name]]) %in% c("POSIXct", "POSIXt", "Date")), "No date variables allowed as predictors")
  }
  mr <- model.response(mf, "numeric")
  # Create Vector y
  if (tv_inputs) {
    y <- attr(data, "orig_response")
    delta <- attr(data, "orig_delta")

    if (inherits(mr, "Surv")) {
      warning("If time-varying inputs are provided, values for delta are extracted from data and the status indictor of the Surv object is ignored",
        immediate. = TRUE
      )
    }

    if (!missing(delta)) {
      warning("If time-varying inputs are provided, values for delta are extracted from data and input delta is ignored",
        immediate. = TRUE
      )
    }
  } else {
    if (inherits(mr, "Surv")) {
      if (!missing(delta)) {
        warning("If a Surv object is provided, values for delta are extracted from the Surv object and input delta is ignored",
          immediate. = TRUE
        )
      }
      y <- mr[, 1]
      delta <- mr[, 2]
    } else {
      y <- mr
      assert(length(delta) == nrow(data), "Status indicator, delta, has to be the same length as data")
    }
  }
  mt <- attr(x = mf, which = "terms")
  # Create Matrix X with dummies and transformations
  z <- model.matrix(object = mt, data = mf)

  colnames(z)[colnames(z) == "(Intercept)"] <- "Intercept"

  # Check that there are no NAs in y and x
  assert(!any(is.na(y)), "No NA values are allowed in survival time")
  assert(all(y > 0), "Survival times must be positive. Zeros are not allowed.")
  assert(!any(is.na(z)), "No NA values are allowed in covariates")

  assert(mod_type %in% c("double", "triple", "ridge"), "Allowed model types are: ridge, double, and triple")
  # default hyperparameter values
  default_hyper <- list(
    c0 = 2.5,
    g0 = 5,
    G0 = 5 / (2.5 - 1),
    e1 = 0.001,
    e2 = 0.001,
    d1 = 0.001,
    d2 = 0.001,
    beta_a_xi = 10,
    beta_a_tau = 10,
    alpha_a_xi = 5,
    alpha_a_tau = 5,
    beta_c_xi = 2,
    beta_c_tau = 2,
    alpha_c_xi = 5,
    alpha_c_tau = 5
  )

  if (!missing(group)) {
    group <- checkvalues(group)
    default_hyper$sigma2_phi <- rep(1, length(unique(group$values)))
  } else {
    group_sort <- c(0)
    default_hyper$sigma2_phi <- c(0)
  }

  # default sv params
  default_hyper_sv <- list(
    Bsigma_sv = 1,
    a0_sv = 5,
    b0_sv = 1.5
  )

  # default tuning parameters
  default_tuning_par <- list(
    a_xi_adaptive = TRUE,
    a_xi_tuning_par = 1,
    a_xi_target_rate = 0.44,
    a_xi_max_adapt = 0.01,
    a_xi_batch_size = 50,
    a_tau_adaptive = TRUE,
    a_tau_tuning_par = 1,
    a_tau_target_rate = 0.44,
    a_tau_max_adapt = 0.01,
    a_tau_batch_size = 50,
    c_xi_adaptive = TRUE,
    c_xi_tuning_par = 1,
    c_xi_target_rate = 0.44,
    c_xi_max_adapt = 0.01,
    c_xi_batch_size = 50,
    c_tau_adaptive = TRUE,
    c_tau_tuning_par = 1,
    c_tau_target_rate = 0.44,
    c_tau_max_adapt = 0.01,
    c_tau_batch_size = 50
  )

  ## Merge user supplied and default parameters
  if (missing(hyperprior_param)) {
    hyperprior_param <- default_hyper
  } else {
    hyperprior_param <- list_merger(default_hyper, hyperprior_param)
  }

  if (missing(sv_param)) {
    sv_param <- default_hyper_sv
  } else {
    sv_param <- list_merger(default_hyper_sv, sv_param)
  }

  if (missing(MH_tuning)) {
    MH_tuning <- default_tuning_par
  } else {
    MH_tuning <- list_merger(default_tuning_par, MH_tuning)
  }
  # Check if all numeric inputs are correct
  to_test_num <- list(
    lambda2_B = lambda2_B,
    kappa2_B = kappa2_B,
    a_xi = a_xi,
    a_tau = a_tau,
    c_xi = c_xi,
    c_tau = c_tau
  )
  if (!missing(group)) {
    cond <- length(hyperprior_param$sigma2_phi) == length(unique(group$values)) && all(!sapply(hyperprior_param$sigma2_phi, numeric_input_bad))
    assert(cond, "all elements of sigma2_phi have to be positive real numbers")
  }
  if (missing(hyperprior_param) == FALSE) {
    to_test_num <- c(to_test_num, hyperprior_param[names(hyperprior_param) != "sigma2_phi"])
  }

  if (missing(sv_param) == FALSE) {
    to_test_num <- c(to_test_num, sv_param)
  }

  if (missing(MH_tuning) == FALSE) {
    to_test_num <- c(to_test_num, MH_tuning[!grepl("(batch|adaptive)", names(MH_tuning))])
  }

  bad_inp <- sapply(to_test_num, numeric_input_bad)


  if (any(bad_inp)) {
    stand_names <- c(names(default_hyper), names(default_hyper_sv), "lambda2_B", "kappa2_B", "a_xi", "a_tau", "c_xi", "c_tau")
    bad_inp_names <- names(to_test_num)[bad_inp]
    bad_inp_names <- bad_inp_names[bad_inp_names %in% stand_names]
    stop(paste0(
      paste(bad_inp_names, collapse = ", "),
      ifelse(length(bad_inp_names) == 1, " has", " have"),
      " to be a real, positive number"
    ))
  }


  # Check the adapt_rates seperately
  if (any(0 > MH_tuning[grepl("rate", names(MH_tuning))] | MH_tuning[grepl("rate", names(MH_tuning))] > 1)) {
    stop("all target_rate parameters in MH_tuning have to be > 0 and < 1")
  }

  # Check if all integer inputs are correct
  to_test_int <- c(
    niter = niter,
    nburn = nburn,
    nthin = nthin,
    MH_tuning[grepl("batch", names(MH_tuning))]
  )

  bad_int_inp <- sapply(to_test_int, int_input_bad)

  if (any(bad_int_inp)) {
    bad_inp_names <- names(to_test_int)[bad_int_inp]
    stop(paste0(
      paste(bad_inp_names, collapse = ", "),
      ifelse(length(bad_inp_names) == 1, " has", " have"),
      " to be a single, positive integer"
    ))
  }

  if ((niter - nburn) < 2) {
    stop("niter has to be larger than or equal to nburn + 2")
  }

  if (nthin == 0) {
    stop("nthin can not be 0")
  }

  if ((niter - nburn) / 2 < nthin) {
    stop("nthin can not be larger than (niter - nburn)/2")
  }

  # Check if all boolean inputs are correct
  to_test_bool <- c(
    learn_lambda2_B = learn_lambda2_B,
    learn_kappa2_B = learn_kappa2_B,
    learn_a_xi = learn_a_xi,
    learn_a_tau = learn_a_tau,
    display_progress = display_progress,
    MH_tuning[grepl("adaptive", names(MH_tuning))]
  )

  bad_bool_inp <- sapply(to_test_bool, bool_input_bad)

  if (any(bad_bool_inp)) {
    bad_inp_names <- names(to_test_bool)[bad_bool_inp]
    stop(paste0(
      paste(bad_inp_names, collapse = ", "),
      ifelse(length(bad_inp_names) == 1, " has", " have"),
      " to be a single logical value"
    ))
  }


  default_phi <- list(
    mod_type_phi = "double",
    learn_a_phi = TRUE,
    a_phi = .1,
    learn_c_phi = TRUE,
    c_phi = .1,
    a_phi_eq_c_phi = FALSE,
    learn_lambda2_B_phi = TRUE,
    lambda2_B_phi = 20,
    e1_phi = 0.001,
    e2_phi = 0.001,
    beta_a_phi = 10,
    alpha_a_phi = 5,
    beta_c_phi = 10,
    alpha_c_phi = 5,
    a_phi_adaptive = TRUE,
    a_phi_tuning_par = 1,
    a_phi_target_rate = 0.44,
    a_phi_max_adapt = 0.01,
    a_phi_batch_size = 50,
    c_phi_adaptive = TRUE,
    c_phi_tuning_par = 1,
    c_phi_target_rate = 0.44,
    c_phi_max_adapt = 0.01,
    c_phi_batch_size = 50
  )

  if (missing(phi_param)) {
    phi_param <- default_phi
  } else {
    phi_param <- list_merger(default_phi, phi_param)
  }

  phi_param$target_rates_phi <- unlist(phi_param[grep("target", names(phi_param))])
  phi_param$max_adapts_phi <- unlist(phi_param[grep("max", names(phi_param))])
  phi_param$batch_sizes_phi <- unlist(phi_param[grep("size", names(phi_param))])
  phi_param$adaptive_phi <- unlist(phi_param[grep("adaptive", names(phi_param))])

  # Extract colnames

  d <- ncol(z)
  if (!is.null(colnames(z))) {
    col_names <- colnames(z)
  } else {
    col_names <- as.character(1:d)
  }

  # Sort objects
  order <- order(y, decreasing = TRUE)
  y_sort <- y[order]

  if (tv_inputs) {
    z_sort <- as.matrix(z)
  } else {
    z_sort <- z[order, ]
    z_sort <- as.matrix(z_sort)
  }
  if (!missing(group)) {
    # Group is guaranteed to be the values here
    group_sort <- group$values[order]
  }

  if (2 %in% delta) {
    delta <- delta - 1
  }
  assert(
    all(delta %in% c(0, 1)),
    "delta must contain only 0/1, 1/2, or TRUE/FALSE"
  )

  delta_sort <- as.matrix(as.integer(delta[order]))


  runtime <- system.time({
    res <- do_shrinkDSM(
      y_sort,
      z_sort,
      mod_type,
      delta_sort,
      S,
      group_sort,
      niter,
      nburn,
      nthin,
      hyperprior_param$d1,
      hyperprior_param$d2,
      hyperprior_param$e1,
      hyperprior_param$e2,
      hyperprior_param$sigma2_phi,
      learn_lambda2_B,
      learn_kappa2_B,
      lambda2_B,
      kappa2_B,
      learn_a_xi,
      learn_a_tau,
      a_xi,
      a_tau,
      learn_c_xi,
      learn_c_tau,
      c_xi,
      c_tau,
      a_eq_c_xi,
      a_eq_c_tau,
      MH_tuning$a_xi_tuning_par,
      MH_tuning$a_tau_tuning_par,
      MH_tuning$c_xi_tuning_par,
      MH_tuning$c_tau_tuning_par,
      hyperprior_param$beta_a_xi,
      hyperprior_param$beta_a_tau,
      hyperprior_param$alpha_a_xi,
      hyperprior_param$alpha_a_tau,
      hyperprior_param$beta_c_xi,
      hyperprior_param$beta_c_tau,
      hyperprior_param$alpha_c_xi,
      hyperprior_param$alpha_c_tau,
      sv_param$Bsigma_sv,
      sv_param$a0_sv,
      sv_param$b0_sv,
      display_progress,
      unlist(MH_tuning[grep("adaptive", names(MH_tuning))]),
      unlist(MH_tuning[grep("target", names(MH_tuning))]),
      unlist(MH_tuning[grep("max", names(MH_tuning))]),
      unlist(MH_tuning[grep("size", names(MH_tuning))]),
      tv_inputs,
      phi_param
    )
  })
  if (display_progress == TRUE) {
    cat("Timing (elapsed): ", file = stderr())
    cat(runtime["elapsed"], file = stderr())
    cat(" seconds.\n", file = stderr())
    cat(round((niter + nburn) / runtime[3]), "iterations per second.\n\n", file = stderr())
    cat("Converting results to coda objects and summarizing draws... ", file = stderr())
  }



  # Remove empty storage elements
  res[sapply(res, function(x) 0 %in% dim(x))] <- NULL
  res$MH_diag[sapply(res$MH_diag, function(x) 0 %in% dim(x))] <- NULL

  # Create object to hold prior values
  res$priorvals <- c(hyperprior_param,
    sv_param,
    a_xi = a_xi,
    a_tau = a_tau,
    lambda2_B = lambda2_B,
    kappa2_B = kappa2_B
  )

  # Add data to output
  res[["model"]] <- list()
  res$model$z <- z
  res$model$y <- y
  res$model$formula <- formula
  res$model$xlevels <- .getXlevels(mt, mf)
  res$model$terms <- mt
  res$model$model <- mf


  # add attributes to the individual objects if they are distributions or individual statistics
  nsave <- floor((niter - nburn) / nthin)
  for (i in names(res)) {
    attr(res[[i]], "type") <- ifelse(nsave %in% dim(res[[i]]), "sample", "stat")

    # Name each individual sample for plotting frontend
    if (attr(res[[i]], "type") == "sample") {
      if (i == "phi") {
        colnames(res[[i]]) <- paste0(i, unique(group_sort))
      } else if (dim(res[[i]])[2] == d) {
        colnames(res[[i]]) <- paste0(i, "_", col_names)
      } else if (dim(res[[i]])[2] == 2 * d) {
        colnames(res[[i]]) <- paste0(i, "_", rep(col_names, 2))
      } else {
        colnames(res[[i]]) <- i
      }
    }

    # Change objects to be coda compatible
    # Only apply to posterior samples
    if (attr(res[[i]], "type") == "sample") {
      # Differentiate between TVP and non TVP
      if (is.na(dim(res[[i]])[3]) == FALSE) {
        # Create a sub list containing an mcmc object for each parameter in TVP case
        dat <- res[[i]]
        res[[i]] <- list()
        for (j in 1:dim(dat)[2]) {
          res[[i]][[j]] <- as.mcmc(t(dat[, j, ]), start = niter - nburn, end = niter, thin = nthin)
          colnames(res[[i]][[j]]) <- paste0(i, "_", j, "_", 1:ncol(res[[i]][[j]]))

          # make it of class mcmc.tvp for custom plotting function
          class(res[[i]][[j]]) <- c("mcmc.dsm.tvp", "mcmc")
          attr(res[[i]][[j]], "S") <- S # c(S, max(res$model$y))
          attr(res[[i]][[j]], "lastsurvtime") <- max(res$model$y)
          attr(res[[i]][[j]], "type") <- "sample"
        }

        if (length(res[[i]]) == 1) {
          res[[i]] <- res[[i]][[j]]
        }

        # Make it of type 'sample' again
        attr(res[[i]], "type") <- "sample"

        # Rename
        if (dim(dat)[2] > 1) {
          names(res[[i]]) <- colnames(dat)
        }
      } else {
        res[[i]] <- as.mcmc(res[[i]], start = niter - nburn, end = niter, thin = nthin)
      }
    }

    # Create summary of posterior
    if (is.list(res[[i]]) == FALSE & attr(res[[i]], "type") == "sample") {
      if (i != "theta_sr" & i != "beta") {
        res$summaries[[i]] <- t(apply(res[[i]], 2, function(x) {
          obj <- as.mcmc(x, start = niter - nburn, end = niter, thin = nthin)
          ESS <- tryCatch(coda::effectiveSize(obj),
            error = function(err) {
              warning("Calculation of effective sample size failed for one or more variable(s). This can happen if the prior placed on the model induces extreme shrinkage.")
              return(NA)
            }, silent = TRUE
          )

          return(c(
            "mean" = mean(obj),
            "sd" = sd(obj),
            "median" = median(obj),
            "HPD" = HPDinterval(obj)[c(1, 2)],
            "ESS" = round(ESS)
          ))
        }))
      } else if (i == "theta_sr") {
        res$summaries[[i]] <- t(apply(res[[i]], 2, function(x) {
          obj <- as.mcmc(abs(x), start = niter - nburn, end = niter, thin = nthin)
          ESS <- tryCatch(coda::effectiveSize(obj),
            error = function(err) {
              warning("Calculation of effective sample size failed for one or more variable(s). This can happen if the prior placed on the model induces extreme shrinkage.")
              return(NA)
            }, silent = TRUE
          )

          return(c(
            "mean" = mean(obj),
            "sd" = sd(obj),
            "median" = median(obj),
            "HPD" = HPDinterval(obj)[c(1, 2)],
            "ESS" = round(ESS)
          ))
        }))
      }
    }
  }

  if (!missing(group)) {
    # Identify factor loadings
    modif <- ifelse(res$phi[, 1] < 0, -1, 1)
    res$phi <- modif * res$phi
    res$f <- modif * res$f

    # add original names to groups
    key <- as.numeric(gsub("phi", "", colnames(res$phi)))
    orig_names <- paste0("phi_", names(group$levels)[key + 1])
    colnames(res$phi) <- orig_names
  }

  if (display_progress == TRUE) {
    cat("Done!\n", file = stderr())
  }

  # add some attributes for the methods and plotting
  attr(res, "class") <- "shrinkDSM"
  attr(res, "S") <- S # c(S, max(res$model$y))
  attr(res, "lastsurvtime") <- max(res$model$y)
  if (!missing(group)) {
    attr(res, "group") <- group
  }
  attr(res, "learn_a_xi") <- learn_a_xi
  attr(res, "learn_a_tau") <- learn_a_tau
  attr(res, "learn_kappa2_B") <- learn_kappa2_B
  attr(res, "learn_lambda2_B") <- learn_lambda2_B
  attr(res, "niter") <- niter
  attr(res, "nburn") <- nburn
  attr(res, "nthin") <- nthin
  attr(res, "colnames") <- col_names
  attr(res, "tv_inputs") <- tv_inputs

  return(res)
}
