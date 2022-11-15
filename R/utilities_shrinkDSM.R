#' Create division points for estimation of a dynamic survival model
#'
#' Create a vector of division points for the model.
#' These points mark the times at which the parameters are allowed to evolve,
#' with the parameters being fixed between division points. The points
#' are generated in a data driven fashion, with a new point being created
#' when \code{events} number of interesting events have been observed since
#' the last division point.
#'
#' @param times A vector of real, positive numbers indicating
#' the survival times. For right censored data, this is the follow up time.
#' @param delta A vector of status indicators, with 0 = alive and 1 = dead.
#' Other choices are \code{TRUE}/\code{FALSE} (\code{TRUE} = death) or 1/2 (2 = death).
#' @param events A positive integer indicating the number of interesting events per interval until
#' a new division is created.
#'
#' @return Returns an integer vector of time points to be used as division points \code{S} in \code{shrinkDSM}.
#'
#' @author Daniel Winkler \email{daniel.winkler@@wu.ac.at}
#' @examples
#' data("gastric")
#'
#' # Create intervals for piecewise exponential model
#' intervals <- divisionpoints(gastric$time, gastric$status, 2)
#' @export
divisionpoints <- function(times, delta, events = 1){
  assert(is.vector(times), "times has to be a vector")
  assert(is.vector(delta), "delta has to be a vector")
  assert(length(times) == length(delta), "times and delta have to be of same length")
  assert(all(times > 0), "all elements of times must be positive.")
  if(2 %in% delta){delta <- delta - 1}
  assert(all(delta %in% c(0,1)),
         "delta must contain only 0/1, 1/2, or TRUE/FALSE")
  assert((events %% 1) == 0, "events has to be an integer")
  assert(events >= 1, "At least one observation per interval is needed (events > 0)")


  S <- sort(unique(times[delta == 1]))[c(rep(FALSE, events-1), TRUE)]
  return(S)
}

summary_tmp <- utils::getFromNamespace("summary.shrinkTVP", "shrinkTVP")

#' @export
summary.shrinkDSM <- function(object, digits = 3, ...) {
  summary_tmp(object, digits, FALSE)
}

# Small convenience function to check if something is a scalar
is.scalar <- function(x) is.atomic(x) && length(x) == 1


plot_tmp <- utils::getFromNamespace("plot.mcmc.tvp", "shrinkTVP")


#' Graphical summary of posterior distribution for a piecewise constant, time-varying parameter
#'
#' \code{plot.mcmc.dsm.tvp} plots empirical posterior quantiles for a piecewise constant, time-varying parameter.
#'
#' @param x \code{mcmc.dsm.tvp} object
#' @param probs vector of boundaries for credible intervals to plot for each point in time, with values in [0,1].
#' The largest and smallest value form the outermost credible interval, the second smallest and second largest the second outermost and so forth.
#' The default value is \code{c(0.025, 0.25, 0.75, 0.975)}. Note that there have to be the same number of probs
#' < 0.5 as there are > 0.5.
#' @param shaded single logical value or a vector of logical values, indicating whether or not to shade the area between the pointwise
#' credible intervals. If a vector is given, the first value given is used to determine if the area between the outermost credible interval
#' is shaded, the second for the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the
#' number of quantile pairs. The default value is \code{TRUE}.
#' @param quantlines single logical value or a vector of logical values, indicating whether or not to draw borders along the pointwise
#' credible intervals. If a vector is given, the first value given is used to determine whether the outermost credible interval
#' is marked by lines, the second for the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the
#' number of credible intervals. The default value is \code{FALSE}.
#' @param shadecol single character string or a vector of character strings. Determines the color of the shaded areas that represent
#' the credible intervals. If a vector is given, the first color given is used for the outermost area,
#' the second for the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the number
#' of shaded areas. Has no effect if \code{shaded = FALSE}. The default value is \code{"skyblue"}.
#' @param shadealpha real number between 0 and 1 or a vector of real numbers between 0 and 1.
#' Determines the level of transparency of the shaded areas that represent
#' the credible intervals. If a vector is given, the first value
#' given is used for the outermost area, the second for the second outermost and so forth. Recycled in the usual fashion
#' if the vector is shorter than the number of shaded areas. Has no effect if \code{shaded = FALSE}.
#' The default value is \code{0.5}.
#' @param quantlty either a single integer in [0,6] or one of the character strings \code{"blank",
#' "solid", "dashed", "dotted", "dotdash", "longdash", "twodash"} or a vector containing these. Determines the line type of the borders
#' drawn around the shaded areas that represent the credible intervals. Note that if a vector is supplied the elements have to either
#' be all integers or all character strings. If a vector is given, the first value given is used for the outermost area, the second for
#' the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the number of shaded areas.
#' Has no effect if \code{quantlines = FALSE}. The default value is \code{2}.
#' @param quantcol single character string or a vector of character strings. Determines the color of the borders drawn around the shaded
#' areas that represent the credible intervals. If a vector is given, the first color given is used for borders of the outermost area,
#' the second for the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the number
#' of shaded areas. Has no effect if \code{quantlines = FALSE}. The default value is \code{"black"}.
#' @param quantlwd single real, positive number or a vector of real, positive numbers. Determines the line width of the borders
#' drawn around the shaded areas that represent the credible intervals. If a vector is given, the first number given is used for
#' the borders of the outermost area, the second for the second outermost and so forth. Recycled in the usual fashion if the vector
#' is shorter than the number of shaded areas. Has no effect if \code{quantlines = FALSE}. The default value is \code{0.5}.
#' @param drawzero single logical value determining whether to draw a horizontal line at zero or not. The default value is \code{TRUE}.
#' @param zerolty single integer in [0,6] or one of the character strings \code{"blank",
#' "solid", "dashed", "dotted", "dotdash", "longdash", "twodash"}. Determines the line type of the horizontal line at zero. Has no effect
#' if \code{drawzero = FALSE}. The default value is \code{2}.
#' @param zerolwd single real, positive number. Determines the line width of the horizontal line at zero. Has no effect
#' if \code{drawzero = FALSE}. The default value is \code{1}.
#' @param zerocol single character string. Determines the color of the horizontal line at zero. Has no effect if \code{drawzero = FALSE}.
#' The default value is \code{"grey"}.
#' @param ... further arguments to be passed to \code{plot}.
#' @return Called for its side effects and returns invisibly.
#'
#' @method plot mcmc.dsm.tvp
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family plotting functions
#' @examples
#' \donttest{
#' set.seed(123)
#' data("gastric")
#'
#' # Create intervals for piecewise exponential model
#' intervals <- divisionpoints(gastric$time, gastric$status, 2)
#'
#' # Estimate model
#' mod <- shrinkDSM(time ~ radiation, gastric,
#'                  delta = gastric$status, S = intervals)
#'
#' # Plot piecewise constant, time-varying parameter
#' plot(mod$beta$beta_radiation)
#' }
#' @export
plot.mcmc.dsm.tvp <- function(x, probs = c(0.025, 0.25, 0.75, 0.975),
                              shaded = TRUE, quantlines = FALSE,
                              shadecol = "skyblue", shadealpha = 0.5,
                              quantlty = 2, quantcol = "black", quantlwd = 0.5,
                              drawzero = TRUE, zerolty = 2, zerolwd = 1, zerocol = "grey", ...){
  S <- attr(x, "S")
  maxy <- attr(x, "lastsurvtime")

  if (S[length(S)] == maxy) S <- S[1:(length(S) - 1)]
  if (S[1] == 0) S <- S[2:length(S)]

  if(ncol(x) > (length(S) + 1)) {
    x_internal <- x[, 2:ncol(x)]
  } else {
    x_internal <- x
  }

  x_vec <- c(0, rep(S, each = 2), maxy)
  x_list <- rep(split(x_internal, rep(1:(ncol(x_internal)), each = nrow(x_internal))), each = 2)
  newobj <- do.call(cbind, x_list)
  attr(newobj, "index") <- x_vec

  plot_tmp(newobj, probs,
           shaded, quantlines,
           shadecol, shadealpha,
           quantlty, quantcol, quantlwd,
           drawzero, zerolty, zerolwd,  zerocol, ...)
}

plot_temp_glob <- utils::getFromNamespace("plot.shrinkTVP", "shrinkTVP")

#' Graphical summary of posterior distribution of fitted dynamic survival model
#'
#' \code{plot.shrinkDSM} generates plots visualizing the posterior distribution estimated
#' as a result from a call to \code{shrinkDSM}.
#'
#' @param x a \code{shrinkDSM} object.
#' @param pars a character vector containing the names of the parameters to be visualized.
#' The names have to coincide with the names of the list elements of the \code{shrinkDSM}
#' object. Throws an error if any element of \code{pars} does not fulfill this criterium.
#' The default is \code{c("beta")}.
#' @param nplot positive integer that indicates the number of tvp plots to display on a single
#' page before a new page is generated. The default value is 3.
#' @param h_borders single real, positive number smaller than 0.5 or a vector containing two such numbers. Determines
#' the relative amount of space (the total amount summing up to 1) left blank on the left and right of the plot, in that order.
#' The default is \code{c(0.05, 0.05)}.
#' @param w_borders single real, positive number smaller than 0.5 or a vector containing two such numbers. Determines
#' the relative amount of space (the total amount summing up to 1) left blank at the top and bottom of the plot, in that order.
#' The default is \code{c(0.02, 0.02)}.
#' @param ... further arguments to be passed to the respective plotting functions.
#' @return Called for its side effects and returns invisibly.
#' @examples
#' \donttest{
#' set.seed(123)
#' data("gastric")
#'
#' # Create intervals for piecewise exponential model
#' intervals <- divisionpoints(gastric$time, gastric$status, 2)
#'
#' # Estimate model
#' mod <- shrinkDSM(time ~ radiation, gastric,
#'                  delta = gastric$status, S = intervals)
#' plot(mod)
#' }
#'
#' # Will produce an error because 'hello' is not a parameter in the model
#' \dontrun{
#' plot(mod, pars = c("beta", "hello"))
#' }
#'
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family plotting functions
#' @export
plot.shrinkDSM <- function(x, pars = c("beta"), nplot = 3, h_borders = c(0.05, 0.05), w_borders = c(0.02, 0.02), ...){
  plot_temp_glob(x, pars, nplot, h_borders, w_borders, ...)
}

## Get functions from shrinkTVP
list_merger <- utils::getFromNamespace("list_merger", "shrinkTVP")
numeric_input_bad <- utils::getFromNamespace("numeric_input_bad", "shrinkTVP")
int_input_bad <- utils::getFromNamespace("int_input_bad", "shrinkTVP")
bool_input_bad <- utils::getFromNamespace("bool_input_bad", "shrinkTVP")

## Create valid grouping variable for calculations
## saving mapping to original names for factors
checkvalues <- function(values){
  values <- as.factor(values)
  vlevels <- levels(values)
  values <- as.numeric(values)-1
  level_map <- list()
  level_map <- mapply(function(g,v) {level_map[[g]] = v},
                      vlevels, as.list((1:length(vlevels))-1))
  return (list(values = values, levels = level_map))
}

## Helper to check for errors and give better messages
assert <- function(cond, msg = ""){
  if (!cond){
    stop(msg,  call. = FALSE)
  }
}


#' Nicer printing of shrinkDSM objects
#'
#' @param x a \code{shrinkDSM} object.
#' @param ... Currently ignored.
#'
#' @return Called for its side effects and returns invisibly.
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @examples
#' \donttest{
#' set.seed(123)
#' data("gastric")
#'
#' # Create intervals for piecewise exponential model
#' intervals <- divisionpoints(gastric$time, gastric$status, 2)
#'
#' # Estimate model
#' mod <- shrinkDSM(time ~ radiation, gastric,
#'                  delta = gastric$status, S = intervals)
#'
#' # Print
#' print(mod)
#'
#' # Alternatively
#' mod
#' }
#' @export
print.shrinkDSM <- function(x, ...){
  ind <- attr(x, "index")
  begin <- min(x$model$y)
  end <- max(x$model$y)
  cat(paste0("Object containing a fitted DSM model with:\n",
             " - ", formatC(length(attr(x, "colnames")), width = 7), " covariates\n",
             " - ", formatC(length(x$model$y), width = 7), " intervals, encompassing time points from ", begin, " to ", end, "\n",
             " - ", formatC(attr(x, "niter"), width = 7), " MCMC draws\n",
             " - ", formatC(attr(x, "nburn"), width = 7), " burn-in\n",
             " - ", formatC(attr(x, "nthin"), width = 7), " thinning\n"))
  invisible(x)
}


#' Prepare time-varying inputs for estimation of a dynamic survival model
#'
#' This function pre-processes time-varying inputs in such a way that \code{shrinkDSM}
#' can work with time-varying inputs. Its main inputs are two data frames, namely
#' \code{surv_data} and \code{covariate_data}. \code{surv_data} contains meta data about
#' each observation (i.e. survival time and censoring indicator), while \code{covariate_data}
#' contains the time-varying covariates (one per observation and time interval) and an
#' index for which time interval each covariate is observed in. The two are merged
#' together via an ID that needs to be unique for each observation and present in both
#' \code{surv_data} and \code{covariate_data}.
#'
#' @param surv_data data frame containing meta-data for each observation (survival time and
#' censoring indicator) as well as an ID for each observation.
#' @param covariate_data data frame containing the time-varying covariates
#' (one per observation and time interval), an index for which time interval each covariate
#' is observed in as well as an ID for each observation.
#' @param id_var character string specifying the column name of the ID variable. If the name
#' is different in \code{surv_data} and \code{covariate_data}, \code{id_var} will be used in
#' \code{surv_data}, whereas \code{covariate_id_var} will be used in \code{covariate_data}.
#' @param surv_var character string specifying the column name of the survival times in
#' \code{surv_data}.
#' @param delta_var character string specifying the column name of the status indicators in
#' \code{surv_data}, 0 = censored or 1 = event observed..
#' @param interval_var character string specifying the column name of the time interval ID in
#' \code{covariate_data}.
#' @param covariate_id_var character string specifying the column name of the ID variable in
#' \code{covariate_data}. Defaults to \code{id_var}.
#'
#' @return Returns an object of class \code{data.frame} and \code{tvsurv} to be used as an input in
#' \code{shrinkDSM}.
#'
#' @author Daniel Winkler \email{daniel.winkler@@wu.ac.at}
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @examples
#' # A toy example with 5 observations and 2 covariates, observed over 3 time periods
#' set.seed(123)
#' n_obs <- 5
#' surv_var <- round(rgamma(n_obs, 1, .1)) + 1
#' delta_var <- sample(size = n_obs, c(0, 1), prob = c(0.2, 0.8), replace = TRUE)
#'
#' surv_data <- data.frame(id_var = 1:n_obs, surv_var, delta_var)
#'
#' # Determine intervals
#' S <- c(3, 11)
#'
#' # Create synthetic observations for each individual
#' covariate_list <- list()
#'
#' for (i in 1:n_obs) {
#'   nr_periods_survived <- sum(surv_var[i] > S) + 1
#'   covariate_list[[i]] <- data.frame(id_var = i,
#'                                     interval_var = 1:nr_periods_survived,
#'                                     x1 = rnorm(nr_periods_survived),
#'                                     x2 = rnorm(nr_periods_survived))
#' }
#'
#' # Bind all individual covariate data frames together
#' # Each observation now has a covariate in each period they
#' # were observed in.
#' covariate_data <- do.call(rbind, covariate_list)
#'
#' # Call prep_tvinput to pre-process for shrinkDSM
#' merged_data <- prep_tvinput(surv_data,
#'                             covariate_data,
#'                             id_var = "id_var",
#'                             surv_var = "surv_var",
#'                             delta_var = "delta_var",
#'                             interval_var = "interval_var")
#'
#' # Can now be used in shrinkDSM
#' # Note that delta is now automatically extracted from merged_data,
#' # providing it will throw a warning
#' mod <- shrinkDSM(surv_var ~ x1 + x2, merged_data, S = S)
#' @export
prep_tvinput <- function(surv_data, covariate_data, id_var, surv_var, delta_var, interval_var, covariate_id_var = id_var) {

  # Input checking
  assert((("data.frame" %in% class(surv_data)) && ("data.frame" %in% class(covariate_data))),
         "both surv_data and covariate_data have to be of class data.frame")

  for (i in c("id_var", "surv_var", "delta_var", "interval_var", "covariate_id_var")) {
    assert(class(get(i)) == "character" & is.scalar(get(i)), paste0(i, " has to be a single character string"))
  }

  assert(id_var %in% names(surv_data), "id_var column not found in surv_data")
  assert(surv_var %in% names(surv_data), "surv_var column not found in surv_data")
  assert(delta_var %in% names(surv_data), "delta_var column not found in surv_data")
  assert(interval_var %in% names(covariate_data), "interval_var column not found in surv_data")
  assert(covariate_id_var %in% names(covariate_data), "covariate_id_var column not found in surv_data")


  # This checks that all IDs in covariate_data are also in surv_data
  unmatched_ids <- unique(covariate_data[!covariate_data[, covariate_id_var] %in% surv_data[, id_var],
                                         covariate_id_var])
  if (length(unmatched_ids) > 0) {
    err <- paste0("following IDs found in covariate_data but not in surv_data: ", paste(unmatched_ids, collapse = ", "))
    stop(err)
  }


  # This piece of code checks whether any individual is missing certain intervals
  index <- unique(covariate_data[, interval_var])
  missing_intervals <- which(tapply(covariate_data[, interval_var], covariate_data[, covariate_id_var], function(x) {
    sorted_intervals <- sort(x)
    !sum(x %in% index) == (sum(sorted_intervals == index[1:length(sorted_intervals)]))
  }))

  if (length(missing_intervals > 0)) {
    err <- paste0("following IDs are missing intervals before last observation: ", paste(missing_intervals, collapse = ", "))
    stop(err)
  }


  # Get time period each observation exits sample
  dies_in <- c(tapply(covariate_data[,covariate_id_var], covariate_data[,covariate_id_var], length))
  surv_data$dies_in <- dies_in

  # Merge time varying (Z) and constant data (survival time, delta)
  # Creates a long data.frame with constant data repeated
  data_long <- merge(surv_data, covariate_data, by.x = id_var,
                     by.y = covariate_id_var, all = TRUE)

  # Sort data first by interval (in increasing order),
  # then by those that exit sample in current time interval,
  # then by overall survival time.
  data_long_order <- order(data_long[,interval_var],
                           data_long$dies_in == data_long[,interval_var],
                           data_long[,surv_var],
                           decreasing = c(FALSE, TRUE, TRUE))
  data_long_sorted <- data_long[data_long_order,]

  # Remove helper column from return
  data_long_sorted$dies_in <- NULL

  # Add class to differentiate in shrinkDSM
  attr(data_long_sorted, "class") <- c("data.frame", "tvsurv")

  # Add original response and delta for init.cpp
  attr(data_long_sorted, "orig_response") <- surv_data[, surv_var]
  attr(data_long_sorted, "orig_delta") <- surv_data[, delta_var]

  return(data_long_sorted)
}
