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


  S <- sort(unique(times[delta == 1])[c(rep(FALSE, events-1), TRUE)])
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
  if(ncol(x) > length(S) + 1) {x <- x[, 2:ncol(x)]}
  x_vec <- c(0, rep(S, each = 2), maxy)
  x_list <- rep(split(x, rep(1:(ncol(x)), each = nrow(x))), each = 2)
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
