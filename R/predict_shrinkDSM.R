
create_lambda <- function(obj, x_test, group = 0){
  # Handle special case where only intercept is estimated
  nrows <- ifelse(length(attr(obj, "colnames")) == 1, nrow(obj$beta), nrow(obj$beta[[1]]))
  ncols <- ifelse(length(attr(obj, "colnames")) == 1, ncol(obj$beta), ncol(obj$beta[[1]]))
  res <- matrix(0, nrow = nrows, ncol = ncols)

  # Again, handle special case where only intercept is estimated
  if (length(attr(obj, "colnames")) == 1) {
    res <- obj$beta * x_test[attr(obj, "colnames")]
  } else {
    lapply(attr(obj, "colnames"), function(nam){
      parent <- parent.frame(2)
      parent$res <- parent$res + obj$beta[[paste0("beta_", nam)]] * x_test[nam]
    })
  }

  if (!is.null(attr(obj, "group"))) {
    res[ ,2:ncol(res)] <- res[ ,2:ncol(res)] + obj$f * obj$phi[, group + 1]
  }
  exp(res)
}

predict_obs <- function(x_row, obj, ndraws, cens) {
  if (!is.null(attr(obj, "group"))) {
    group <- x_row[length(x_row)]
  } else {
    group <- NA
  }
  lams <- create_lambda(obj, x_row, group)
  time_points <- c(0, attr(obj, "S"))

  if (max(obj$model$y) != max(attr(obj, "S"))) {
    time_points <- c(time_points, max(obj$model$y))
  }

  U <- diff(time_points)

  time_per <- 1:length(U)
  par_surv <- matrix(0, nrow = ndraws, ncol = length(U))
  selector <- rep(TRUE, ndraws)
  for (j in time_per){
    if (sum(selector) == 0) break

    par_surv[selector, j] <- rexp(sum(selector), lams[selector, j + 1])
    selector <- par_surv[, j] > U[j]
    if (j < length(U)) {
      par_surv[selector, j] <- U[j]
    } else if (cens == TRUE) {
      par_surv[selector, j] <- U[j]
    }
  }

  y <- apply(par_surv, 1, sum)
  return(y)
}



#' Draw from posterior predictive density of a fitted time-varying parameter survival model
#'
#' Draws from the posterior predictive distribution of survival times based on a fitted
#' time-varying parameter survival model resulting from a call to
#' \code{shrinkDSM}.
#'
#' @param object an object of class \code{shrinkDSM}, containing the fitted model.
#' @param newdata a data frame containing the covariates used for the prediction.
#' The names of the covariates
#' have to match the names used during model estimation in the call to \code{shrinkDSM}.
#' @param cens logical value indicating whether the predictions should be censored at the
#' largest survival time in the data used for estimation.
#' The default value is \code{TRUE}.
#' @param ... included for S3 method consistency and currently ignored.
#'
#' @return The value returned is a list object of class \code{shrinkTVP_pred} containing the
#' samples from the
#' posterior predictive density.
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @author Daniel Winkler \email{dwinkler@@wu.ac.at}
#' @family prediction functions
#'
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
#' # Draw from posterior predictive distribution
#' newdata <- data.frame(radiation = c(0, 1))
#' pred <- predict(mod, newdata = newdata)
#' }
#'@export
predict.shrinkDSM <- function(object, newdata, cens = TRUE, ...){

  assert(!attr(object, "tv_inputs"), "Predictions are not possible for models with time-varying covariates.")

  terms <- delete.response(object$model$terms)
  m <- model.frame(terms, data = newdata, xlev = object$model$xlevels)
  x_test <- model.matrix(terms, m)
  if(!is.null(attr(object, "group"))) {
    x_test <- cbind(x_test, newdata[, attr(object, "group")])
  }

  colnames(x_test)[colnames(x_test) == "(Intercept)"] <- "Intercept"

  # Handle special case where only an intercept is estimated
  ndraws <- ifelse(length(attr(object, "colnames")) == 1, nrow(object$beta), nrow(object$beta[[1]]))

  res <- apply(x_test, 1, function(x) predict_obs(x, object, ndraws, cens))
  per_surv <- apply(res, 2, function(y) sum(y >= max(object$model$y))/length(y))

  class(res) <- "shrinkDSM_pred"
  attr(res, "censored") <- cens
  attr(res, "xlim") <- c(min(object$model$y), max(object$model$y))
  attr(res, "per_surv") <- per_surv
  return(res)
}

#' Graphical summary of posterior predictive density
#'
#' \code{plot.shrinkDSM_pred} generates plots visualizing the posterior predictive density generated by \code{predict.shrinkDSM}.
#'
#' @param x a \code{shrinkDSM_pred} object.
#' @param dens_args \emph{optional} named list containing arguments passed to \code{density}.
#' @param legend logical value inidicating whether a legend should be added. Defaults to \code{TRUE}.
#' @param ... further arguments to be passed to \code{plots}.
#'
#' @return Called for its side effects and returns invisibly.
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
#' # Draw from posterior predictive distribution
#' newdata <- data.frame(radiation = c(0, 1))
#' pred <- predict(mod, newdata = newdata)
#'
#' # Plot predictions
#' plot(pred)
#' }
#'@export
plot.shrinkDSM_pred <- function(x, dens_args, legend = TRUE, ...){
  assert(!bool_input_bad(legend), "legend has to be a single logical value")
  # Fuse user and standard args for density
  xlim <- attr(x, "xlim")
  standard_plot_args <- list(from = xlim[1], to = xlim[2])
  if(!missing(dens_args)){
    assert(is.list(dens_args), "dens_args has to be a list")
    missing_args <- names(standard_plot_args)[!names(standard_plot_args) %in% names(dens_args)]
    dens_args[missing_args] <- standard_plot_args[missing_args]
  } else{
    dens_args <- standard_plot_args
  }
  dens_args$x <- x[,1]

  dens_estim <- list()
  dens_estim[[1]] <- do.call(density, dens_args)

  if (ncol(x) > 1){
    for (i in 2:ncol(x)){
      dens_args$x <- x[,i]
      dens_estim[[i]] <- do.call(density, dens_args)
    }
  }
  # Fuse user and standard args for plot
  ## Extract all user specified args
  plot_args <- list(...)
  standard_plot_args <- list(
    ylim = c(0, max(sapply(dens_estim, function(x) max(x$y * 1.1)))),
    main = "Predictive Density of Survival Time"
   )
  missing_args <- names(standard_plot_args)[!names(standard_plot_args) %in% names(plot_args)]
  plot_args[missing_args] <- standard_plot_args[missing_args]
  plot_args$x <- dens_estim[[1]]
  do.call(plot, plot_args)

  if (ncol(x) > 1){
    for (i in 2:ncol(x)){
      lines(dens_estim[[i]], col = i)
    }
  }

  if(legend){
  lty <- ifelse("lty" %in% names(plot_args), plot_args$lty, 1)
  legend("topright",
    legend = paste("Percent survived:", round(attr(x, "per_surv")*100, 2)),
    col = 1:ncol(x), lty = lty, bty = "n")
  }
}

