% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_shrinkDSM.R
\name{plot.shrinkDSM_pred}
\alias{plot.shrinkDSM_pred}
\title{Graphical summary of posterior predictive density}
\usage{
\method{plot}{shrinkDSM_pred}(x, dens_args, legend = TRUE, ...)
}
\arguments{
\item{x}{a \code{shrinkDSM_pred} object.}

\item{dens_args}{\emph{optional} named list containing arguments passed to \code{density}.}

\item{legend}{logical value inidicating whether a legend should be added. Defaults to \code{TRUE}.}

\item{...}{further arguments to be passed to \code{plots}.}
}
\value{
Called for its side effects and returns invisibly.
}
\description{
\code{plot.shrinkDSM_pred} generates plots visualizing the posterior predictive density generated by \code{predict.shrinkDSM}.
}
\examples{
\donttest{
set.seed(123)
data("gastric")

# Create intervals for piecewise exponential model
intervals <- divisionpoints(gastric$time, gastric$status, 2)

# Estimate model
mod <- shrinkDSM(time ~ radiation, gastric,
                 delta = gastric$status, S = intervals)

# Draw from posterior predictive distribution
newdata <- data.frame(radiation = c(0, 1))
pred <- predict(mod, newdata = newdata)

# Plot predictions
plot(pred)
}
}
\seealso{
Other plotting functions: 
\code{\link{plot.mcmc.dsm.tvp}()},
\code{\link{plot.shrinkDSM}()}
}
\concept{plotting functions}
