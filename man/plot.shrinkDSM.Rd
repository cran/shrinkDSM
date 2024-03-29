% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities_shrinkDSM.R
\name{plot.shrinkDSM}
\alias{plot.shrinkDSM}
\title{Graphical summary of posterior distribution of fitted dynamic survival model}
\usage{
\method{plot}{shrinkDSM}(
  x,
  pars = c("beta"),
  nplot = 3,
  h_borders = c(0.05, 0.05),
  w_borders = c(0.02, 0.02),
  ...
)
}
\arguments{
\item{x}{a \code{shrinkDSM} object.}

\item{pars}{a character vector containing the names of the parameters to be visualized.
The names have to coincide with the names of the list elements of the \code{shrinkDSM}
object. Throws an error if any element of \code{pars} does not fulfill this criterium.
The default is \code{c("beta")}.}

\item{nplot}{positive integer that indicates the number of tvp plots to display on a single
page before a new page is generated. The default value is 3.}

\item{h_borders}{single real, positive number smaller than 0.5 or a vector containing two such numbers. Determines
the relative amount of space (the total amount summing up to 1) left blank on the left and right of the plot, in that order.
The default is \code{c(0.05, 0.05)}.}

\item{w_borders}{single real, positive number smaller than 0.5 or a vector containing two such numbers. Determines
the relative amount of space (the total amount summing up to 1) left blank at the top and bottom of the plot, in that order.
The default is \code{c(0.02, 0.02)}.}

\item{...}{further arguments to be passed to the respective plotting functions.}
}
\value{
Called for its side effects and returns invisibly.
}
\description{
\code{plot.shrinkDSM} generates plots visualizing the posterior distribution estimated
as a result from a call to \code{shrinkDSM}.
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
plot(mod)
}

# Will produce an error because 'hello' is not a parameter in the model
\dontrun{
plot(mod, pars = c("beta", "hello"))
}

}
\seealso{
Other plotting functions: 
\code{\link{plot.mcmc.dsm.tvp}()},
\code{\link{plot.shrinkDSM_pred}()}
}
\author{
Peter Knaus \email{peter.knaus@wu.ac.at}
}
\concept{plotting functions}
