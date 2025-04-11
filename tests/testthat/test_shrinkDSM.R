library(shrinkDSM)
library(survival)
set.seed(42)
data("gastric")
intervals <- divisionpoints(gastric$time, gastric$status, 2)
test_that("divisionpoints", {
    expect_equal(divisionpoints(gastric$time, gastric$status, 1), unique(gastric$time[gastric$status == 1]))
    expect_equal(intervals, unique(gastric$time[gastric$status == 1])[c(FALSE, TRUE)])
    expect_equal(divisionpoints(gastric$time, gastric$status, 3), unique(gastric$time[gastric$status == 1])[c(FALSE, FALSE, TRUE)])
})

mod <- shrinkDSM(time ~ radiation, gastric,
    delta = gastric$status, S = intervals, niter = 5
)

mod_surv <- shrinkDSM(Surv(time, status) ~ radiation, gastric,
    S = intervals, niter = 5
)
mod2 <- shrinkDSM(time ~ radiation, gastric,
    delta = gastric$status, S = intervals,
    mod_type = "triple", niter = 5
)
# Change some of the hyperparameters
mod3 <- shrinkDSM(time ~ radiation, gastric,
    delta = gastric$status, S = intervals,
    mod_type = "triple",
    hyperprior_param = list(
        beta_a_xi = 5,
        alpha_a_xi = 10
    ), niter = 5
)

mods <- list(mod, mod_surv, mod2, mod3)
# Estimate baseline model
test_that("Main shrinkDSM function", {
    for (idx in seq_along(mods)) {
        m <- mods[[idx]]
        expect_s3_class(m, "shrinkDSM")
        expect_visible(m)
    }
})

test_that("print", {
    for (idx in seq_along(mods)) {
        m <- mods[[idx]]
        expect_match(capture.output(print(m))[1], "Object containing a fitted DSM model with:")
    }
})

test_that("Summary", {
    for (idx in seq_along(mods)) {
        m <- mods[[idx]]
        expect_match(capture.output(summary(m))[2], "Summary of 3 MCMC draws after burn-in of 2")
    }
})

test_that("Plot parameters", {
    # Plot piecewise constant, time-varying parameter
    for (idx in seq_along(mods)) {
        m <- mods[[idx]]
        expect_invisible(plot(m))
        expect_invisible(plot(m$beta$beta_radiation))
        expect_invisible(plot(m$beta[[1]]))
        expect_invisible(plot(m$beta[["beta_radiation"]]))
    }
})

test_that("Prediction and predictive plots", {
    newdata <- data.frame(radiation = c(0, 1))
    for (idx in seq_along(mods)) {
        m <- mods[[idx]]
        p <- predict(m, newdata)
        expect_s3_class(p, "shrinkDSM_pred")
        expect_visible(p)
        expect_invisible(plot(p))
    }
})
