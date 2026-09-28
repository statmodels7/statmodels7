# A coordinate at the edge of its chart can look stationary on the free scale
# without being at a maximum. On MASS::geyser, `waiting ~ regime(k = 2)`
# stopped with a transition probability at 2.8e-12 and a log-likelihood of
# -1134.01, while moving that probability inside raises it; the maximum is
# -1099.63.

geyser_df <- function() data.frame(waiting = geyser_waiting)
G <- distributions7::gaussian1_distrib()

# the regime term's start before modelterms7 0.81.0: zeros on every chart
zero_start <- function(term, ...) {
  nm <- modelterms7::term_params(term)
  stats::setNames(numeric(length(nm)), nm)
}

with_zero_start <- function(expr) {
  gen <- modelterms7::term_start
  RT <- modelterms7::RegimeTerm
  keep <- S7::method(gen, RT)
  S7::method(gen, RT) <- zero_start
  on.exit(S7::method(gen, RT) <- keep, add = TRUE)
  suppressMessages(expr)
}

# the design of a fit, moved to the point the old start stopped at
false_point <- function(fit) {
  design <- statmod_design(fit@spec)
  sst <- statmod_structural_state(design)
  tn <- names(fit@structural)[1]
  z <- sst$zeta[[tn]]
  z[c("level1", "gap2", "alr1.1", "alr2.1")] <-
    c(0, 3.015777, -2.105419, 26.600657)
  sst$zeta[[tn]] <- z
  cf <- fit@coefficients
  cf$mu[[1]] <- 62.7206
  cf$sigma[[1]] <- 2.2411
  list(design = design, coef = cf)
}

test_that("the false point is flagged and a true boundary maximum is not", {
  fit <- statmod(waiting ~ regime(k = 2), G, geyser_df())
  fp <- false_point(fit)
  v <- edge_violations(fit@spec, fp$design, fp$coef)
  expect_equal(nrow(v), 1L)
  expect_identical(v$name, "alr2.1")
  expect_equal(v$target, 8)
  # moving P22 from 2.8e-12 to 3.4e-4 gains about 0.06
  expect_gt(v$gain, 0.03)
  expect_lt(v$gain, 0.1)
  # at the global maximum P11 sits at its edge and points outward
  f0 <- statmod(waiting ~ 0 + regime(k = 2), G, geyser_df())
  expect_gt(abs(f0@structural[[1]]$unconstrained[["alr1.1"]]), 8)
  expect_equal(nrow(edge_violations(f0@spec, statmod_design(f0@spec),
                                    f0@coefficients)), 0L)
})

test_that("the fit restarts from the edge and reaches the maximum", {
  # negative control: with the old start and the check switched off, the fit
  # stops at the false point, so the case really reaches the code
  local_mocked_bindings(edge_violations = function(...) {
    data.frame(kind = character(0), param = character(0),
               term = character(0), name = character(0), eta = numeric(0),
               target = numeric(0), gain = numeric(0))
  })
  old <- with_zero_start(statmod(waiting ~ regime(k = 2), G, geyser_df()))
  expect_lt(as.numeric(logLik(old)), -1130)
})

test_that("with the check, the old start still ends at the maximum", {
  fit <- with_zero_start(statmod(waiting ~ regime(k = 2), G, geyser_df()))
  expect_gt(as.numeric(logLik(fit)), -1099.7)
  expect_equal(nrow(fit@methods$kkt), 0L)
  expect_identical(statmod_certificate(fit)$state, "converged")
  expect_no_error(capture.output(summary(fit)))
})

test_that("the quantile start reaches the maximum by itself", {
  fit <- statmod(waiting ~ regime(k = 2), G, geyser_df())
  expect_gt(as.numeric(logLik(fit)), -1099.7)
  expect_no_error(capture.output(summary(fit)))
})

test_that("a violation that survives is reported by the certificate", {
  fit <- statmod(waiting ~ 0 + regime(k = 2), G, geyser_df())
  expect_identical(statmod_certificate(fit)$state, "converged")
  kk <- data.frame(kind = "structural", param = "mu",
                   term = names(fit@structural)[1], name = "alr2.1",
                   eta = 26.6, target = 8, gain = 0.06)
  m <- fit@methods
  m$kkt <- kk
  bad <- S7::set_props(fit, methods = m)
  cer <- statmod_certificate(bad)
  expect_identical(cer$state, "not converged")
  expect_match(cer$reason[1], "alr2.1", fixed = TRUE)
})
