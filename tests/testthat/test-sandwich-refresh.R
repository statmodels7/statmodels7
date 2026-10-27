g1r <- distributions7::gaussian1_distrib()
quiet_r <- function(e) suppressWarnings(suppressMessages(e))

# The log-likelihood of each observation as a function of the whole vector of
# coefficients, with the block of the moving term rebuilt at every point, and
# its numerical Jacobian. It shares no code with the analytic scores.
obs_loglik <- function(fit) {
  spec <- fit@spec
  design <- statmod_design(spec)
  params <- spec@distrib@params
  lens <- vapply(fit@coefficients[params], length, integer(1))
  grp <- factor(rep(params, lens), levels = params)
  function(beta) {
    cf <- split(beta, grp)
    cf <- Map(function(v, p) stats::setNames(v, names(fit@coefficients[[p]])),
              cf, params)
    d <- statmod_design_at(spec, cf, design)
    ev <- statmod_eta(spec, d, cf)
    spec@weights * distributions7::distrib_pdf(spec@distrib, spec@response,
                                              ev$theta, log = TRUE)
  }
}

hand_hc0 <- function(fit, keep) {
  f <- obs_loglik(fit)
  b0 <- unlist(fit@coefficients, use.names = FALSE)
  G <- numDeriv::jacobian(f, b0)[, keep, drop = FALSE]
  V <- vcov(fit, vcov_type = bayesian_vcov(
    information = distributions7::observed_info()), readable = FALSE)
  Vb <- V[keep, keep, drop = FALSE]
  Vb %*% crossprod(G) %*% Vb
}

fsw_r <- function(sw) {
  frequentist_vcov(sandwich = sw, information = distributions7::observed_info())
}

test_that("HC0 of a segmented fit is the sandwich of the numerical scores", {
  skip_if_not_installed("numDeriv")
  set.seed(8)
  n <- 150
  d <- data.frame(x = runif(n))
  d$y <- ifelse(d$x < 0.5, d$x, 0.5 + 2 * (d$x - 0.5)) + rnorm(n, sd = 0.1)
  f <- quiet_r(statmod(y ~ seg(x, psi = 0.4), g1r, d))
  expect_false(is.null(attr(statmod_design(f@spec), "refresh")))
  V <- vcov(f, vcov_type = fsw_r(hc0_sandwich()), readable = FALSE)
  keep <- !is.na(diag(V))
  expect_gt(sum(keep), 3L)
  ref <- hand_hc0(f, keep)
  expect_equal(unname(V[keep, keep]), unname(ref), tolerance = 1e-5)
  # the break-point's own coordinate has a finite positive variance
  expect_true(all(diag(V)[keep] > 0))
  # the types and the readable quantities take the same route
  for (sw in list(hc1_sandwich(), hc2_sandwich(), hc3_sandwich(),
                  hc0_sandwich(weights = "sampling"))) {
    Vs <- vcov(f, vcov_type = fsw_r(sw), readable = FALSE)
    expect_equal(is.na(diag(Vs)), is.na(diag(V)))
  }
  expect_true(all(is.finite(confint(f, vcov_type = fsw_r(hc1_sandwich()))$lower)))
})

test_that("HC0 of a nonlinear term is the sandwich of the numerical scores", {
  skip_if_not_installed("numDeriv")
  set.seed(3)
  n <- 120
  d <- data.frame(x = runif(n))
  d$y <- 2 * exp(1.2 * d$x) + rnorm(n, sd = 0.2)
  f <- quiet_r(statmod(y ~ 0 + nl(~ a * exp(r * x), start = list(a = 1.5, r = 1)),
                       g1r, d))
  V <- vcov(f, vcov_type = fsw_r(hc0_sandwich()), readable = FALSE)
  keep <- !is.na(diag(V))
  expect_gt(sum(keep), 2L)
  ref <- hand_hc0(f, keep)
  expect_equal(unname(V[keep, keep]), unname(ref), tolerance = 1e-5)
})

test_that("a sharp step whose position is not held is rejected", {
  set.seed(8)
  n <- 150
  d <- data.frame(x = runif(n))
  d$y <- ifelse(d$x < 0.5, 0, 1) + rnorm(n, sd = 0.2)
  fj <- quiet_r(statmod(y ~ jump(x, psi = 0.4), g1r, d))
  rf <- attr(statmod_design(fj@spec), "refresh")
  frozen <- any(vapply(rf, function(r) isTRUE(r$frozen), logical(1)))
  if (frozen) {
    expect_error(vcov(fj, vcov_type = fsw_r(hc0_sandwich())), "sharp")
  } else {
    expect_no_error(vcov(fj, vcov_type = fsw_r(hc0_sandwich())))
  }
})

test_that("a segmented fit takes a HAC variance with a time variable", {
  set.seed(8)
  n <- 150
  d <- data.frame(x = runif(n), t = seq_len(n))
  d$y <- ifelse(d$x < 0.5, d$x, 0.5 + 2 * (d$x - 0.5)) + rnorm(n, sd = 0.1)
  f <- quiet_r(statmod(y ~ seg(x, psi = 0.4), g1r, d))
  V <- vcov(f, vcov_type = frequentist_vcov(
    sandwich = hac_sandwich(time = ~ t, bandwidth = 0)), readable = FALSE)
  V0 <- vcov(f, vcov_type = frequentist_vcov(sandwich = hc0_sandwich()),
             readable = FALSE)
  # a bandwidth of zero leaves the diagonal terms, which are those of HC0
  expect_equal(V, V0, tolerance = 1e-8)
})
