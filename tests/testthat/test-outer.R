# Estimating the hyperparameters by a marginal criterion.

set.seed(21)
n <- 300
ds <- data.frame(x = runif(n, -2, 2))
ds$y <- sin(1.4 * ds$x) + stats::rnorm(n, sd = 0.3)

test_that("the criterion is the Laplace formula, assembled independently", {
  skip_if_not_installed("numDeriv")
  # the observed information is asked for so that numDeriv's Hessian of the
  # penalized objective is the same matrix the criterion uses
  fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), ds,
                 outer_method = reml(hessian = "observed"))
  spec <- fit@spec
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, fit@hyper, design, expected = FALSE)
  b <- unlist(fit@coefficients[spec@distrib@params], use.names = FALSE)

  ll <- as.numeric(stats::logLik(fit))
  rho <- statmod_penalty_at(spec, fit@coefficients, fit@hyper, design, "value")
  M <- numDeriv::hessian(obj$fn, b)
  by_hand <- ll - rho + length(b) / 2 * log(2 * pi) -
    determinant(M, logarithm = TRUE)$modulus[[1L]] / 2

  m <- statmod_marginal(spec, design, fit@coefficients, fit@hyper,
                        reml(hessian = "observed"))
  expect_equal(m$value, as.numeric(by_hand), tolerance = 1e-5)
  expect_equal(m$value, fit@criterion, tolerance = 1e-8)
  expect_identical(m$q, length(b))
})

test_that("REML is Wood's criterion, reached by the other route", {
  # ell - rho + (p/2) log 2pi - (1/2) log|H+S| must equal
  # ell - (lambda/2) b'Pb + (1/2) log|lambda P|_+ - (1/2) log|H+S|
  #      + (M_p/2) log 2pi,   M_p = p - rank
  # the two differing only by how the penalty's normalizing constant is
  # written. Keeping that constant is what makes the first form correct with
  # nothing added by hand, so the identity is worth asserting rather than
  # assuming. The sign of the last term is the whole content of the check:
  # written the other way round the two forms are apart by exactly
  # (p - r) log 2pi, which here is the intercept, the smooth's linear
  # component and the scale's intercept.
  fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), ds,
                 outer_method = reml())
  spec <- fit@spec
  design <- statmod_design(spec)
  nm <- names(spec@terms$mu)[vapply(spec@terms$mu, function(t)
    !is.null(modelterms7::term_penalty(t)), logical(1))]
  pen <- modelterms7::term_penalty(spec@terms$mu[[nm]])
  th <- as.list(fit@hyper$mu[[nm]])
  bt <- fit@coefficients$mu[design$mu$blocks[[nm]]]

  P <- penalties7::penalty_matrix(pen, th)
  quad <- as.numeric(crossprod(bt, P %*% bt)) / 2
  lpd <- penalties7::penalty_logpdet(pen, th)$value
  r <- penalties7::penalty_rank(pen)
  p_tot <- length(unlist(fit@coefficients[spec@distrib@params]))

  m <- statmod_marginal(spec, design, fit@coefficients, fit@hyper, reml())
  H <- statmod_information_at(spec, fit@coefficients, design, TRUE)
  S <- statmod_penalty_at(spec, fit@coefficients, fit@hyper, design, "hessian")
  logdet <- determinant(H + S, logarithm = TRUE)$modulus[[1L]]
  wood <- as.numeric(stats::logLik(fit)) - quad + lpd / 2 - logdet / 2 +
    (p_tot - r) / 2 * log(2 * pi)

  expect_equal(m$value, wood, tolerance = 1e-8)
  # and the count of unpenalized directions is what that term carries
  expect_equal(p_tot - r, 3)
})

test_that("the reported hyperparameter is where the criterion is best", {
  fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), ds,
                 outer_method = reml())
  nm <- names(fit@hyper$mu)[1L]
  lam <- fit@hyper$mu[[nm]][["lambda"]]
  at <- function(v) {
    h <- fit@hyper
    h$mu[[nm]][["lambda"]] <- v
    statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), ds,
            hyper = list(mu = stats::setNames(list(c(lambda = v)), nm)))
  }
  crit <- function(f) {
    statmod_marginal(f@spec, statmod_design(f@spec), f@coefficients,
                     f@hyper, reml())$value
  }
  best <- fit@criterion
  expect_gt(best, crit(at(lam * 3)))
  expect_gt(best, crit(at(lam / 3)))
  expect_gt(best, crit(at(lam * 1.4)))
  expect_gt(best, crit(at(lam / 1.4)))
})

test_that("a smoothing parameter estimated by REML lands on a sane edf", {
  fit <- statmod(y ~ s(x, k = 12), distributions7::gaussian1_distrib(), ds,
                 outer_method = reml())
  e <- fit@edf$edf[fit@edf$term != "linpar" & fit@edf$parameter == "mu"]
  # a sine over four periods of the covariate is neither a straight line nor
  # an interpolation of 300 points
  expect_gt(e, 3)
  expect_lt(e, 11)
  # and the fit is close to the truth it was drawn from
  expect_lt(sqrt(mean((predict(fit, "mu") - sin(1.4 * ds$x))^2)), 0.06)
  expect_true(is.finite(fit@criterion))
  expect_gt(nrow(fit@history$outer), 3)
})

test_that("ML puts less variance on a random effect than REML", {
  # the classic downward bias: profiling the fixed effects rather than
  # integrating them leaves the variance component too small
  set.seed(22)
  m <- 25
  g <- factor(rep(paste0("g", seq_len(m)), each = 10))
  u <- stats::rnorm(m, sd = 0.8)
  dv <- data.frame(g = g, x = runif(250))
  dv$y <- 1 + 2 * dv$x + u[as.integer(g)] + stats::rnorm(250, sd = 0.5)

  r <- statmod(y ~ x + random(~ 1 | g),
               distributions7::gaussian1_distrib(), dv, outer_method = reml())
  l <- statmod(y ~ x + random(~ 1 | g),
               distributions7::gaussian1_distrib(), dv, outer_method = ml())
  nm <- "random(~1 | g)"
  sd_r <- r@hyper$mu[[nm]][[1L]]
  sd_l <- l@hyper$mu[[nm]][[1L]]

  expect_lt(sd_l, sd_r)
  # both near the 0.8 they were drawn from
  expect_gt(sd_r, 0.5)
  expect_lt(sd_r, 1.2)
  # and the fixed effects are recovered
  expect_equal(r@coefficients$mu[2L], 2, tolerance = 0.25)
})

test_that("a kinked penalty keeps the hyperparameter it was given", {
  # a Laplace approximation asks for a second derivative, and a lasso has none
  # at zero, so its lambda is not something a marginal criterion can estimate
  set.seed(23)
  dl <- ds
  for (j in 1:3) dl[[paste0("n", j)]] <- stats::rnorm(n)
  fit <- statmod(y ~ s(x, k = 10) + lasso(~ n1 + n2 + n3),
                 distributions7::gaussian1_distrib(), dl,
                 hyper = list(mu = list(lasso = c(lambda = 40))),
                 outer_method = reml())
  expect_equal(unname(fit@hyper$mu[["lasso(~n1 + n2 + n3)"]][["lambda"]]), 40)
  # while the smooth's did move
  expect_false(isTRUE(all.equal(
    unname(fit@hyper$mu[["s(x, k = 10)"]][["lambda"]]), 1)))
})

test_that("an outer method with nothing to estimate says so", {
  expect_error(statmod(y ~ x, distributions7::gaussian1_distrib(), ds,
                       outer_method = reml()),
               "no hyperparameter to", fixed = TRUE)
  expect_error(statmod(y ~ s(x, k = 8), distributions7::gaussian1_distrib(),
                       ds, outer_method = "reml"),
               "reml(), ml() or NULL", fixed = TRUE)
})

test_that("ml refuses a penalty whose null space it cannot read", {
  # an anisotropic tensor smooth carries an additive penalty, which reports a
  # rank but exposes no null basis; guessing which directions are profiled
  # would give a criterion that is arithmetic without a meaning
  set.seed(24)
  dt <- data.frame(x1 = runif(200, -1, 1), x2 = runif(200, -1, 1))
  dt$y <- dt$x1^2 + dt$x2 + stats::rnorm(200, sd = 0.3)
  expect_error(statmod(y ~ te(x1, x2, k = 4),
                       distributions7::gaussian1_distrib(), dt,
                       outer_method = ml()),
               "cannot read the null space")
  # while reml integrates everything and needs no such basis
  fit <- statmod(y ~ te(x1, x2, k = 4),
                 distributions7::gaussian1_distrib(), dt,
                 outer_method = reml())
  expect_true(is.finite(fit@criterion))
  # one smoothing parameter per margin, both moved off the probe value
  th <- fit@hyper$mu[["te(x1, x2, k = 4)"]]
  expect_gte(length(th), 2L)
})

test_that("the summary marks an estimated hyperparameter as estimated", {
  fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), ds,
                 outer_method = reml())
  s <- summary(fit)
  kinds <- vapply(s@tables$mu, `[[`, character(1), "kind")
  b <- s@tables$mu[[which(kinds == "smooth")]]
  expect_identical(b$table$role[b$table$name == "lambda"], "estimated")
  expect_output(print(s), "(estimated)", fixed = TRUE)
  expect_true(any(grepl("REML", s@notes)))
  expect_output(print(fit), "REML")
})

test_that("verbose names the outer loop among its switches", {
  expect_true(verbosity(1)$outer)
  expect_false(verbosity(0)$outer)
  expect_true(verbosity(c(outer = TRUE))$outer)
  expect_false(verbosity(c(outer = TRUE))$blocks)
  expect_error(verbosity(c(wrong = TRUE)), "'outer'", fixed = TRUE)
})
