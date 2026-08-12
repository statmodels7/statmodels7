# A structural term rewrites the predictor instead of contributing columns.
# What has to be right is the gradient: the level a filter produces was
# driven by scores read at predictors the coefficients also enter, so the
# derivative of the objective in a coefficient is not the block times the
# score, and a fit built on that converges to the wrong point.

sim_gas <- function(n, omega, a1, b1, sd = 1, seed = 21) {
  set.seed(seed)
  f <- numeric(n)
  y <- numeric(n)
  s <- 0
  f0 <- omega / (1 - b1)
  for (t in seq_len(n)) {
    f[t] <- omega + a1 * (if (t > 1) s else 0) +
      b1 * (if (t > 1) f[t - 1] else f0)
    y[t] <- f[t] + stats::rnorm(1, sd = sd)
    s <- (y[t] - f[t]) / sd^2
  }
  data.frame(t = seq_len(n), y = y, f = f)
}

test_that("a structural term is left out of the design and routed", {
  dd <- sim_gas(80, 0.3, 0.4, 0.7)
  dd$x <- stats::rnorm(80)
  spec <- statmod_spec(y ~ x + gas(p = 1, q = 1, time = t) - 1,
                       distributions7::gaussian1_distrib(), dd)
  su <- statmod_structural(spec)
  expect_length(su, 1L)
  expect_identical(su[[1L]]$param, "mu")
  expect_identical(su[[1L]]$kind, "filter")
  # it contributes no columns: the design is x alone, the intercept
  # being removed because a filter carries a level of its own
  des <- statmod_design(spec)
  expect_identical(des$mu$npar, 1L)
  expect_false(su[[1L]]$term %in% names(des$mu$blocks))
  # and the shape is read off the methods the term registers
  expect_identical(structural_kind(
    modelterms7::term_build(modelterms7::regime(2), dd)), "loglik")
  expect_identical(structural_kind(
    modelterms7::term_build(modelterms7::linpar(~x), dd)), "")
})

test_that("the gradient in the coefficients carries the filter's feedback", {
  dd <- sim_gas(200, 0.3, 0.4, 0.7, seed = 5)
  dd$x <- stats::rnorm(200)
  spec <- statmod_spec(y ~ x + gas(p = 1, q = 1, time = t) - 1,
                       distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  sst <- statmod_structural_state(des)
  sst$zeta[[1L]] <- c(omega = 0.2, a1 = 0.3, pacf1 = 0.6)
  cf <- list(mu = 0.8, sigma = 0.1)

  L <- function(v) {
    sst$key <- NULL
    statmod_loglik_at(spec, list(mu = v[1], sigma = v[2]), des)
  }
  an <- unlist(statmod_score_at(spec, cf, des), use.names = FALSE)
  sst$key <- NULL
  num <- numDeriv::grad(L, c(cf$mu, cf$sigma))
  expect_equal(an, num, tolerance = 1e-6)

  # the naive gradient -- the block times the score, which is what every
  # term without a recursion gives -- is wrong here, and by a lot
  ev <- statmod_eta(spec, des, cf)
  g <- distributions7::distrib_gradient(spec@distrib, spec@response,
                                        ev$theta, scale = "link")
  naive <- c(as.numeric(crossprod(des$mu$X, g$mu)),
             as.numeric(crossprod(des$sigma$X, g$sigma)))
  expect_gt(max(abs(naive - num)), 1)
})

test_that("the gradient in the term's own parameters is exact", {
  dd <- sim_gas(200, 0.3, 0.4, 0.7, seed = 5)
  dd$x <- stats::rnorm(200)
  spec <- statmod_spec(y ~ x + gas(p = 1, q = 1, time = t) - 1,
                       distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  sst <- statmod_structural_state(des)
  cf <- list(mu = 0.8, sigma = 0.1)
  nm <- c("omega", "a1", "pacf1")
  z0 <- c(0.2, 0.3, 0.6)

  L <- function(z) {
    sst$zeta[[1L]] <- stats::setNames(z, nm)
    sst$key <- NULL
    statmod_loglik_at(spec, cf, des)
  }
  sst$zeta[[1L]] <- stats::setNames(z0, nm)
  sst$key <- NULL
  an <- statmod_structural_score(spec, cf, des)[[1L]]
  expect_equal(unname(an), numDeriv::grad(L, z0), tolerance = 1e-6)
})

test_that("statmod recovers a score-driven model it simulated", {
  truth <- c(omega = 0.3, a1 = 0.4, pacf1 = 0.7)
  dd <- sim_gas(2000, truth[["omega"]], truth[["a1"]], truth[["pacf1"]])
  fit <- statmod(y ~ gas(p = 1, q = 1, time = t) - 1,
                 distributions7::gaussian1_distrib(), dd)
  expect_true(fit@converged)
  est <- fit@structural[[1L]]$parameter
  expect_named(est, c("omega", "a1", "pacf1"))
  expect_equal(unname(est), unname(truth), tolerance = 0.15)
  expect_equal(exp(fit@coefficients$sigma), 1, tolerance = 0.1,
               ignore_attr = TRUE)
  # the unconstrained values are carried too, and map back
  z <- fit@structural[[1L]]$unconstrained
  expect_equal(unname(unlist(structural_psi(fit@spec@terms$mu[[1L]], z))),
               unname(est), tolerance = 1e-12)
  expect_identical(fit@structural[[1L]]$equation, "mu")

  # the objective fell at every sweep
  h <- fit@history$blocks
  expect_true(all(h$change >= -1e-8))
  expect_true(any(h$block == "structural"))
})

test_that("the fitted parameters are carried into a prediction", {
  dd <- sim_gas(300, 0.3, 0.4, 0.7, seed = 9)
  fit <- statmod(y ~ gas(p = 1, q = 1, time = t) - 1,
                 distributions7::gaussian1_distrib(), dd)
  # the specification the fit returns carries the parameters it arrived at,
  # so a design built from it starts there. Without them the filter would be
  # rebuilt at its starting values, which is the model without the term.
  expect_named(fit@spec@structural, names(fit@structural))
  z <- statmod_structural_state(statmod_design(fit@spec))$zeta[[1L]]
  expect_equal(z, fit@structural[[1L]]$unconstrained)
  expect_false(all(z == 0))

  # and predicting on the fitting rows returns the fitted values
  expect_equal(predict(fit, "mu", newdata = dd), fit@fitted$mu,
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("a model with no structural term carries none of this", {
  dd <- data.frame(y = stats::rnorm(60), x = stats::runif(60))
  spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
  expect_length(statmod_structural(spec), 0L)
  expect_null(statmod_structural_state(statmod_design(spec)))
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
  expect_length(fit@structural, 0L)
})

test_that("two structural terms are refused, and named", {
  dd <- sim_gas(80, 0.3, 0.4, 0.7)
  err <- tryCatch(
    statmod_spec(y ~ gas(1, 1, time = t) + gas(2, 1, time = t) - 1,
                 distributions7::gaussian1_distrib(), dd),
    error = conditionMessage)
  expect_match(err, "at most one structural term", fixed = TRUE)
  expect_match(err, "gas(1, 1, time = t)", fixed = TRUE)
  expect_match(err, "gas(2, 1, time = t)", fixed = TRUE)

  # in different equations too: each is driven by a score that depends on
  # the other's state, which is one recursion written as two
  err2 <- tryCatch(
    statmod_spec(y ~ gas(1, 1, time = t) - 1 | sigma ~ regime(2),
                 distributions7::gaussian1_distrib(), dd),
    error = conditionMessage)
  expect_match(err2, "at most one structural term", fixed = TRUE)
})

# --- the likelihood shape ---------------------------------------------------

sim_regime <- function(n, mu, P, sd = 1, seed = 3) {
  set.seed(seed)
  k <- length(mu)
  s <- integer(n)
  s[1L] <- 1L
  for (t in 2:n) s[t] <- sample.int(k, 1L, prob = P[s[t - 1L], ])
  data.frame(t = seq_len(n), y = stats::rnorm(n, mu[s], sd), state = s)
}
reg_P <- matrix(c(0.93, 0.07, 0.06, 0.94), 2, 2, byrow = TRUE)

test_that("a regime term's likelihood is the term's, not the density's", {
  dd <- sim_regime(400, c(0, 3), reg_P)
  spec <- statmod_spec(y ~ regime(2, time = t) - 1,
                       distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  sst <- statmod_structural_state(des)
  nmz <- names(sst$zeta[[1L]])
  sst$zeta[[1L]] <- stats::setNames(c(0, log(3), 2, -2), nmz)
  sst$key <- NULL
  cf <- list(mu = numeric(0), sigma = 0.1)

  ll <- statmod_loglik_at(spec, cf, des)
  # the density at the posterior-weighted predictor is a DIFFERENT number,
  # and a larger one: a mixture is not the density at the mixed point
  ev <- statmod_eta(spec, des, cf)
  naive <- sum(distributions7::distrib_pdf(spec@distrib, spec@response,
                                           ev$theta, log = TRUE))
  expect_false(isTRUE(all.equal(ll, naive)))
  expect_length(ev$regimes, 1L)
  expect_identical(dim(ev$regimes[[1L]]$gamma), c(400L, 2L))
  expect_equal(rowSums(ev$regimes[[1L]]$gamma), rep(1, 400),
               tolerance = 1e-12)
})

test_that("the regime gradient is Fisher's identity, and is exact", {
  dd <- sim_regime(400, c(0, 3), reg_P)
  dd$x <- stats::rnorm(400)
  spec <- statmod_spec(y ~ x + regime(2, time = t) - 1 | sigma ~ 1,
                       distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  sst <- statmod_structural_state(des)
  nmz <- names(sst$zeta[[1L]])
  z0 <- stats::setNames(c(0, log(3), 2, -2), nmz)
  setz <- function(z) {
    sst$zeta[[1L]] <- stats::setNames(z, nmz)
    sst$key <- NULL
  }
  cf <- list(mu = 0.5, sigma = 0.05)

  L <- function(v) {
    setz(z0)
    statmod_loglik_at(spec, list(mu = v[1L], sigma = v[2L]), des)
  }
  setz(z0)
  an <- unlist(statmod_score_at(spec, cf, des), use.names = FALSE)
  setz(z0)
  expect_equal(an, numDeriv::grad(L, c(cf$mu, cf$sigma)), tolerance = 1e-5)

  Lz <- function(z) {
    setz(z)
    statmod_loglik_at(spec, cf, des)
  }
  setz(z0)
  anz <- statmod_structural_score(spec, cf, des)[[1L]]
  setz(z0)
  expect_equal(unname(anz), numDeriv::grad(Lz, as.numeric(z0)),
               tolerance = 1e-6)
})

test_that("statmod recovers a regime model it simulated", {
  dd <- sim_regime(600, c(0, 3), reg_P)
  fit <- statmod(y ~ regime(2, time = t) - 1,
                 distributions7::gaussian1_distrib(), dd)
  expect_true(fit@converged)
  est <- fit@structural[[1L]]$parameter
  lev <- cumsum(c(est[["level1"]], est[["gap2"]]))
  expect_equal(unname(lev), c(0, 3), tolerance = 0.3)
  expect_equal(exp(fit@coefficients$sigma), 1, tolerance = 0.15,
               ignore_attr = TRUE)
  # the levels are ordered by construction, so the answer can be reported
  expect_lt(lev[1L], lev[2L])
})

test_that("a regime's first level is held by a linear intercept too", {
  dd <- sim_regime(600, c(0, 3), reg_P)
  dd$x <- stats::rnorm(600)
  # level1 shifts every state, so it is an intercept by another name; the
  # GAPS are not, being differences between regimes, and stay free
  fit <- statmod(y ~ x + regime(2, time = t),
                 distributions7::gaussian1_distrib(), dd)
  expect_true(fit@converged)
  expect_identical(fit@structural[[1L]]$held, "level1")
  est <- fit@structural[[1L]]$parameter
  expect_identical(est[["level1"]], 0)
  expect_equal(unname(est[["gap2"]]), 3, tolerance = 0.3)
  expect_equal(unname(fit@coefficients$mu[1L]), 0, tolerance = 0.3)

  # with no constant in the equation the level is free again
  d2 <- statmod_design(statmod_spec(y ~ x + regime(2, time = t) - 1,
                                    distributions7::gaussian1_distrib(), dd))
  expect_length(statmod_structural_state(d2)$held[[1L]], 0L)
})

# --- the observed information ----------------------------------------------

test_that("a linear intercept carries the level, and the term's is held", {
  dd <- sim_gas(80, 0.3, 0.4, 0.7)
  dd$x <- stats::rnorm(80)
  dd$g <- factor(rep(c("a", "b"), 40))
  held_of <- function(fm) {
    d <- statmod_design(statmod_spec(fm, distributions7::gaussian1_distrib(),
                                     dd))
    statmod_structural_state(d)$held[[1L]]
  }
  # exactly confounded, so one of the two goes: the linear intercept wins
  expect_identical(held_of(y ~ x + gas(1, 1, time = t)), "omega")
  # the question is asked of the SPAN and not of a column's name, so a
  # factor coded without an intercept holds it too
  expect_identical(held_of(y ~ g + gas(1, 1, time = t) - 1), "omega")
  # and where the equation spans no constant the level is free
  expect_length(held_of(y ~ x + gas(1, 1, time = t) - 1), 0L)
  # an intercept in ANOTHER equation is not this equation's
  expect_length(held_of(y ~ x + gas(1, 1, time = t) - 1 | sigma ~ 1), 0L)
})

test_that("the two parametrizations are the same model", {
  dd <- sim_gas(400, 0.3, 0.4, 0.7, seed = 11)
  dd$x <- stats::rnorm(400)
  a <- statmod(y ~ x + gas(p = 1, q = 1, time = t),
               distributions7::gaussian1_distrib(), dd)
  b <- statmod(y ~ x + gas(p = 1, q = 1, time = t) - 1,
               distributions7::gaussian1_distrib(), dd)
  expect_true(a@converged)
  expect_true(b@converged)
  expect_identical(a@structural[[1L]]$held, "omega")
  expect_identical(a@structural[[1L]]$parameter[["omega"]], 0)

  # holding the level and letting the intercept carry it reaches the same
  # maximum as leaving the level free with no intercept
  expect_equal(a@loglik, b@loglik, tolerance = 1e-6)
  # and the intercept IS the level the other fit found, through the
  # stationary mean of the recursion, omega / (1 - b)
  est <- b@structural[[1L]]$parameter
  expect_equal(unname(a@coefficients$mu[1L]),
               unname(est[["omega"]] / (1 - est[["pacf1"]])),
               tolerance = 1e-3)

  # and the variance exists, where inverting a singular matrix would not
  V <- vcov(a)
  expect_true(all(is.finite(V)))
  expect_true(all(diag(V) > 0))
})

test_that("the observed information carries the recursion's curvature", {
  dd <- sim_gas(150, 0.3, 0.4, 0.7, seed = 5)
  dd$x <- stats::rnorm(150)
  dd$z <- stats::runif(150)
  spec <- statmod_spec(y ~ x + gas(p = 1, q = 1, time = t) - 1 | sigma ~ z,
                       distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  sst <- statmod_structural_state(des)
  nmz <- names(sst$zeta[[1L]])
  params <- spec@distrib@params
  npar <- vapply(des, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  nb <- sum(npar)
  np <- length(nmz)

  split_u <- function(u) {
    list(coef = stats::setNames(lapply(seq_along(params), function(a)
           u[offs[a] + seq_len(npar[a])]), params),
         zeta = u[nb + seq_len(np)])
  }
  setz <- function(z) {
    sst$zeta[[1L]] <- stats::setNames(z, nmz)
    sst$key <- NULL
  }
  grad_full <- function(u) {
    s <- split_u(u)
    setz(s$zeta)
    c(unlist(statmod_score_at(spec, s$coef, des), use.names = FALSE),
      statmod_structural_score(spec, s$coef, des)[[1L]])
  }
  # AWAY from the optimum: near it a wrong hessian and a right one both look
  # plausible against a gradient that is nearly zero
  u0 <- c(0.8, -0.3, 0.25, 0.2, 0.3, 0.6)

  setz(u0[nb + seq_len(np)])
  I <- statmod_full_information(spec, split_u(u0)$coef, des)
  setz(u0[nb + seq_len(np)])
  num <- numDeriv::jacobian(grad_full, u0)
  expect_identical(dim(I), c(nb + np, nb + np))
  expect_lt(max(abs(I + num)) / max(abs(num)), 1e-7)
  expect_identical(I, t(I))

  # and the naive curvature, which is what the scoring step inverts, is a
  # long way from it: that gap is the whole reason this exists
  setz(u0[nb + seq_len(np)])
  naive <- statmod_information_at(spec, split_u(u0)$coef, des,
                                  expected = FALSE)
  blk <- I[seq_len(nb), seq_len(nb)]
  expect_gt(max(abs(naive - blk)) / max(abs(blk)), 0.1)
})

test_that("vcov inverts the joint matrix, not the coefficient block", {
  dd <- sim_gas(400, 0.3, 0.4, 0.7, seed = 11)
  dd$x <- stats::rnorm(400)
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t) - 1,
                 distributions7::gaussian1_distrib(), dd)
  expect_true(fit@converged)
  V <- vcov(fit)
  des <- statmod_design(fit@spec)
  nb <- sum(vapply(des, function(d) d$npar, integer(1)))
  expect_identical(dim(V), c(nb, nb))
  expect_true(all(is.finite(V)))
  expect_true(all(diag(V) > 0))

  # the joint route accounts for estimating the filter's parameters, so it
  # cannot be the inverse of the coefficient block alone
  naive <- solve(statmod_information_at(fit@spec, fit@coefficients, des,
                                        FALSE))
  expect_gt(max(abs(diag(V) - diag(naive))) / max(abs(diag(V))), 1e-3)

  # An estimated nuisance costs information, and what that says is a
  # statement about the JOINT matrix and not about the naive one: the
  # coefficient block of the inverse is (A - B C^-1 B')^-1, which dominates
  # A^-1 in the positive semidefinite order. The naive curvature is not A,
  # so it enters no such inequality -- a claim that it did was wrong here
  # before the test refused it.
  I <- statmod_full_information(fit@spec, fit@coefficients, des)
  nbb <- seq_len(nb)
  expect_true(all(diag(V) >= diag(solve(I[nbb, nbb, drop = FALSE])) - 1e-10))
})

test_that("a regime model reports the observed information, not the EM one", {
  dd <- sim_regime(300, c(0, 1.6), reg_P, seed = 12)
  dd$x <- stats::rnorm(300)
  spec <- statmod_spec(y ~ x + regime(2, time = t) - 1 | sigma ~ 1,
                       distributions7::gaussian1_distrib(), dd)
  fit <- statmod(y ~ x + regime(2, time = t) - 1 | sigma ~ 1,
                 distributions7::gaussian1_distrib(), dd)
  expect_true(fit@converged)

  spec <- fit@spec
  design <- statmod_design(spec)
  su <- attr(design, "structural")[[1L]]
  zeta <- statmod_structural_state(design)$zeta
  tn <- su$term
  zn <- modelterms7::term_params(spec@terms[[su$param]][[tn]])
  free <- setdiff(seq_along(zn),
                  match(statmod_structural_state(design)$held[[tn]], zn))
  nb <- sum(vapply(design, function(d) d$npar, integer(1)))
  b0 <- unlist(fit@coefficients[spec@distrib@params], use.names = FALSE)

  put <- function(u) {
    cf <- fit@coefficients
    i <- 0L
    for (p in spec@distrib@params) {
      k <- design[[p]]$npar
      if (k > 0L) { cf[[p]] <- u[i + seq_len(k)]; i <- i + k }
    }
    zz <- zeta[[tn]]
    zz[free] <- u[nb + seq_along(free)]
    dg <- design
    st <- statmod_structural_state(dg)
    st$zeta[[tn]] <- zz
    attr(dg, "structure") <- st
    list(coef = cf, design = dg)
  }
  gexact <- function(u) {
    a <- put(u)
    c(unlist(statmod_score_at(spec, a$coef, a$design), use.names = FALSE),
      statmod_structural_score(spec, a$coef, a$design)[[tn]][free])
  }

  # AWAY from the optimum, where a wrong Hessian and a right one differ;
  # at the optimum both look plausible against a vanishing gradient
  u0 <- c(b0, zeta[[tn]][free]) + 0.3
  a <- put(u0)
  H <- numDeriv::jacobian(gexact, u0)
  H <- (H + t(H)) / 2
  I <- statmod_full_information(spec, a$coef, a$design)
  expect_equal(dim(I), c(length(u0), length(u0)))
  expect_equal(-I, H, tolerance = 1e-4)
  expect_identical(I, t(I))

  # and the complete-data matrix, which is what a scoring step inverts,
  # is a different matrix
  Ic <- statmod_information_at(spec, a$coef, a$design, expected = FALSE)
  kb <- seq_len(nb)
  expect_gt(max(abs(Ic - I[kb, kb])) / max(abs(I[kb, kb])), 0.05)

  # at the fitted point it EXCEEDS the observed one, by the
  # missing-information principle, so its standard errors are too small
  a0 <- put(c(b0, zeta[[tn]][free]))
  Io <- statmod_full_information(spec, a0$coef, a0$design)
  Ic0 <- statmod_information_at(spec, a0$coef, a0$design, expected = FALSE)
  D <- Ic0 - Io[kb, kb]
  expect_gt(min(eigen(D, symmetric = TRUE, only.values = TRUE)$values), 0)
  se_o <- sqrt(diag(solve(Io))[kb])
  se_c <- sqrt(diag(solve(Ic0)))
  expect_true(all(se_c <= se_o + 1e-12))
  # vcov() reports the observed one, over the coefficient block of the
  # JOINT inverse
  expect_equal(unname(sqrt(diag(vcov(fit)))), unname(se_o), tolerance = 1e-8)
})

test_that("a structural term is counted and reported", {
  dd <- sim_regime(400, c(0, 3), reg_P, seed = 21)
  dd$x <- stats::rnorm(400)
  # the intercept-wins spelling: no '- 1', so the equation keeps its
  # constant and the term's level is held
  fit <- statmod(y ~ x + regime(2, time = t), gaussian1_distrib(), dd)
  expect_true(fit@converged)
  expect_identical(fit@structural[[1L]]$held, "level1")

  # the term's own free parameters are COUNTED. Its design block is empty,
  # so counting columns gave it zero and every criterion built on the total
  # was that much too generous.
  zn <- modelterms7::term_params(fit@spec@terms$mu[[2L]])
  free <- length(zn) - 1L
  row <- fit@edf[fit@edf$term != "linpar" & fit@edf$parameter == "mu", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$edf, as.numeric(free))
  # mu's linpar carries the intercept and x, sigma's the intercept
  expect_equal(sum(fit@edf$edf), 2 + free + 1)

  # and they are REPORTED, with a standard error from the joint information
  s <- summary(fit)
  st <- s@structural
  expect_false(is.null(st))
  expect_setequal(st$name, zn)
  expect_true(st$held[st$name == "level1"])
  expect_true(is.na(st$se[st$name == "level1"]))
  expect_true(all(is.finite(st$se[!st$held])))
  # an interval built on the unconstrained scale and mapped back keeps a
  # gap positive, which one built on the parameter scale need not
  expect_true(all(st$lower[st$name == "gap2"] > 0))
  expect_true(all(st$estimate[!st$held] > st$lower[!st$held] &
                    st$estimate[!st$held] < st$upper[!st$held]))
  out <- utils::capture.output(print(s))
  expect_true(any(grepl("structural", out)))
  expect_true(any(grepl("gap2", out)))
  expect_true(any(grepl("conditional log-likelihood", out)))
})

test_that("logLik reports which likelihood it is", {
  dd <- sim_regime(300, c(0, 3), reg_P, seed = 22)
  fit <- statmod(y ~ regime(2, time = t) - 1, gaussian1_distrib(), dd)
  lc <- logLik(fit)
  expect_equal(as.numeric(lc), fit@loglik)
  expect_equal(attr(lc, "df"), sum(fit@edf$edf))
  # a model whose hyperparameters were never chosen by a marginal criterion
  # has no marginal likelihood, and says so rather than returning one
  expect_error(logLik(fit, type = "marginal"), "marginal log-likelihood")
})

test_that("the two likelihoods of a penalized fit are both available", {
  set.seed(23)
  m <- 40L; per <- 12L; n <- m * per
  dd <- data.frame(g = factor(rep(seq_len(m), each = per)),
                   x = stats::rnorm(n))
  b <- stats::rnorm(m, sd = 0.8)
  dd$y <- 1 + 0.7 * dd$x + b[as.integer(dd$g)] + stats::rnorm(n)
  fit <- statmod(y ~ x + random(~ 1 | g), gaussian1_distrib(), dd,
                 outer_method = ml())

  lc <- logLik(fit)
  lm_ <- logLik(fit, type = "marginal")
  # the conditional one is read at the fitted coefficients and is the
  # LARGER; the marginal integrates the random effects away
  expect_gt(as.numeric(lc), as.numeric(lm_))
  # and the counts belong to their own likelihood: an effective count for
  # the conditional, the estimated parameters for the marginal
  expect_gt(attr(lc, "df"), 10)
  expect_lt(attr(lm_, "df"), 10)
  expect_equal(as.numeric(lm_), as.numeric(fit@criterion))
})
