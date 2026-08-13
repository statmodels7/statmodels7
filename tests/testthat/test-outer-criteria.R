# The prediction-error criteria and their exact derivatives.

set.seed(61)
n <- 250
dc <- data.frame(x = runif(n, -2, 2), z = runif(n))
dc$y <- sin(1.4 * dc$x) + stats::rnorm(n, sd = 0.3)

pe_handles <- function(formula, data, method) {
  distrib <- distributions7::gaussian1_distrib()
  # The reference differentiates quantities read at the penalized MODE, so
  # whatever the inner fit leaves short of stationarity is noise the
  # numerical derivative amplifies; the default is the tightest usable
  # setting, the stall guard on the objective firing first below it (see
  # test-outer-hessian.R).
  inner <- iwls()
  # the probe value, not the estimated optimum: see test-outer-gradient.R
  fit0 <- statmod(formula, distrib, data, outer_criterion = NULL,
                  inner_optimizer = inner)
  spec <- fit0@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  at <- function(eta) {
    hy <- eta_to_hyper(eta, idx, fit0@hyper)
    hl <- lapply(hy, function(pp)
      lapply(pp, function(v) stats::setNames(as.numeric(v), names(v))))
    hl <- hl[lengths(hl) > 0L]
    list(hy = hy, fit = statmod(formula, distrib, data, hyper = hl,
                                inner_optimizer = inner))
  }
  list(
    idx = idx,
    eta0 = hyper_to_eta(fit0@hyper, idx),
    fn = function(eta) {
      a <- at(eta)
      statmod_pe(spec, design, a$fit@coefficients, a$hy, method)$value
    },
    gr = function(eta) {
      a <- at(eta)
      statmod_pe_derivs(spec, design, a$fit@coefficients, a$hy, method, idx,
                        1L)$grad
    },
    he = function(eta) {
      a <- at(eta)
      statmod_pe_derivs(spec, design, a$fit@coefficients, a$hy, method, idx,
                        2L)$hess
    })
}

test_that("the criterion is minus twice the log-likelihood plus k times edf", {
  fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dc,
                 outer_criterion = aic())
  spec <- fit@spec
  design <- statmod_design(spec)
  H <- statmod_information_at(spec, fit@coefficients, design, FALSE)
  S <- statmod_penalty_at(spec, fit@coefficients, fit@hyper, design, "hessian")
  tau <- sum(solve(H + S) * t(H))
  ll <- as.numeric(stats::logLik(fit))
  expect_equal(fit@criterion, -2 * ll + 2 * tau, tolerance = 1e-8)

  # an unpenalized model spends exactly its coefficients, which is what says
  # the trace is the right quantity
  plain <- statmod(y ~ x, distributions7::gaussian1_distrib(), dc)
  sp2 <- plain@spec
  d2 <- statmod_design(sp2)
  H2 <- statmod_information_at(sp2, plain@coefficients, d2, FALSE)
  S2 <- statmod_penalty_at(sp2, plain@coefficients, plain@hyper, d2, "hessian")
  expect_equal(outer_tau(H2 + S2, H2), 3, tolerance = 1e-10)
})

test_that("bic prices a degree of freedom at log n", {
  a <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dc,
               outer_criterion = aic())
  b <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dc,
               outer_criterion = bic())
  expect_equal(outer_k(aic(), n), 2)
  expect_equal(outer_k(bic(), n), log(n))
  # the dearer degree of freedom buys a smoother fit
  expect_gt(b@hyper$mu[["s(x, k = 10)"]][["lambda"]],
            a@hyper$mu[["s(x, k = 10)"]][["lambda"]])
  expect_lt(sum(b@edf$edf), sum(a@edf$edf))
})

test_that("the gradient of aic matches numDeriv", {
  skip_if_not_installed("numDeriv")
  h <- pe_handles(y ~ s(x, k = 10), dc, aic())
  for (shift in c(0.5, -1.2)) {
    eta <- h$eta0 + shift
    expect_equal(h$gr(eta), numDeriv::grad(h$fn, eta), tolerance = 1e-5)
  }
})

test_that("the Hessian of aic matches numDeriv", {
  skip_if_not_installed("numDeriv")
  h <- pe_handles(y ~ s(x, k = 10), dc, aic())
  eta <- h$eta0 + 0.4
  expect_equal(as.numeric(h$he(eta)),
               as.numeric(numDeriv::jacobian(h$gr, eta)), tolerance = 1e-4)
})

test_that("the derivatives match numDeriv with the scale modelled", {
  # a second distribution parameter makes the third derivative three-index and
  # gives the trace a block it would otherwise never see
  skip_if_not_installed("numDeriv")
  h <- pe_handles(y ~ s(x, k = 8) | sigma ~ z, dc, bic())
  eta <- h$eta0 + 0.2
  expect_equal(h$gr(eta), numDeriv::grad(h$fn, eta), tolerance = 1e-5)
  expect_equal(as.numeric(h$he(eta)),
               as.numeric(numDeriv::jacobian(h$gr, eta)), tolerance = 1e-4)
})

test_that("the derivatives match numDeriv with two smoothing parameters", {
  skip_if_not_installed("numDeriv")
  set.seed(62)
  n2 <- 250
  d2 <- data.frame(a = runif(n2, -2, 2), b = runif(n2, -2, 2))
  d2$y <- sin(1.4 * d2$a) + d2$b^2 + stats::rnorm(n2, sd = 0.3)
  h <- pe_handles(y ~ s(a, k = 8) + s(b, k = 8), d2, aic())
  eta <- h$eta0 + c(0.4, -0.5)
  expect_equal(h$gr(eta), numDeriv::grad(h$fn, eta), tolerance = 1e-5)
  got <- h$he(eta)
  expect_equal(got, t(got), tolerance = 1e-12)
  expect_equal(as.numeric(got),
               as.numeric(numDeriv::jacobian(h$gr, eta)), tolerance = 1e-4)
})

test_that("a variance component is covered by aic too", {
  skip_if_not_installed("numDeriv")
  set.seed(63)
  m <- 20
  dr <- data.frame(g = factor(rep(paste0("g", seq_len(m)), each = 10)),
                   x = runif(200))
  u <- stats::rnorm(m, sd = 0.7)
  dr$y <- 1 + 2 * dr$x + u[as.integer(dr$g)] + stats::rnorm(200, sd = 0.4)
  h <- pe_handles(y ~ x + random(~ 1 | g), dr, aic())
  # measured across four points, the exact gradient against numDeriv:
  #   eta  0.000   0.4785699 / 0.4785074   1.3e-4
  #   eta  0.200   0.4376933 / 0.4376832   2.3e-5
  #   eta  0.600   0.2571270 / 0.2571267   9.2e-7
  #   eta -0.500  -0.9987500 / -0.9989900  2.4e-4
  # the criterion is strongly curved here -- the gradient goes from 0.48 to
  # -1.00 over half a unit -- and that is what the reference's own step is
  # fighting, so the tolerance is set to what it can support. It moved when
  # the fits started from the intercept-only MLE instead of from zeros: the
  # reference refits the mode at every perturbation, so where the mode sits
  # is part of its accuracy.
  for (shift in c(0, 0.2, 0.6, -0.5)) {
    eta <- h$eta0 + shift
    expect_equal(h$gr(eta), numDeriv::grad(h$fn, eta), tolerance = 5e-4)
  }
  f <- statmod(y ~ x + random(~ 1 | g),
               distributions7::gaussian1_distrib(), dr, outer_criterion = aic())
  expect_true(is.finite(f@criterion))
})

test_that("the search minimizes a prediction-error criterion", {
  # the sign is the method's to declare, and getting it wrong would take the
  # search to the roughest fit rather than the best one
  expect_true(outer_minimize(aic()))
  expect_true(outer_minimize(bic()))
  expect_false(outer_minimize(reml()))

  fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dc,
                 outer_criterion = aic())
  spec <- fit@spec
  design <- statmod_design(spec)
  nm <- "s(x, k = 10)"
  at <- function(v) {
    f <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dc,
                 hyper = list(mu = stats::setNames(list(c(lambda = v)), nm)))
    statmod_pe(spec, design, f@coefficients, f@hyper, aic())$value
  }
  lam <- fit@hyper$mu[[nm]][["lambda"]]
  expect_lt(fit@criterion, at(lam * 3))
  expect_lt(fit@criterion, at(lam / 3))
  expect_lt(fit@criterion, at(lam * 1.3))
  expect_lt(fit@criterion, at(lam / 1.3))
})

test_that("aic and reml need not agree, and both are stationary", {
  # they estimate different things, so a difference is not a defect; what has
  # to hold is that each is at rest where it stopped
  a <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dc,
               outer_criterion = aic())
  r <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dc,
               outer_criterion = reml(hessian = "observed"))
  spec <- a@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  ga <- statmod_pe_derivs(spec, design, a@coefficients, a@hyper, aic(), idx,
                          1L)$grad
  gr <- statmod_marginal_grad(spec, design, r@coefficients, r@hyper,
                              reml(hessian = "observed"), idx)
  expect_lt(abs(ga), 1e-4)
  expect_lt(abs(gr), 1e-4)
  # and they do NOT land in the same place: measured on this data, 12.9
  # against 3.7. A first version of this asserted they would agree within
  # half, from a guess, and the measurement refused it -- which is the whole
  # reason a caller is given the choice.
  expect_gt(a@hyper$mu[[1L]][["lambda"]], r@hyper$mu[[1L]][["lambda"]])
})

test_that("aic is validated and prints", {
  expect_error(aic(k = -1), "non-negative")
  expect_error(aic(k = c(2, 3)), "single")
  expect_output(print(aic()), "AIC")
  expect_output(print(bic()), "BIC")
})

test_that("the gradient never asks for a second derivative", {
  # a penalty that supplies penalty_dhessian() and nothing beyond it must
  # still give an exact gradient: asking it for a derivative the gradient does
  # not use would reject it for a quantity nobody wanted
  spec <- statmod_spec(y ~ s(x, k = 8), distributions7::gaussian1_distrib(),
                       dc)
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar

  first <- outer_pieces(spec, design, statmod(y ~ s(x, k = 8),
                                              distributions7::gaussian1_distrib(),
                                              dc)@coefficients,
                        statmod_hyper_start(spec), idx, offs, sum(npar), 1L)
  expect_named(first, c("S", "c"))
  expect_null(first$S2)
})

test_that("the edf correction reproduces mgcv's smoothing-parameter term", {
  # The ordinary effective degrees of freedom read the smoothing parameter
  # as known; it was estimated from the same data. mgcv sums two terms into
  # edf2 and this is the first, J V_theta J' contracted with the
  # information. The second is mgcv's correction for the Gaussian scale,
  # which it profiles out and we model, so it has no counterpart here.
  skip_if_not_installed("mgcv")
  skip_on_cran()

  mgcv_vc1 <- function(g) {
    XWX <- crossprod(g$R)
    db <- g$db.drho
    if (!is.null(g$L)) db <- db %*% g$L[1:ncol(db), , drop = FALSE]
    ev <- eigen(g$outer.info$hess, symmetric = TRUE)
    d <- ev$values
    ind <- d <= 0
    d[ind] <- 0
    d[!ind] <- 1 / sqrt(d[!ind])
    rV <- (d * t(ev$vectors))[, 1:ncol(db), drop = FALSE]
    sum(rowSums(crossprod(rV %*% t(db)) * XWX) / g$sig2)
  }

  for (n in c(200L, 400L)) {
    set.seed(20)
    dd <- data.frame(x = stats::runif(n))
    dd$y <- sin(2 * pi * dd$x) + 0.4 * dd$x + stats::rnorm(n, sd = 0.3)
    # mgcv's s() cannot be namespace-qualified inside a formula, and
    # attaching mgcv would mask OUR s() for the statmod call below, so the
    # gam formula is given an environment where s is mgcv's and nothing else
    # is disturbed
    genv <- new.env(parent = globalenv())
    genv$s <- mgcv::s
    g <- mgcv::gam(stats::as.formula("y ~ s(x, bs = 'bs', k = 15)",
                                     env = genv),
                   data = dd, method = "REML")
    u <- statmod(y ~ s(x, k = 15), distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
    got <- summary(u, correct = TRUE)@df - summary(u)@df
    # ABSOLUTE, because the quantity is a fraction of a parameter and what
    # is left between the two packages is the difference between their
    # bases: the naive counts differ by the same amount as the corrected
    # ones, so the correction itself is what agrees
    expect_lt(abs(got - mgcv_vc1(g)), 2e-3)
    # and it is an ADDITION: the naive count is too generous
    expect_gt(got, 0)
  }
})

test_that("the edf correction is zero where nothing was estimated", {
  # A kinked penalty's hyperparameter is not chosen by a differentiable
  # criterion -- outer_hyper_index() skips it -- so there is nothing to
  # propagate. That is a refusal by construction, not an approximation.
  set.seed(23)
  p <- 12L
  X <- matrix(stats::rnorm(150 * p), 150, p)
  dd <- data.frame(y = X[, 1] * 2 - X[, 2] * 1.5 + stats::rnorm(150),
                   X = I(X))
  u <- statmod(y ~ lasso(X), distributions7::gaussian1_distrib(), dd,
               hyper = list(mu = list("lasso(X)" = c(lambda = 20))))
  expect_equal(summary(u, correct = TRUE)@df, summary(u)@df)
  expect_true(any(grepl("nothing for the correction",
                        summary(u, correct = TRUE)@notes)))

  # and a model with no penalty at all is unmoved
  d2 <- data.frame(x = stats::rnorm(100))
  d2$y <- 1 + 2 * d2$x + stats::rnorm(100)
  f2 <- statmod(y ~ x, distributions7::gaussian1_distrib(), d2)
  expect_equal(summary(f2, correct = TRUE)@df, summary(f2)@df)
})

test_that("the correction applies to a random effect, not only a smooth", {
  # a random effect IS a ridge penalty here, so the same code path covers it
  # with no term-specific branch anywhere
  set.seed(22)
  m <- 40L
  per <- 10L
  n <- m * per
  dd <- data.frame(g = factor(rep(seq_len(m), each = per)),
                   x = stats::rnorm(n))
  b <- stats::rnorm(m, sd = 0.8)
  dd$y <- 1 + 0.7 * dd$x + b[as.integer(dd$g)] + stats::rnorm(n)
  u <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(),
               dd, outer_criterion = reml())
  naive <- summary(u)@df
  corr <- summary(u, correct = TRUE)@df
  expect_gt(corr, naive)
  expect_lt(corr - naive, 5)         # a correction, not a second model
  expect_true(any(grepl("degrees of freedom carry",
                        summary(u, correct = TRUE)@notes)))
})

test_that("a term's edf is its share of the WHOLE model's smoother", {
  # The definition is the trace of the term's diagonal block of
  # F = (H + S)^-1 H over the coefficients of EVERY equation. A block-wise
  # tr[(H_bb + S_b)^-1 H_bb] drops the coupling with the rest of the model,
  # and the case that shows it is a parameter other than the mean: the
  # Demmler-Reinsch basis is orthogonalized against the constant in the
  # unweighted metric, so where the weights vary -- which they do in the
  # mean's equation as soon as the scale is modelled -- the orthogonality
  # the construction arranged does not survive.
  ref_total <- function(fit) {
    des <- statmod_design(fit@spec)
    H <- statmod_information_at(fit@spec, fit@coefficients, des, TRUE,
                                "bartlett")
    S <- statmod_penalty_at(fit@spec, fit@coefficients, fit@hyper, des,
                            "hessian")
    S[!is.finite(S)] <- 0
    list(d = diag(solve(H + S, H)), des = des)
  }
  per_term <- function(fit, r) {
    npar <- vapply(r$des, function(d) d$npar, integer(1))
    offs <- cumsum(npar) - npar
    out <- numeric(0)
    for (a in seq_along(fit@spec@distrib@params)) {
      p <- fit@spec@distrib@params[a]
      for (nm in names(fit@spec@terms[[p]])) {
        cols <- r$des[[p]]$blocks[[nm]]
        if (!length(cols)) next
        out[paste(p, nm)] <- sum(r$d[offs[a] + cols])
      }
    }
    out
  }

  set.seed(40)
  n <- 400L
  dd <- data.frame(x = stats::runif(n), z = stats::runif(n))
  dd$y <- sin(2 * pi * dd$x) +
    stats::rnorm(n, sd = exp(-1 + 0.8 * sin(2 * pi * dd$z)))

  cases <- list(
    statmod(y ~ s(x, k = 10) | sigma ~ s(z, k = 10),
            distributions7::gaussian1_distrib(), dd, outer_criterion = reml()),
    statmod(y ~ 1 | sigma ~ s(z, k = 10),
            distributions7::gaussian1_distrib(), dd, outer_criterion = reml())
  )
  for (fit in cases) {
    r <- ref_total(fit)
    got <- stats::setNames(fit@edf$edf, paste(fit@edf$parameter,
                                              fit@edf$term))
    want <- per_term(fit, r)
    expect_equal(got[names(want)], want, tolerance = 1e-10)
    expect_equal(sum(fit@edf$edf), sum(r$d), tolerance = 1e-10)
  }

  # a model with no penalty at all still counts one per coefficient
  d2 <- data.frame(x = stats::rnorm(80))
  d2$y <- 1 + 2 * d2$x + stats::rnorm(80)
  f2 <- statmod(y ~ x, distributions7::gaussian1_distrib(), d2)
  expect_equal(sum(f2@edf$edf), 3)

  # and a kinked penalty is still counted by its survivors, not by a trace
  # that has no curvature to read at the kink
  set.seed(23)
  p <- 12L
  X <- matrix(stats::rnorm(150 * p), 150, p)
  d3 <- data.frame(y = X[, 1] * 2 - X[, 2] * 1.5 + stats::rnorm(150),
                   X = I(X))
  f3 <- statmod(y ~ lasso(X), distributions7::gaussian1_distrib(), d3,
                hyper = list(mu = list("lasso(X)" = c(lambda = 20))))
  row <- f3@edf[f3@edf$term == "lasso(X)", ]
  expect_equal(row$edf, sum(coef(f3)$mu[-1] != 0))
})
