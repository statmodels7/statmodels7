# A term whose design block depends on its own coefficients, fitted by
# statmod(). The block is refreshed inside the objective, so the scoring
# step's increment is the Gauss-Newton one and the predictor is the term's
# contribution rather than its linearization.

set.seed(1)
n <- 300
dseg <- data.frame(x = sort(stats::runif(n, 0, 10)))
dseg$y <- 1 + 0.5 * dseg$x + 2 * pmax(dseg$x - 6, 0) + stats::rnorm(n, sd = 0.3)

# modelterms7's own iteration, which shares no code with the fitting layer:
# refresh the block, fit the working model, take the increment where the
# block is a Jacobian and the coefficients where it is not
seg_iterate <- function(built, y, iters = 60) {
  b <- built@blueprint$coef
  cur <- built
  for (it in seq_len(iters)) {
    cur <- modelterms7::term_refresh(cur, b)
    X <- modelterms7::term_matrix(cur)
    if (cur@kind == "seg") {
      b <- b + as.numeric(qr.coef(qr(X), y - modelterms7::term_value(cur)))
    } else {
      b <- as.numeric(qr.coef(qr(X), y))
    }
    if (modelterms7::seg_converged(cur)) break
  }
  list(coef = b, psi = modelterms7::seg_psi(cur), iterations = it)
}

fitted_term <- function(fit, what) {
  nm <- grep(what, names(fit@spec@terms$mu), fixed = TRUE, value = TRUE)[1L]
  fit@spec@terms$mu[[nm]]
}

test_that("statmod finds a break-point, from every start, and agrees", {
  ref <- seg_iterate(modelterms7::term_build(modelterms7::seg(x, psi = 5),
                                             dseg), dseg$y - 1)
  psis <- numeric(0)
  rss <- numeric(0)
  for (p0 in c(2, 3, 5, 8)) {
    fm <- stats::as.formula(sprintf("y ~ seg(x, psi = %g)", p0))
    fit <- statmod(fm, distributions7::gaussian1_distrib(), dseg)
    expect_true(fit@converged)
    psis <- c(psis, modelterms7::seg_psi(fitted_term(fit, "seg(")))
    rss <- c(rss, sum((dseg$y - fit@fitted$mu)^2))
  }
  # the same answer whatever the start, the truth being 6
  expect_equal(psis, rep(psis[1L], length(psis)), tolerance = 1e-6)
  expect_equal(rss, rep(rss[1L], length(rss)), tolerance = 1e-6)
  expect_equal(psis[1L], 6, tolerance = 0.1)
  # and the same answer the term's own iteration reaches
  expect_equal(psis[1L], ref$psi, tolerance = 0.01)
})

test_that("a continuous construction is fitted without a step in it", {
  # the defect this replaces: with the block frozen, the psi column became an
  # ordinary regressor and the fitted mean of a CONTINUOUS term carried a
  # step, while the run reported convergence
  fit <- statmod(y ~ seg(x, psi = 3), distributions7::gaussian1_distrib(),
                 dseg)
  tm <- fitted_term(fit, "seg(")
  psi <- modelterms7::seg_psi(tm)
  mu <- fit@fitted$mu

  # the fit either side of the break-point, over the smallest gap in x that
  # straddles it: a step would show as a slope far outside the two segments'
  lo <- max(dseg$x[dseg$x < psi])
  hi <- min(dseg$x[dseg$x > psi])
  slope <- (mu[dseg$x == hi] - mu[dseg$x == lo]) / (hi - lo)
  b <- fit@coefficients$mu
  expect_gt(slope, min(b[2L], b[2L] + b[3L]) - 1e-6)
  expect_lt(slope, max(b[2L], b[2L] + b[3L]) + 1e-6)

  # and the contribution is the segmented function itself, not the block
  # times the coefficients
  val <- modelterms7::term_value(tm, coef = b[-1L])
  lin <- as.numeric(modelterms7::term_matrix(tm) %*% b[-1L])
  expect_equal(mu, b[1L] + val, tolerance = 1e-10)
  expect_gt(max(abs(val - lin)), 1)
})

test_that("a discontinuous term is fitted and stopped on its own rule", {
  set.seed(4)
  dj <- data.frame(x = sort(stats::runif(400, 0, 10)))
  dj$y <- 3 * (dj$x > 6.5) + stats::rnorm(400, sd = 0.4)
  got <- list()
  for (p0 in c(3, 5, 8)) {
    fm <- stats::as.formula(sprintf("y ~ jump(x, psi = %g)",
                                    p0))
    fit <- statmod(fm, distributions7::gaussian1_distrib(), dj)
    expect_true(fit@converged)
    got[[length(got) + 1L]] <- c(
      psi = modelterms7::seg_psi(fitted_term(fit, "jump(")),
      kappa = unname(fit@coefficients$mu[2L]))
  }
  g <- do.call(rbind, got)
  expect_equal(g[, "psi"], rep(6.5, 3), tolerance = 0.1)
  expect_equal(g[, "kappa"], rep(3, 3), tolerance = 0.2)

  # the verdict is the term's, not the score's. The block of a discontinuous
  # construction is a linearization with a frozen weight, so the gradient it
  # gives belongs to the working model: at the answer it sits at 0.18 per
  # observation and never vanishes, while the break-point has stopped moving
  fit <- statmod(y ~ jump(x, psi = 5),
                 distributions7::gaussian1_distrib(), dj)
  expect_true(modelterms7::term_converged(fitted_term(fit, "jump(")))
  expect_true(modelterms7::term_converged(
    modelterms7::term_build(modelterms7::linpar(~x), dj)))
})

test_that("nl() reaches nls's answer", {
  set.seed(7)
  dn <- data.frame(x = seq(0, 3, length.out = 120))
  dn$y <- 2 * exp(-1.3 * dn$x) + stats::rnorm(120, sd = 0.05)
  fit <- statmod(y ~ nl(~ a * exp(-r * x), start = list(a = 1, r = 1)) - 1,
                 distributions7::gaussian1_distrib(), dn)
  expect_true(fit@converged)
  ref <- stats::nls(y ~ a * exp(-r * x), dn, start = list(a = 1, r = 1))
  # nls's own convergence tolerance is 1e-8 relative, so that is the
  # agreement worth asserting
  expect_equal(unname(fit@coefficients$mu), unname(stats::coef(ref)),
               tolerance = 1e-5)
})

test_that("prediction reapplies the fitted term rather than rebuilding it", {
  fit <- statmod(y ~ seg(x, psi = 5), distributions7::gaussian1_distrib(),
                 dseg)
  rows <- c(3L, 40L, 150L, 220L, 299L)
  # the identity that pins it: predicting on rows the model was fitted to
  # returns the fitted values there. The contribution is not the block times
  # the coefficients, so a predictor reading only the block fails this
  expect_equal(predict(fit, "mu", newdata = dseg[rows, , drop = FALSE]),
               fit@fitted$mu[rows], tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(predict(fit, "mu", newdata = dseg), fit@fitted$mu,
               tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("a break-point term composes with a penalized one", {
  set.seed(11)
  dm <- data.frame(x = sort(stats::runif(250, 0, 10)),
                   z1 = stats::rnorm(250), z2 = stats::rnorm(250))
  dm$y <- 1 + 0.5 * dm$x + 2 * pmax(dm$x - 6, 0) + 0.8 * dm$z1 +
    stats::rnorm(250, sd = 0.3)
  fit <- statmod(y ~ seg(x, psi = 5) + ridge(~ z1 + z2),
                 distributions7::gaussian1_distrib(), dm,
                 outer_criterion = reml())
  expect_true(fit@converged)
  expect_equal(modelterms7::seg_psi(fitted_term(fit, "seg(")), 6,
               tolerance = 0.3)
  # the hyperparameter of the ridge is still estimated, and the term with
  # the signal is not shrunk away
  expect_true(is.finite(fit@hyper$mu[["ridge(~z1 + z2)"]][["sigma"]]))
  des <- statmod_design(fit@spec)
  cols <- des$mu$blocks[["ridge(~z1 + z2)"]]
  z1 <- fit@coefficients$mu[cols][grep("z1", des$mu$coef_names[cols])]
  expect_gt(abs(z1), 0.4)
})
