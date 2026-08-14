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

test_that("a break-point term has its own section, reporting psi", {
  # The block is a working linearization, so its coefficients are not the
  # quantities of the model: a reader wants the changes and the positions,
  # not the auxiliary pair a position is read off.
  set.seed(41)
  n <- 400
  dd <- data.frame(x = sort(stats::runif(n, 0, 10)))
  dd$y <- 1 + 0.4 * dd$x + 2 * pmax(dd$x - 6, 0) + stats::rnorm(n, sd = 0.3)
  fit <- statmod(y ~ seg(x), distributions7::gaussian1_distrib(), dd)
  s <- summary(fit)

  blk <- Filter(function(b) identical(b$kind, "breakpoint"), s@tables$mu)
  expect_length(blk, 1L)
  tb <- blk[[1L]]$table
  expect_identical(tb$name, c("beta", "gamma1", "psi1"))
  # the position agrees with the term's own report, and is not a
  # coefficient's name
  nm <- grep("seg(", names(fit@spec@terms$mu), fixed = TRUE, value = TRUE)[1L]
  expect_equal(tb$estimate[tb$name == "psi1"],
               as.numeric(modelterms7::seg_psi(fit@spec@terms$mu[[nm]])),
               tolerance = 1e-10)
  expect_equal(tb$estimate[tb$name == "psi1"], 6, tolerance = 0.1)
  expect_true(is.finite(tb$se[tb$name == "psi1"]))
  # a break-point gets an interval and no test: the null a z would report
  # on is that the position is zero, and under the null a reader cares
  # about the position is a nuisance parameter that vanishes
  expect_true(is.na(tb$statistic[tb$name == "psi1"]))
  expect_true(is.na(tb$p_value[tb$name == "psi1"]))
  expect_true(all(is.finite(unlist(tb[tb$name == "psi1", c("lower", "upper")]))))

  out <- paste(utils::capture.output(print(s)), collapse = "\n")
  expect_match(out, "Break-points")
  expect_match(out, "psi1")
  # and the auxiliary pair is nowhere in it
  expect_false(grepl("seg.psi1", out, fixed = TRUE))

  # a discontinuous term reports the position it holds, read off two
  # coefficients, and its standard error by the delta method
  dj <- dd
  dj$y <- 1 + 2 * (dj$x > 6) + stats::rnorm(n, sd = 0.3)
  fj <- statmod(y ~ jump(x), distributions7::gaussian1_distrib(), dj)
  tj <- Filter(function(b) identical(b$kind, "breakpoint"),
               summary(fj)@tables$mu)[[1L]]$table
  expect_identical(tj$name, c("delta1", "psi1"))
  expect_equal(tj$estimate[tj$name == "psi1"], 6, tolerance = 0.1)
  expect_true(is.finite(tj$se[tj$name == "psi1"]))
  expect_false(any(grepl("^g[0-9]", tj$name)))
})

test_that("a break-point term starts where it asks to, not at zero", {
  # zero is degenerate for a jump: psi = -g/delta is 0/0, every
  # break-point lands on the same clamped position and the block is
  # singular, so the fit used to return every coefficient exactly zero
  set.seed(42)
  y <- stats::rbinom(200, 1, rep(c(0.1, 0.7, 0.3, 0.9), each = 50))
  dd <- data.frame(x = seq_len(200) / 200, y = y)
  fit <- statmod(y ~ jump(x, npsi = 3), distributions7::bernoulli_distrib(),
                 dd)
  expect_true(fit@converged)
  expect_false(all(fit@coefficients$mu[-1L] == 0))
  nm <- grep("jump(", names(fit@spec@terms$mu), fixed = TRUE, value = TRUE)[1L]
  psi <- as.numeric(modelterms7::seg_psi(fit@spec@terms$mu[[nm]]))
  expect_equal(sort(psi) * 200, c(50, 100, 150), tolerance = 6)

  # asked for explicitly, the origin is still the origin
  z <- statmod_spec(y ~ jump(x, npsi = 3),
                    distributions7::bernoulli_distrib(), dd)
  d <- statmod_design(z)
  expect_true(all(start_at(start_origin(), z, d, NULL)$mu == 0))
})

test_that("a break-point term whose psi is not named starts on a grid", {
  # The objective has local optima in the break-point and the iteration
  # converges from within a basin, so the term's own default -- the
  # interior quantiles of the covariate, which look at the covariate and
  # not at the response -- is a conventional start and reaches one of
  # them. From the grid this model recovers the truth; from the median it
  # converges to a genuine local minimum with the jump the WRONG SIGN.
  set.seed(2)
  n <- 400
  dj <- data.frame(x = sort(stats::runif(n, 0, 10)))
  dj$y <- 1 + 0.4 * dj$x + 2 * pmax(dj$x - 6, 0) + 1.5 * (dj$x > 6) +
    stats::rnorm(n, sd = 0.3)

  psi_of <- function(fit) {
    nm <- grep("jseg(", names(fit@spec@terms$mu), fixed = TRUE,
               value = TRUE)[1L]
    as.numeric(modelterms7::seg_psi(fit@spec@terms$mu[[nm]]))
  }
  fit <- statmod(y ~ jseg(x), distributions7::gaussian1_distrib(), dj)
  expect_equal(psi_of(fit), 6, tolerance = 0.15)
  # the coefficients are a bare vector: the names are the design's
  cf <- stats::setNames(fit@coefficients$mu,
                        statmod_design(fit@spec)$mu$coef_names)
  expect_gt(cf[["jseg.delta1"]], 0)
  expect_equal(unname(cf[["jseg.gamma1"]]), 2, tolerance = 0.2)

  # a caller who names psi has said where to begin, and is left there
  held <- statmod(y ~ jseg(x, psi = 5), distributions7::gaussian1_distrib(),
                  dj)
  expect_lt(abs(psi_of(held) - 5), 0.2)
  # and it is the worse optimum, which is the point of the grid
  expect_gt(as.numeric(logLik(fit)), as.numeric(logLik(held)))

  # the rule reads the response, so it is skipped where the response is
  # not plain numbers rather than being given a reading of its own
  dc <- dj
  dc$ev <- rep(1, n)
  spec <- modelterms7::jseg(x)
  expect_identical(seg_grid_start(spec, dc,
                                  modelterms7::cens(dc$y, dc$ev)), spec)
  expect_identical(seg_grid_start(modelterms7::linpar(~x), dc, dj$y),
                   modelterms7::linpar(~x))
  # and a named psi turns it off at the source
  named <- modelterms7::jseg(x, psi = 5)
  expect_identical(seg_grid_start(named, dj, dj$y), named)
})
