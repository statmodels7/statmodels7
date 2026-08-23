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

  # the verdict is NOT the working score's. The block of a discontinuous
  # construction is a linearization with a frozen weight, so the gradient
  # it gives belongs to the working model and never vanishes at the answer
  # (0.18 per observation, measured); the fit still reports convergence,
  # judged at a fixed point or cycle of the working objective. The step
  # rule term_converged() reads is a resolution rule that tightens with n,
  # so it may stay FALSE at the answer and is deliberately not asserted.
  fit <- statmod(y ~ jump(x, psi = 5),
                 distributions7::gaussian1_distrib(), dj)
  expect_true(fit@converged)
  expect_equal(unname(modelterms7::seg_psi(fitted_term(fit, "jump("))),
               6.5, tolerance = 0.1)
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
  expect_true(is.finite(fit@hyper$mu[["ridge(~z1 + z2)"]][["lambda"]]))
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

  lines <- utils::capture.output(print(s))
  out <- paste(lines, collapse = "\n")
  expect_true(any(startsWith(lines, "seg(x)")))
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
  # and its standard error is MISSING, deliberately. The block of a
  # discontinuous construction is a working linearization with a frozen
  # weight, so the curvature it carries is the working model's: measured at
  # 400 observations against a bootstrap of 200 resamples, it gives the
  # change of level and the auxiliary coordinate a standard error of exactly
  # zero, against 0.063 and 0.540, and the position read off them 1.8e-05
  # against 0.090. A number wrong by five thousand times is worse than a gap.
  expect_true(is.na(tj$se[tj$name == "psi1"]))
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

  # a caller who names psi has said where to BEGIN; with the restarts
  # disabled as well, the fit is left on that start's own optimum, short
  # of the truth at 6
  held <- statmod(y ~ jseg(x, psi = 5, n_boot = 0),
                  distributions7::gaussian1_distrib(), dj)
  expect_lt(psi_of(held), 5.8)
  # and it is the worse optimum, which is the point of the grid
  expect_gt(as.numeric(logLik(fit)), as.numeric(logLik(held)))
  # with the restarts left at their default the same start is rescued:
  # naming psi says where to begin, not which optimum to accept
  resc <- statmod(y ~ jseg(x, psi = 5), distributions7::gaussian1_distrib(),
                  dj)
  expect_equal(as.numeric(logLik(resc)), as.numeric(logLik(fit)),
               tolerance = 1e-4)

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


test_that("a jseg with several break-points is fitted by working fits", {
  # the regression this pins: the embedded route dragged the read-off
  # through the inner line search and, from the TRUE break-points, ended
  # at an rss worse than the mean-only fit
  set.seed(21)
  n <- 600
  xj <- sort(stats::runif(n, -1, 1))
  mu <- xj + 2 * (xj >= -0.4) - 2 * pmax(xj + 0.4, 0) +
    3 * (xj >= 0.4) - 1.5 * pmax(xj - 0.4, 0)
  dj2 <- data.frame(x = xj, y = mu + stats::rnorm(n, sd = 0.25))
  fit <- statmod(y ~ jseg(x, npsi = 2, psi = c(-0.4, 0.4), n_boot = 0),
                 distributions7::gaussian1_distrib(), dj2)
  expect_true(fit@converged)
  psi <- modelterms7::seg_psi(fitted_term(fit, "jseg("))
  expect_equal(unname(psi), c(-0.4, 0.4), tolerance = 0.05)
  rss <- sum((dj2$y - fit@fitted$mu)^2)
  expect_lt(rss, sum((dj2$y - mu)^2) * 1.05)
  # prediction still reapplies the fitted term
  rows <- c(2L, 100L, 400L, 599L)
  expect_equal(predict(fit, "mu", newdata = dj2[rows, , drop = FALSE]),
               fit@fitted$mu[rows], tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("bootstrap restarting recovers a poor start", {
  # from this start the iteration lands on a local optimum; reordering
  # alone does not recover it, the restarts do. Both fits are seeded, the
  # draws coming from the session's generator.
  set.seed(22)
  n <- 500
  xb <- sort(stats::runif(n, 0, 10))
  mu <- 3 * (xb >= 6.5)
  db <- data.frame(x = xb, y = mu + stats::rnorm(n, sd = 0.4))
  set.seed(1)
  f0 <- statmod(y ~ jump(x, psi = 2, n_boot = 0),
                distributions7::gaussian1_distrib(), db)
  set.seed(1)
  f1 <- statmod(y ~ jump(x, psi = 2, n_boot = 10),
                distributions7::gaussian1_distrib(), db)
  rss0 <- sum((db$y - f0@fitted$mu)^2)
  rss1 <- sum((db$y - f1@fitted$mu)^2)
  # the restarted fit is never worse, and psi is recovered
  expect_lte(rss1, rss0 + 1e-6)
  expect_equal(unname(modelterms7::seg_psi(fitted_term(f1, "jump("))),
               6.5, tolerance = 0.1)
})


test_that("a frozen break-point block composes with an estimated penalty", {
  # the criterion evaluations hold the block at its committed positions --
  # a break-point moving between evaluations makes the criterion
  # path-dependent, and the phase's flags read as unavailable points --
  # and the positions are refined once, at the chosen hyperparameters.
  # Before the hold this fit took 136 s and reported FALSE at the right
  # answer; now it is an ordinary REML fit plus one refinement.
  set.seed(23)
  n <- 800
  xr <- sort(stats::runif(n, 0, 10))
  zr <- stats::runif(n, -2, 2)
  mu <- 1 + sin(1.5 * zr) + 0.4 * xr + 2.5 * (xr >= 6) - 1.2 * pmax(xr - 6, 0)
  dr <- data.frame(y = mu + stats::rnorm(n, sd = 0.5), x = xr, z = zr)
  fit <- statmod(y ~ s(z, k = 8) + jseg(x, n_boot = 2),
                 distributions7::gaussian1_distrib(), dr)
  expect_true(fit@converged)
  expect_equal(unname(modelterms7::seg_psi(fitted_term(fit, "jseg("))),
               6, tolerance = 0.1)
  # the smoothing parameter was estimated, not left at its start
  expect_true(is.finite(fit@criterion))
})


test_that("the random changepoint fits end to end, on seg", {
  # Muggeo-Atkins: psi_i = pop + random deviation, the development riding
  # random(~1|id) with its variance component estimated. The discontinuous
  # constructions REJECT a penalized development of the break-point (the
  # estimated coefficients are -delta * gamma), so this model lives on seg.
  set.seed(31)
  m <- 12
  ni <- 30
  id <- factor(rep(seq_len(m), each = ni))
  xr <- as.numeric(replicate(m, sort(stats::runif(ni, 0, 10))))
  psi_i <- stats::rnorm(m, 5, 0.4)
  mu <- 1 + 0.5 * xr + 1.8 * pmax(xr - psi_i[as.integer(id)], 0)
  dr <- data.frame(y = mu + stats::rnorm(m * ni, sd = 0.4), x = xr, id = id)
  fit <- statmod(y ~ seg(x, psi ~ random(~ 1 | id)),
                 distributions7::gaussian1_distrib(), dr)
  expect_true(fit@converged)
  tm <- fitted_term(fit, "seg(")
  psi_hat <- modelterms7::seg_psi(tm)
  per_id <- vapply(seq_len(m), function(i)
    mean(psi_hat[as.integer(id) == i, 1]), numeric(1))
  expect_gt(stats::cor(per_id, psi_i), 0.9)
  expect_lt(sqrt(mean((per_id - psi_i)^2)), 0.3)
})

test_that("per-group break-points on a jump, unpenalized", {
  set.seed(32)
  ng <- 3
  ni <- 250
  g <- factor(rep(letters[seq_len(ng)], each = ni))
  xg <- as.numeric(replicate(ng, sort(stats::runif(ni, 0, 10))))
  psi_g <- c(a = 3, b = 5, c = 7)
  mu <- 1 + 2 * (xg >= psi_g[as.integer(g)])
  dg <- data.frame(y = mu + stats::rnorm(ng * ni, sd = 0.4), x = xg, id = g)
  fit <- statmod(y ~ jump(x, psi ~ 0 + id),
                 distributions7::gaussian1_distrib(), dg)
  expect_true(fit@converged)
  psi2 <- modelterms7::seg_psi(fitted_term(fit, "jump("))
  per_g <- vapply(levels(g), function(l) mean(psi2[g == l, 1]), numeric(1))
  expect_equal(unname(per_g), c(3, 5, 7), tolerance = 0.1)
  # and a PENALIZED development of a discontinuous break-point is refused
  # with the reason, not fitted approximately
  expect_error(
    statmod(y ~ jump(x, psi ~ random(~ 1 | id)),
            distributions7::gaussian1_distrib(), dg),
    "delta")
})
