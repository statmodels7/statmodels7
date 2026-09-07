## The restricted fit and the four likelihood statistics.
##
## Three of the four -- the likelihood ratio, Rao's score and Terrell's
## gradient -- are functions of the model refitted with one coefficient HELD
## at a value, so everything here rests on that refit being the constrained
## maximum and not merely a fit started from the constraint.  Every check
## below is against a route that shares no arithmetic with ours: a glm() with
## the held term as an OFFSET, which is the same restricted problem written
## another way; glm()'s own deviance difference for the likelihood ratio; the
## score statistic's general form U'K^-1 U, of which the package computes only
## the one product it reduces to; and, for the inverted interval, a profile
## bisection written here over glm()'s deviance.
##
## The profile interval another package would supply is not used: an external
## package is not a dependency of this toolkit even in Suggests, so the
## reference is implemented here instead.

pois_data <- function(n = 300, seed = 4) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n, -1, 1), z = stats::runif(n, -1, 1))
  d$y <- stats::rpois(n, exp(0.4 + 0.8 * d$x - 0.5 * d$z))
  ## drawn LAST so that x, z and y are the same numbers they always were
  d$g <- factor(sample(20L, n, replace = TRUE))
  d
}

pois_fit <- function(d) {
  statmod(y ~ x + z, distributions7::poisson_distrib(), d)
}

## a smooth wants data with a curve in it. Fitted to a LINEAR truth the
## smoothing parameter correctly runs to 3.8e+08, every rotated coordinate
## sits at 1e-08 and the block is a straight line: a degenerate object to
## read inference off, and the place where the inner flag reports FALSE at a
## mode located to 3.6e-10. That regime is a control below, not the setting.
smooth_data <- function(n = 300L, seed = 4) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n, -1, 1), z = stats::runif(n, -1, 1))
  d$y <- stats::rpois(n, exp(0.4 + 0.8 * d$x + sin(3 * d$z)))
  d
}

smooth_fit <- function(d = smooth_data()) {
  statmod(y ~ x + s(z, k = 8), distributions7::poisson_distrib(), d)
}

lasso_fit <- function(m = 200L, seed = 7) {
  set.seed(seed)
  X <- matrix(stats::rnorm(m * 8L), m, 8L)
  d <- data.frame(y = stats::rpois(m, exp(0.3 + X %*% c(1, -0.8, rep(0, 6)))))
  d$X <- X
  statmod(y ~ lasso(X, lambda = 2), distributions7::poisson_distrib(), d)
}


test_that("a restricted fit holds one coefficient and maximizes over the rest", {
  d <- pois_data()
  fit <- pois_fit(d)
  nms <- statmod_design(fit@spec)$mu$coef_names
  bx <- fit@coefficients$mu[[match("x", nms)]]

  ## held AT the unrestricted value the constraint binds nothing, so the
  ## refit must reproduce the fit rather than merely come close to it
  r0 <- statmod_restrict(fit, "mu", "x", bx)
  expect_identical(r0$coefficients$mu[[2L]], bx)
  expect_equal(r0$loglik, as.numeric(logLik(fit)), tolerance = 1e-8)
  expect_equal(unlist(r0$coefficients), unlist(fit@coefficients),
               tolerance = 1e-6)

  ## held AWAY from it the free coordinates are at their own optimum and the
  ## held one is not: a score of zero everywhere would mean the hold had not
  ## been enforced, and a score of zero at the held coordinate too would mean
  ## the point is the unrestricted maximum
  r <- statmod_restrict(fit, "mu", "x", 0.4)
  expect_true(r$converged)
  expect_lt(max(abs(r$score[-r$at])), 1e-5)
  expect_gt(abs(r$score[[r$at]]), 1)
  expect_lt(r$loglik, as.numeric(logLik(fit)))
})


test_that("the restricted fit is the same problem glm() solves with an offset", {
  d <- pois_data()
  fit <- pois_fit(d)
  for (b in c(0, 0.4, 1.2)) {
    r <- statmod_restrict(fit, "mu", "x", b)
    g <- stats::glm(y ~ z + offset(b * x), family = stats::poisson(),
                    data = d)
    expect_equal(r$loglik, as.numeric(stats::logLik(g)), tolerance = 1e-8)
    expect_equal(c(r$coefficients$mu[[1L]], r$coefficients$mu[[3L]]),
                 unname(stats::coef(g)), tolerance = 1e-6)
  }
})


test_that("the four statistics agree with the references that have one", {
  d <- pois_data()
  fit <- pois_fit(d)
  g <- stats::glm(y ~ x + z, family = stats::poisson(), data = d)

  ## Wald against glm's own z^2
  w <- statmod_stat_at(fit, "mu", "x", 0, "wald")
  expect_equal(w$statistic,
               (summary(g)$coefficients["x", "z value"])^2,
               tolerance = 1e-4)

  ## the likelihood ratio against the deviance difference, which is the same
  ## quantity assembled from two unpenalized fits
  g0 <- stats::glm(y ~ z, family = stats::poisson(), data = d)
  lr <- statmod_stat_at(fit, "mu", "x", 0, "lr")
  expect_equal(lr$statistic,
               as.numeric(stats::deviance(g0) - stats::deviance(g)),
               tolerance = 1e-6)

  ## the score statistic's GENERAL form, U'K^-1 U over every coordinate,
  ## against the single product the package computes. They coincide only
  ## because every free component of U vanishes at the restricted maximum,
  ## so the agreement is also a statement that the refit converged.
  for (b in c(0, 0.4, 1.2)) {
    r <- statmod_restrict(fit, "mu", "x", b)
    V <- solve(r$information())
    full <- as.numeric(t(r$score) %*% V %*% r$score)
    s <- statmod_stat_at(fit, "mu", "x", b, "score")
    expect_equal(s$statistic, full, tolerance = 1e-6)
  }

  ## every p-value is the chi-squared one on one degree of freedom
  for (t in c("wald", "lr", "score", "gradient")) {
    s <- statmod_stat_at(fit, "mu", "x", 0, t)
    expect_identical(s$df, 1L)
    expect_equal(s$p.value,
                 stats::pchisq(s$statistic, 1L, lower.tail = FALSE))
  }
})


test_that("at the estimate itself all four statistics vanish", {
  d <- pois_data()
  fit <- pois_fit(d)
  nms <- statmod_design(fit@spec)$mu$coef_names
  bx <- fit@coefficients$mu[[match("x", nms)]]
  for (t in c("wald", "lr", "score", "gradient")) {
    expect_lt(abs(statmod_stat_at(fit, "mu", "x", bx, t)$statistic), 1e-6)
  }
})


test_that("inverting the Wald statistic returns the closed-form interval", {
  d <- pois_data()
  fit <- pois_fit(d)
  nms <- statmod_design(fit@spec)$mu$coef_names
  bx <- fit@coefficients$mu[[match("x", nms)]]
  se <- sqrt(stats::vcov(fit, readable = FALSE)[2L, 2L])
  zq <- stats::qnorm(0.975)
  ci <- statmod_invert(fit, "mu", "x", 0.95, "wald")
  ## the Wald statistic is a parabola in b, so the search must land on the
  ## two points the closed form gives and the interval must be symmetric
  expect_equal(as.numeric(ci), c(bx - zq * se, bx + zq * se),
               tolerance = 1e-7)
  expect_equal((ci[["upper"]] - bx) - (bx - ci[["lower"]]), 0,
               tolerance = 1e-9)
})


test_that("the likelihood-ratio interval is the profile one", {
  d <- pois_data()
  fit <- pois_fit(d)
  ## the reference: bisect the deviance difference of glm() with the held
  ## term as an offset. It shares no arithmetic with statmod_invert(), which
  ## brackets from the Wald interval and inverts our own restricted fit.
  g <- stats::glm(y ~ x + z, family = stats::poisson(), data = d)
  q <- stats::qchisq(0.95, 1L)
  gap <- function(b) {
    gb <- stats::glm(y ~ z + offset(b * x), family = stats::poisson(),
                     data = d)
    (stats::deviance(gb) - stats::deviance(g)) - q
  }
  bx <- unname(stats::coef(g)[["x"]])
  se <- sqrt(stats::vcov(g)[2L, 2L])
  ref <- c(stats::uniroot(gap, c(bx - 6 * se, bx), tol = 1e-10)$root,
           stats::uniroot(gap, c(bx, bx + 6 * se), tol = 1e-10)$root)

  ci <- statmod_invert(fit, "mu", "x", 0.95, "lr")
  expect_equal(as.numeric(ci), ref, tolerance = 1e-5)

  ## and it is NOT the Wald interval: the whole reason to invert a
  ## likelihood is that the two differ where the likelihood is not a
  ## parabola, so a check that could be satisfied by returning the Wald
  ## ends would say nothing
  w <- statmod_invert(fit, "mu", "x", 0.95, "wald")
  expect_gt(max(abs(as.numeric(ci) - as.numeric(w))), 1e-5)
})


test_that("the three restricted intervals are not symmetric and Wald's is", {
  ## small n and a steep predictor, where the likelihood is visibly skewed
  set.seed(11)
  m <- 25
  d <- data.frame(x = stats::runif(m, -1, 1))
  d$y <- stats::rpois(m, exp(0.2 + 1.4 * d$x))
  fit <- statmod(y ~ x, distributions7::poisson_distrib(), d)
  b <- fit@coefficients$mu[[2L]]
  asym <- function(t) {
    ci <- statmod_invert(fit, "mu", "x", 0.95, t)
    (ci[["upper"]] - b) - (b - ci[["lower"]])
  }
  expect_lt(abs(asym("wald")), 1e-9)
  expect_gt(asym("lr"), 1e-3)
  expect_gt(asym("gradient"), 1e-3)
})


test_that("what cannot be held is refused, and what can is not", {
  d <- pois_data()
  fit <- pois_fit(d)

  expect_error(statmod_restrict(fit, "sigma", "x", 0), "not a parameter")
  expect_error(statmod_restrict(fit, "mu", "q", 0), "names no coefficient")

  ## A penalty that is twice differentiable does NOT put a coordinate out of
  ## reach: the hold is enforced in the solve, and what the statistic then
  ## means is the credible reading vcov(type = "bayesian") already gives that
  ## row. The smooth's linear column is the case that motivated it -- s() is
  ## penalized by diag(0, 1, ..., 1), so `lin` is not shrunk at all.
  fs <- smooth_fit()
  oks <- names(testable_coords(fs))
  expect_true(all(c("mu:(Intercept)", "mu:x", "mu:s(z).lin", "mu:s(z).z1")
                  %in% oks))
  r <- statmod_restrict(fs, "mu", "s(z).lin", 0)
  expect_identical(r$coefficients$mu[[3L]], 0)
  ## the free coordinates are at their own optimum, which is the property
  ## that matters; the flag is the inner rule's and reports FALSE at an
  ## extreme smoothing parameter even where the mode is located to 1e-10
  expect_lt(max(abs(r$score[-r$at])), 1e-4)
  expect_gt(abs(r$score[[r$at]]), 1)

  ## and a random effect's own coefficients likewise
  fr <- statmod(y ~ x + random(~1 | g), distributions7::poisson_distrib(), d)
  nr <- statmod_design(fr@spec)$mu$coef_names
  expect_true(paste0("mu:", nr[[3L]]) %in% names(testable_coords(fr)))

  ## a KINKED penalty is the one that is refused, and for two reasons that
  ## point the same way -- see statmod_restrict()'s page
  fl <- lasso_fit()
  nl <- statmod_design(fl@spec)$mu$coef_names
  expect_error(statmod_restrict(fl, "mu", nl[[2L]], 0), "kinked penalty")
  expect_identical(names(testable_coords(fl)), "mu:(Intercept)")

  ## testable_coords says the same thing once, for the whole vector
  ok <- testable_coords(fit)
  expect_identical(unname(ok), seq_len(3L))
  expect_identical(names(ok), c("mu:(Intercept)", "mu:x", "mu:z"))

  ## a model whose contribution is a recursion's state and not X beta
  set.seed(3)
  n <- 200
  dg <- data.frame(t = seq_len(n))
  dg$y <- as.numeric(stats::arima.sim(list(ar = 0.5), n))
  fg <- statmod(y ~ gas(p = 1, q = 1),
                distributions7::gaussian1_distrib(), dg)
  expect_identical(unname(testable_coords(fg)), integer(0))
  expect_error(statmod_restrict(fg, "mu", "(Intercept)", 0), "structural")
})


test_that("the coordinate descent refuses a hold instead of ignoring it", {
  ## statmod_restrict() turns such a coordinate away before it reaches here,
  ## so this is what stops that refusal from being the only guard. Measured
  ## before it existed: held at 3, the sweep returned 2.7348.
  fl <- lasso_fit()
  nl <- statmod_design(fl@spec)$mu$coef_names
  spec <- S7::set_props(fl@spec,
                        held_coef = list(mu = stats::setNames(3, nl[[2L]])))
  des <- statmod_design(spec)
  blk <- statmod_blocks(spec, des)
  obj <- statmod_objective(spec, fl@hyper, des, TRUE, "opg")
  beta <- unlist(fl@coefficients, use.names = FALSE)
  expect_error(
    coord_fit(obj, beta, blk$sparse[[1L]], fl@hyper, spec, des, TRUE, "opg"),
    "cannot leave one where it was")

  ## and a block with nothing held runs as it did
  obj0 <- statmod_objective(fl@spec, fl@hyper, statmod_design(fl@spec),
                            TRUE, "opg")
  blk0 <- statmod_blocks(fl@spec, statmod_design(fl@spec))
  expect_silent(coord_fit(obj0, beta, blk0$sparse[[1L]], fl@hyper, fl@spec,
                          statmod_design(fl@spec), TRUE, "opg"))
})


test_that("a penalized coordinate is tested on the penalized objective", {
  fs <- smooth_fit()
  nms <- statmod_design(fs@spec)$mu$coef_names

  ## THE STATISTIC IS NON-NEGATIVE BY CONSTRUCTION, which reading the
  ## log-likelihood alone does not give: the restricted point maximizes the
  ## penalized objective, so it can sit above the unrestricted point on the
  ## likelihood and the difference come out negative. Measured on `s(z).z3`
  ## at -0.0379 and -0.1075 before the change.
  for (nm in c("x", "s(z).lin", "s(z).z1", "s(z).z3")) {
    b0 <- fs@coefficients$mu[[match(nm, nms)]]
    v <- b0 + c(-0.5, -0.1, 0.1, 0.5) * max(abs(b0), 0.2)
    st <- vapply(v, function(b)
      statmod_stat_at(fs, "mu", nm, b, "lr")$statistic, numeric(1))
    expect_true(all(st >= 0))
    ## and exactly zero at the estimate
    expect_lt(abs(statmod_stat_at(fs, "mu", nm, b0, "lr")$statistic), 1e-6)
  }

  ## THE CHECK that the reading is the one claimed: the interval by inversion
  ## and the Wald interval from the bayesian variance are the same to second
  ## order, so on a nearly quadratic objective they nearly coincide. A route
  ## reading the wrong objective would not.
  V <- stats::vcov(fs, readable = FALSE)
  zq <- stats::qnorm(0.975)
  for (nm in c("s(z).lin", "s(z).z1")) {
    key <- paste0("mu:", nm)
    b0 <- fs@coefficients$mu[[match(nm, nms)]]
    se <- sqrt(V[key, key])
    ci <- confint(fs, key, test = "lr", readable = FALSE)
    expect_equal(c(ci$lower[[1L]], ci$upper[[1L]]),
                 c(b0 - zq * se, b0 + zq * se), tolerance = 0.02,
                 ignore_attr = TRUE)
    ## and the inversion is not merely returning them: it is asymmetric
    expect_gt(abs((ci$upper[[1L]] - b0) - (b0 - ci$lower[[1L]])), 1e-4)
  }

  ## THE OPPOSITE LIMIT, which is the same statement read the other way: at
  ## an extreme smoothing parameter the objective in a rotated coordinate IS
  ## the quadratic penalty, so the inverted interval is exactly the Wald one
  ## and exactly symmetric. Measured, `s(z).z1` sits at 1.2e-08 there.
  fl <- statmod(y ~ x + s(z, k = 8), distributions7::poisson_distrib(),
                pois_data())
  expect_gt(fl@hyper$mu[["s(z, k = 8)"]][["lambda"]], 1e6)
  ## its own names: the two models share a formula, but reading `nms` here
  ## would make that a coincidence rather than a fact of this fit
  nml <- statmod_design(fl@spec)$mu$coef_names
  b1 <- fl@coefficients$mu[[match("s(z).z1", nml)]]
  s1 <- sqrt(stats::vcov(fl, readable = FALSE)["mu:s(z).z1", "mu:s(z).z1"])
  c1 <- confint(fl, "mu:s(z).z1", test = "lr", readable = FALSE)
  expect_equal(c(c1$lower[[1L]], c1$upper[[1L]]),
               c(b1 - zq * s1, b1 + zq * s1), tolerance = 1e-5,
               ignore_attr = TRUE)
  expect_lt(abs((c1$upper[[1L]] - b1) - (b1 - c1$lower[[1L]])), 1e-6)
})


test_that("with no penalty the penalized objective IS the log-likelihood", {
  ## the control that keeps every reference check above valid: the change of
  ## objective must be invisible wherever there is nothing to be penalized
  d <- pois_data()
  fit <- pois_fit(d)
  expect_equal(fit@objective, -as.numeric(logLik(fit)))
  for (b in c(0, 0.4, 1.2)) {
    r <- statmod_restrict(fit, "mu", "x", b)
    expect_equal(2 * (r$objective - fit@objective),
                 2 * (as.numeric(logLik(fit)) - r$loglik), tolerance = 1e-8)
  }
})


test_that("confint(test =) inverts the test it names", {
  d <- pois_data()
  fit <- pois_fit(d)
  for (m in c("lr", "score", "gradient")) {
    a <- confint(fit, "mu:x", test = m, readable = FALSE)
    b <- statmod_invert(fit, "mu", "x", 0.95, m)
    expect_equal(c(a$lower[[1L]], a$upper[[1L]]), as.numeric(b),
                 ignore_attr = TRUE)
    ## the estimate and the standard error are the fit's own whatever the
    ## test is, and only the two limits move
    expect_equal(a$se[[1L]],
                 confint(fit, "mu:x", readable = FALSE)$se[[1L]])
  }
  ## the default is unchanged
  expect_equal(confint(fit, "mu:x", test = "wald", readable = FALSE),
               confint(fit, "mu:x", readable = FALSE))

  ## a test is about ONE coefficient, so the readable quantities -- which are
  ## functions of several at once -- cannot be asked for
  expect_error(confint(fit, "mu:x", test = "lr"), "readable = FALSE")

  ## a row under a KINKED penalty keeps NA rather than falling back on the
  ## Wald limits, and a request for nothing but such rows is refused
  fl <- lasso_fit()
  nl <- statmod_design(fl@spec)$mu$coef_names
  tb <- confint(fl, test = "lr", readable = FALSE)
  kink <- rownames(tb) != "mu:(Intercept)"
  expect_true(all(is.na(tb$lower[kink])))
  expect_true(is.finite(tb$lower[!kink]))
  expect_error(confint(fl, paste0("mu:", nl[[2L]]), test = "lr",
                       readable = FALSE), "nothing to invert")

  ## while a smooth's own coordinates are inverted like any other
  cs <- confint(smooth_fit(), "mu:s(z).lin", test = "lr",
                readable = FALSE)
  expect_true(all(is.finite(c(cs$lower[[1L]], cs$upper[[1L]]))))
})


test_that("summary(test =) reports the signed root of the test it names", {
  d <- pois_data()
  fit <- pois_fit(d)
  sw <- summary(fit)@tables$mu[[1L]]$table
  ## the default is Wald's z, the estimate over its standard error
  expect_equal(sw$statistic, sw$estimate / sw$se)

  for (t in c("lr", "score", "gradient")) {
    tb <- summary(fit, test = t)@tables$mu[[1L]]$table
    for (i in seq_len(nrow(tb))) {
      s <- statmod_stat_at(fit, "mu", tb$name[[i]], 0, t)
      ## the signed root of the chi-squared value, with the sign of the
      ## estimate, and the p-value of the value itself
      expect_equal(tb$statistic[[i]],
                   sign(tb$estimate[[i]]) * sqrt(s$statistic))
      expect_equal(tb$p_value[[i]], s$p.value)
    }
    ## the intervals stay Wald's: only the statistic and the p-value move
    expect_equal(tb$lower, sw$lower)
    expect_equal(tb$upper, sw$upper)
  }

  ## the summary says which test produced the column
  expect_identical(summary(fit)@test, "wald")
  expect_identical(summary(fit, test = "score")@test, "score")
  out <- utils::capture.output(print(summary(fit, test = "lr")))
  expect_true(any(grepl("lr test, signed root", out, fixed = TRUE)))
  expect_false(any(grepl("lr test", utils::capture.output(
    print(summary(fit))), fixed = TRUE)))
})


test_that("a smooth's linear column is tested and a kinked row is not", {
  s <- summary(smooth_fit(), test = "lr")
  ## the parametric block, and the smooth's own coefficient rows with it --
  ## the linear column is what this test exists for, being unshrunk
  par_tb <- s@tables$mu[[1L]]$table
  expect_true(all(is.finite(par_tb$statistic)))
  sm_tb <- s@tables$mu[[2L]]$table
  cf <- sm_tb$role == "coefficient"
  expect_true(any(cf))
  expect_true(all(is.finite(sm_tb$statistic[cf])))
  ## a hyperparameter row still carries none: the null a test would report
  ## on is that it is zero, which is the edge of its range
  expect_true(all(is.na(sm_tb$statistic[!cf])))
  expect_true(any(grepl("likelihood-ratio", s@notes)))
  expect_false(any(grepl("no statistic", s@notes)))

  ## a KINKED block is where the column has holes, and the note says so
  sl <- summary(lasso_fit(), test = "lr")
  tb <- do.call(rbind, lapply(sl@tables$mu, function(b) b$table))
  cfl <- tb$role == "coefficient"
  expect_true(all(is.na(tb$statistic[cfl & tb$name != "(Intercept)"])))
  expect_true(any(grepl("no statistic", sl@notes)))
})


test_that("the old name of confint's argument is refused rather than ignored", {
  fit <- pois_fit(pois_data())
  ## `method` was this argument's name in 0.103.0. It reaches the DOTS here
  ## rather than any formal -- no formal of confint() begins with it -- and
  ## from there it would go to vcov(), which has no such argument either, so
  ## without the guard the request would be answered with the Wald limits in
  ## silence. That is the failure this exists for, not tidiness.
  expect_error(confint(fit, "mu:x", method = "lr", readable = FALSE),
               "'method' is now 'test'", fixed = TRUE)
  ## the same request under the name summary() already used
  a <- confint(fit, "mu:x", test = "lr", readable = FALSE)
  expect_true(all(is.finite(c(a$lower[[1L]], a$upper[[1L]]))))
})


test_that("statmod_test() carries the null a summary cannot be asked about", {
  d <- pois_data()
  fit <- pois_fit(d)
  bx <- fit@coefficients$mu[[2L]]

  for (t in c("wald", "lr", "score", "gradient")) {
    tt <- statmod_test(fit, "mu", "x", 0.5, t)
    expect_true(S7::S7_inherits(tt, StatmodTest))
    s <- statmod_stat_at(fit, "mu", "x", 0.5, t)
    expect_equal(tt@statistic, s$statistic)
    expect_equal(tt@p.value, s$p.value)
    expect_identical(tt@df, 1L)
    expect_identical(tt@test, t)
    ## what the object says about itself: the estimate, the null it was put
    ## against, and which coefficient of which equation
    expect_equal(tt@estimate, bx)
    expect_identical(tt@null_value, 0.5)
    expect_identical(tt@parameter, "mu")
    expect_identical(tt@coefficient, "x")
  }

  ## tested against its own estimate every statistic vanishes
  for (t in c("lr", "score", "gradient")) {
    tt <- statmod_test(fit, "mu", "x", bx, t)
    expect_lt(abs(tt@statistic), 1e-6)
    expect_gt(tt@p.value, 0.99)
  }

  ## and against zero it is what the summary reports for the same row, which
  ## is where the two surfaces have to agree
  z <- statmod_test(fit, "mu", "x", 0, "lr")
  st <- summary(fit, test = "lr")@tables$mu[[1L]]$table
  j <- match("x", st$name)
  expect_equal(st$statistic[[j]]^2, z@statistic)
  expect_equal(st$p_value[[j]], z@p.value)

  out <- utils::capture.output(print(statmod_test(fit, "mu", "x", 1, "lr")))
  expect_true(any(grepl("Likelihood-ratio test", out, fixed = TRUE)))
  expect_true(any(grepl("mu:x", out, fixed = TRUE)))
  expect_true(any(grepl("is not equal to 1", out, fixed = TRUE)))
})


test_that("what a test says about the refit is the mode error, not the flag", {
  ## The two answer different questions. The flag says whether a stopping
  ## rule fired; whether the point is usable is a matter of distance, and
  ## `restricted_mode_error()` measures it in log-likelihood units. Measured,
  ## they come apart BACKWARDS: over nineteen held values of this fit's slope
  ## every refit is at its mode -- the worst at 1.04e-11 against a limit of
  ## 1e-03 -- while four report the flag as FALSE, and those four are the four
  ## whose mode is located BEST. Where a refit lands on its mode in one step
  ## the objective does not move and the stall guard fires.
  ##
  ## WHICH points the flag rejects is last-bit arithmetic and is not asserted;
  ## that every one of them is at its mode is asserted, with six orders of
  ## room, and so is the rule the print method reads.
  d <- pois_data()
  fit <- pois_fit(d)
  vals <- seq(0.4, 1.2, by = 0.1)
  me <- vapply(vals, function(b) {
    statmod_restrict(fit, "mu", "x", b)$mode_error()
  }, numeric(1))
  expect_true(all(is.finite(me)))
  expect_true(all(me < mode_error_limit()))
  expect_lt(max(me), 1e-6)

  ## the mode error is read only when it is asked for, so the loops that call
  ## the helper in quantity -- the summary and the interval -- pay no Hessian
  expect_true(is.na(statmod_stat_at(fit, "mu", "x", 0.5, "lr")$mode_error))
  asked <- statmod_stat_at(fit, "mu", "x", 0.5, "lr", mode_error = TRUE)
  expect_true(is.finite(asked$mode_error))
  expect_equal(asked$statistic,
               statmod_stat_at(fit, "mu", "x", 0.5, "lr")$statistic)
  ## and statmod_test(), which a reader reads one row of, does ask
  expect_true(is.finite(statmod_test(fit, "mu", "x", 0.5, "lr")@mode_error))
  expect_true(is.na(statmod_test(fit, "mu", "x", 0.5, "wald")@mode_error))

  ## the print rule, pinned on objects built by hand so that neither branch
  ## depends on which points an optimizer happens to flag
  warned <- function(...) {
    any(grepl("stopped above its own mode",
              utils::capture.output(print(StatmodTest(
                test = "lr", statistic = 1, df = 1L, p.value = 0.3,
                estimate = 1, null_value = 0, parameter = "mu",
                coefficient = "x", ...))), fixed = TRUE))
  }
  ## a located mode says nothing, WHATEVER the flag says -- this is the case
  ## the change is about, and reading the flag would warn here
  expect_false(warned(mode_error = 1e-20, converged = FALSE))
  expect_false(warned(mode_error = 1e-20, converged = TRUE))
  ## a mode that is NOT located says so, whatever the flag says
  expect_true(warned(mode_error = 1, converged = TRUE))
  expect_true(warned(mode_error = NA_real_, converged = TRUE))
  ## and Wald, which refits nothing, has no mode to be above
  expect_false(warned(mode_error = NA_real_, converged = NA))
})


test_that("a held value of another length is refused where it is passed", {
  fit <- pois_fit(pois_data())
  ## it used to reach the hold as a named vector of that length and fail
  ## several frames down, inside the objective, on a names assignment -- an
  ## error naming neither the argument nor the mistake
  expect_error(statmod_restrict(fit, "mu", "x", NULL), "'value'", fixed = TRUE)
  expect_error(statmod_restrict(fit, "mu", "x", c(1, 2)), "'value'",
               fixed = TRUE)
  expect_error(statmod_restrict(fit, "mu", "x", NA_real_), "'value'",
               fixed = TRUE)
  expect_error(statmod_test(fit, "mu", "x", numeric(0), "lr"), "'value'",
               fixed = TRUE)
})
