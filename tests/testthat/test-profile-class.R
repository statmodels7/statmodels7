# The object that statmod_profile() returns, its intervals and its graphs.
# The reference of the intervals is the inversion of the likelihood ratio
# test in confint.StatmodFit(), an independent route to the same limits.

lin_fit <- function() {
  set.seed(1)
  dd <- data.frame(x = runif(80))
  dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
  statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
}

smooth_fit <- function() {
  set.seed(3)
  d <- data.frame(z = runif(200), x = rnorm(200))
  d$y <- rpois(200, exp(1 + sin(2 * pi * d$z) + 0.4 * d$x))
  statmod(y ~ s(z) + x, distributions7::poisson_distrib(), d)
}

test_that("the object records the labels, their kinds and the profile column", {
  fit <- smooth_fit()
  p <- statmod_profile(fit, c("mu/s(z)/lambda", "mu:s(z).z3", "mu:x"),
                       values = profile_range(width = 1, n = 3), joint = FALSE)
  expect_true(S7::S7_inherits(p, StatmodProfile))
  expect_identical(p@profiled_labels, c("mu/s(z)/lambda", "mu:s(z).z3", "mu:x"))
  q <- p@profiled_quantities
  expect_true(S7::S7_inherits(q[["mu/s(z)/lambda"]], ProfiledHyperparameter))
  expect_true(S7::S7_inherits(q[["mu:s(z).z3"]], ProfiledPenalizedCoefficient))
  expect_true(S7::S7_inherits(q[["mu:x"]], ProfiledCoefficient))
  expect_identical(unname(p@profile_column),
                   c("outer_criterion", "profile_pen_loglik", "profile_pen_loglik"))
  expect_true(all(p@is_likelihood_ratio))
  expect_false(p@joint)
  expect_identical(names(p@points)[1:2], c("profiled_label", "value"))
  expect_identical(nrow(p@points), 9L)
  expect_identical(p@unrestricted_point$profiled_label, p@profiled_labels)
  expect_identical(as.data.frame(p), p@points)
  expect_output(print(p), "penalized coefficient")
  expect_identical(p@hyper_handling, "held")
  expect_identical(p@outer_method@kind, "reml")
})

test_that("a separate profile equals the profile of each label by itself", {
  fit <- lin_fit()
  v <- profile_range(width = 2, n = 5)
  p <- statmod_profile(fit, c("mu:x", "mu:(Intercept)"), values = v, joint = FALSE)
  for (lab in p@profiled_labels) {
    one <- statmod_profile(fit, lab, values = v)
    rows <- p@points$profiled_label == lab
    expect_equal(p@points$value[rows], one@points[[lab]])
    expect_equal(p@points$profile_pen_loglik[rows], one@points$profile_pen_loglik,
                 tolerance = 1e-10)
  }
  expect_error(statmod_profile(fit, c("mu:x", "mu:(Intercept)"),
                               values = list("mu:z" = 1:3), joint = FALSE),
               "not in 'which'")
})

test_that("confint() inverts the likelihood ratio like confint() of the fit", {
  fit <- lin_fit()
  p <- statmod_profile(fit, "mu:x", values = profile_range(width = 3, n = 13))
  ci <- confint(p)
  ref <- confint(fit, parm = "mu:x", test = lr_test(), readable = FALSE)
  expect_equal(ci$lower, ref$lower, tolerance = 1e-6)
  expect_equal(ci$upper, ref$upper, tolerance = 1e-6)
  expect_identical(rownames(ci), "mu:x")
  ci9 <- confint(p, level = 0.9)
  ref9 <- confint(fit, parm = "mu:x", level = 0.9, test = lr_test(), readable = FALSE)
  expect_equal(c(ci9$lower, ci9$upper), c(ref9$lower, ref9$upper), tolerance = 1e-6)
})

test_that("simultaneous intervals use the chi-squared quantile with k degrees of freedom", {
  fit <- lin_fit()
  p <- statmod_profile(fit, c("mu:x", "mu:(Intercept)"),
                       values = profile_range(width = 4, n = 17), joint = FALSE)
  s <- confint(p, type = "simultaneous")
  lev <- stats::pchisq(stats::qchisq(0.95, 2), 1)
  ref <- confint(fit, parm = c("mu:x", "mu:(Intercept)"), level = lev,
                 test = lr_test(), readable = FALSE)
  expect_equal(s[rownames(ref), "lower"], ref$lower, tolerance = 1e-5)
  expect_equal(s[rownames(ref), "upper"], ref$upper, tolerance = 1e-5)
  m <- confint(p)
  expect_true(all(s$lower < m$lower & s$upper > m$upper))
})

test_that("a joint profile gives intervals from the maximum over the other labels", {
  fit <- lin_fit()
  j <- statmod_profile(fit, c("mu:x", "mu:(Intercept)"),
                       values = profile_range(width = 3, n = 9))
  ci <- confint(j)
  ref <- confint(fit, parm = c("mu:x", "mu:(Intercept)"), test = lr_test(),
                 readable = FALSE)
  # the grid of the other label resolves the maximum to a few per cent of the
  # width of the interval
  w <- ref$upper - ref$lower
  expect_lt(max(abs(ci[rownames(ref), "lower"] - ref$lower) / w), 0.05)
  expect_lt(max(abs(ci[rownames(ref), "upper"] - ref$upper) / w), 0.05)
})

test_that("a limit outside the values is missing, with a warning", {
  fit <- lin_fit()
  p <- statmod_profile(fit, "mu:x", values = profile_range(width = 0.5, n = 5))
  expect_warning(ci <- confint(p), "Widen the values")
  expect_true(is.na(ci$lower) && is.na(ci$upper))
  expect_error(confint(p, parm = "mu:z"), "does not hold")
  expect_error(confint(p, level = 1.2), "between 0 and 1")
})

test_that("confint() rejects a profile whose drop is not a likelihood ratio", {
  fit <- smooth_fit()
  p <- statmod_profile(fit, c("mu/s(z)/lambda", "mu:x"),
                       values = list("mu/s(z)/lambda" = c(1, 1.3),
                                     "mu:x" = c(0.35, 0.4)))
  expect_false(any(p@is_likelihood_ratio))
  expect_identical(unname(p@profile_column), rep("outer_criterion", 2))
  expect_error(confint(p), "does not give a likelihood ratio")
  expect_output(print(p), "no likelihood ratio")
})

test_that("each kind of label is drawn by its own graph", {
  fit <- smooth_fit()
  p <- statmod_profile(fit, c("mu/s(z)/lambda", "mu:x"),
                       values = profile_range(width = 1, n = 3), joint = FALSE)
  g <- autoplot(p)
  expect_s3_class(g, "patchwork")
  one <- statmod_profile(fit, "mu:x", values = profile_range(width = 1, n = 3))
  g1 <- autoplot(one)
  expect_s3_class(g1, "ggplot")
  expect_identical(g1$labels$y, "penalized log-likelihood")
  expect_identical(g1$labels$title, "mu: x")
  gh <- autoplot(statmod_profile(fit, "mu/s(z)/lambda",
                                 values = profile_range(width = 1, n = 3)))
  expect_identical(gh$labels$y, "REML criterion")
  expect_identical(gh$labels$x, "log lambda")
  set.seed(5)
  w <- runif(120, 0, 10)
  db <- data.frame(w = w, y = 1 + 0.8 * pmax(w - 4, 0) + rnorm(120))
  fb <- statmod(y ~ seg(w), distributions7::gaussian1_distrib(), db)
  pb <- statmod_profile(fb, "mu:seg.psi1", values = profile_range(width = 1, n = 3))
  expect_true(S7::S7_inherits(pb@profiled_quantities[[1L]], ProfiledBreakpoint))
  expect_length(pb@profiled_quantities[[1L]]@covariate_values, 120L)
  gb <- autoplot(pb)
  expect_identical(gb$labels$y, "log-likelihood")
  expect_true(any(vapply(gb$layers, function(l) inherits(l$geom, "GeomRug"), TRUE)))
})

test_that("a joint profile is drawn as contours, one panel for each pair", {
  fit <- lin_fit()
  v <- profile_range(width = 2, n = 3)
  j2 <- statmod_profile(fit, c("mu:x", "mu:(Intercept)"), values = v)
  g2 <- autoplot(j2)
  expect_s3_class(g2, "ggplot")
  expect_true(any(vapply(g2$layers, function(l) inherits(l$geom, "GeomContour"), TRUE)))
  j3 <- statmod_profile(fit, c("mu:x", "mu:(Intercept)", "sigma:(Intercept)"),
                        values = v)
  g3 <- autoplot(j3)
  expect_s3_class(g3, "patchwork")
  expect_length(g3$patches$plots, 2L)
})

bp_data <- function() {
  set.seed(5)
  n <- 300
  w <- round(runif(n, 0, 10), 3)
  d <- data.frame(w = w)
  d$y <- 1 + 0.2 * w + 1.2 * pmax(w - 4, 0) + rnorm(n, sd = 0.8)
  d$yj <- 1 + 0.2 * w + 1.5 * (w > 4) + rnorm(n, sd = 0.8)
  d$ys <- 1 + 0.2 * w + 1.2 * (w > 4) + 0.8 * pmax(w - 4, 0) + rnorm(n, sd = 0.8)
  d
}

test_that("the position of a sharp jump() or jseg() is held by moving the term", {
  d <- bp_data()
  fam <- distributions7::gaussian1_distrib()
  fj <- statmod(yj ~ w + jump(w), fam, d, outer_criterion = NULL)
  expect_error(statmod_profile(fj, "mu:jump.g1"), "Hold the position itself, 'mu:jump.psi1'")
  p <- statmod_profile(fj, "mu:jump.psi1", values = profile_range(n = 7))
  q <- p@profiled_quantities[[1L]]
  expect_true(S7::S7_inherits(q, ProfiledBreakpoint))
  expect_true(q@piecewise_constant)
  v <- p@points[["mu:jump.psi1"]]
  u <- sort(unique(d$w))
  # the default values are midpoints between consecutive observations
  expect_true(all(abs(v - (u[findInterval(v, u)] + u[findInterval(v, u) + 1L]) / 2) < 1e-12))
  ref <- vapply(v, function(q) as.numeric(logLik(stats::lm(yj ~ w + I(w > q), d))), 1)
  expect_equal(p@points$profile_loglik, ref, tolerance = 1e-9)
  ci <- confint(p)
  expect_true(all(c(ci$lower, ci$upper) %in% u))
  expect_s3_class(autoplot(p), "ggplot")
  fs <- statmod(ys ~ jseg(w), fam, d, outer_criterion = NULL)
  ps <- statmod_profile(fs, "mu:jseg.psi1", values = profile_range(n = 5))
  vs <- ps@points[["mu:jseg.psi1"]]
  refs <- vapply(vs, function(q)
    as.numeric(logLik(stats::lm(ys ~ w + I(w > q) + pmax(w - q, 0), d))), 1)
  expect_equal(ps@points$profile_loglik, refs, tolerance = 1e-9)
})

test_that("a seg() position held next to an observation refits from its neighbour", {
  d <- bp_data()
  f <- statmod(y ~ seg(w), distributions7::gaussian1_distrib(), d, outer_criterion = ml())
  est <- coef(f)$mu[["seg.psi1"]]
  u <- sort(unique(d$w))
  i <- findInterval(est, u)
  v <- sort(unique(c(seq(u[i - 1L], u[i + 2L], length.out = 25), u[(i - 1L):(i + 2L)])))
  p <- statmod_profile(f, "mu:seg.psi1", values = list(v))
  expect_false(anyNA(p@points$profile_loglik))
  expect_true(all(p@points$converged))
  ref <- vapply(v, function(q)
    as.numeric(logLik(stats::lm(y ~ w + pmax(w - q, 0), d), REML = FALSE)), 1)
  expect_equal(p@points$profile_loglik, ref, tolerance = 1e-9)
  # the default values include the observations inside the window
  pd <- statmod_profile(f, "mu:seg.psi1", values = profile_range(width = 2, n = 5))
  expect_gt(sum(pd@points[["mu:seg.psi1"]] %in% u), 0L)
})

test_that("under reml() each point reaches the criterion's own optimum of the dispersion", {
  d <- bp_data()
  f <- statmod(y ~ seg(w), distributions7::gaussian1_distrib(), d)
  ctx <- profile_context(f, "mu:seg.psi1", "held")
  # the REML criterion with the position and the dispersion both held
  crit_at <- function(psi, s) {
    sp <- profile_spec(ctx, psi)
    hc <- sp@held_coef
    hc$sigma <- c("(Intercept)" = s)
    sp <- S7::set_props(sp, held_coef = hc)
    dd <- statmod_design(sp)
    b <- ctx$beta_fit
    b[ctx$at] <- psi
    b[length(b)] <- s
    rr <- statmod_alternate(sp, dd, statmod_blocks(sp, dd), ctx$hyper0, ctx$inner,
                            b, ctx$cfg$expected, ctx$cfg$emethod, ctx$cfg$maxit,
                            ctx$cfg$tol, verbosity(0))
    statmod_marginal(sp, dd, rr$obj$split(rr$par), ctx$hyper0, ctx$outer,
                     basis = integrated_basis(sp, dd, ctx$outer@kind))$value
  }
  # the point next to the estimate stopped at the dispersion of the fit,
  # 0.0008 REML units below its optimum
  v <- ctx$estimate + c(-0.0025, 0.0005, 0.0035)
  p <- statmod_profile(f, "mu:seg.psi1", values = list(v))
  s0 <- utils::tail(ctx$beta_fit, 1L)
  ref <- vapply(v, function(psi) optimize(function(s) crit_at(psi, s),
                                          s0 + c(-0.04, 0.04), maximum = TRUE,
                                          tol = 1e-8)$objective, 1)
  expect_equal(p@points$outer_criterion, ref, tolerance = 1e-7)
})

test_that("the lambda of a path is profiled on the path's own criterion", {
  set.seed(7)
  n <- 150
  X <- matrix(rnorm(n * 6), n, 6, dimnames = list(NULL, paste0("x", 1:6)))
  dl <- data.frame(X, y = X[, 1] - 0.5 * X[, 2] + rnorm(n))
  fam <- distributions7::gaussian1_distrib()
  fb <- statmod(y ~ lasso(~ x1 + x2 + x3 + x4 + x5 + x6), fam, dl)
  key <- paste0("mu/", names(fb@hyper$mu)[1L], "/lambda")
  h <- fb@history$outer
  # the path's own values away from its top, where the active set is settled
  keep <- h$value < 0.5 * max(h$value)
  p <- statmod_profile(fb, key, values = list(log(h$value[keep])))
  expect_identical(unname(p@profile_column), "path_criterion")
  expect_false(p@is_likelihood_ratio[[1L]])
  got <- p@points$path_criterion[match(log(h$value[keep]), p@points[[key]])]
  expect_equal(got, h$criterion[keep], tolerance = 1e-5)
  expect_equal(p@unrestricted_point$path_criterion, fb@criterion[[1L]], tolerance = 1e-8)
  q <- p@profiled_quantities[[1L]]
  expect_identical(q@criterion_name, "BIC")
  expect_identical(nrow(q@path_record), nrow(h))
  g <- autoplot(p)
  expect_identical(g$labels$y, "BIC")
  set.seed(11)
  fc <- statmod(y ~ lasso(~ x1 + x2 + x3 + x4 + x5 + x6), fam, dl,
                sparse_criterion = cv(nfolds = 5))
  pc <- statmod_profile(fc, key, values = profile_range(n = 3))
  expect_true(all(is.na(pc@points$path_criterion)))
  expect_true(any(is.finite(pc@profiled_quantities[[1L]]@path_record$se)))
  expect_identical(autoplot(pc)$labels$y, "cross-validated criterion")
})

test_that("a profile of an approximate method gives intervals and graphs", {
  data(quine, package = "MASS")
  fit <- statmod(Days ~ Eth + Lrn, distrib = distributions7::poisson_distrib(), data = quine)
  pg <- statmod_profile(fit, which = c("mu:EthN", "mu:LrnSL"), joint = FALSE)
  pt <- statmod_profile(fit, which = c("mu:EthN", "mu:LrnSL"), joint = FALSE,
                        method = profile_taylor())
  expect_true(all(is.na(pt@points$converged)))
  ct <- confint(pt)
  cg <- confint(pg)
  expect_equal(ct$lower, cg$lower, tolerance = 1e-3)
  expect_equal(ct$upper, cg$upper, tolerance = 1e-3)
  expect_s3_class(autoplot(pt), "ggplot")
  pj <- statmod_profile(fit, which = c("mu:EthN", "mu:LrnSL"),
                        values = profile_range(width = 3, n = 9),
                        method = profile_taylor())
  expect_s3_class(autoplot(pj), "ggplot")
  expect_true(all(is.finite(unlist(confint(pj, type = "simultaneous")[, c("lower", "upper")]))))
})
