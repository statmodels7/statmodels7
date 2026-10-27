g1 <- distributions7::gaussian1_distrib()

# a gaussian random-intercept fit and the data it came from
ri_data <- function(G = 12, m = 4, seed = 21) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(G), each = m)), x = rnorm(G * m))
  d$y <- 1 + 0.5 * d$x + 1.1 * rnorm(G)[d$g] + rnorm(G * m)
  d
}

# the draws of a cluster bootstrap, as the engine makes them from set.seed(s)
cluster_picks <- function(s, B, G) {
  set.seed(s)
  lapply(seq_len(B), function(i) sample.int(G, G, TRUE))
}

# the replicate data: the groups drawn, each draw a new group
stack_picks <- function(d, pick) {
  idx <- split(seq_len(nrow(d)), match(d$g, unique(d$g)))
  rows <- unlist(idx[pick], use.names = FALSE)
  inst <- rep(seq_along(pick), lengths(idx[pick]))
  db <- d[rows, , drop = FALSE]
  db$newid <- factor(paste(d$g[rows], inst, sep = "#"))
  db
}

# the fixed coefficients and log sigma of a gaussian random intercept whose
# effects have standard deviation tau: sigma maximizes the restricted
# likelihood of V = sigma^2 I + tau^2 Z Z' and the coefficients are the
# generalized least squares ones at that V
ri_closed <- function(db, tau) {
  X <- cbind(1, db$x)
  Z <- stats::model.matrix(~ 0 + newid, db)
  ZZ <- tcrossprod(Z)
  y <- db$y
  gls <- function(ls) {
    V <- exp(2 * ls) * diag(nrow(X)) + tau^2 * ZZ
    Vi <- solve(V)
    A <- crossprod(X, Vi %*% X)
    beta <- solve(A, crossprod(X, Vi %*% y))
    r <- y - X %*% beta
    ll <- -0.5 * (as.numeric(determinant(V)$modulus) +
                    as.numeric(determinant(A)$modulus) +
                    sum(r * (Vi %*% r)))
    list(beta = as.numeric(beta), ll = ll)
  }
  ls <- stats::optimize(function(z) -gls(z)$ll, c(-4, 2), tol = 1e-13)$minimum
  c(gls(ls)$beta, ls)
}

test_that("the regrouped cluster bootstrap is the cluster bootstrap with new groups", {
  d <- ri_data()
  fit <- statmod(y ~ x + random(~ 1 | g), g1, d)
  tau <- fit@hyper$mu[[1L]][["sigma"]]
  B <- 10
  # the closed form is the model: at the data it reproduces the fit
  d$newid <- d$g
  expect_equal(ri_closed(d, tau), unname(c(fit@coefficients$mu[1:2],
                                           fit@coefficients$sigma)),
               tolerance = 1e-6)
  set.seed(31)
  fit <- add_bootstrap(fit, n = B, resample = cluster_resample(~ g))
  e <- statmod_bootstrap(fit, resample = cluster_resample(~ g))
  expect_identical(e@random, "regrouped")
  expect_equal(sum(e@status), 0)
  picks <- cluster_picks(31, B, 12)
  ref <- t(vapply(picks, function(p) ri_closed(stack_picks(d, p), tau),
                  numeric(3)))
  keep <- c(1:2, ncol(coef(e)))
  expect_equal(unname(coef(e)[, keep]), ref, tolerance = 1e-6)
  # the variance of the fixed coefficients
  V <- vcov(fit, vcov_type = bootstrap_vcov(resample = cluster_resample(~ g)))
  expect_equal(unname(V[c(1, 2), c(1, 2)]), unname(stats::cov(ref[, 1:2])),
               tolerance = 1e-6)
})

test_that("a replicate refitted in full is the fit of the stacked data", {
  d <- ri_data(G = 10, m = 4, seed = 4)
  fit <- statmod(y ~ x + random(~ 1 | g), g1, d)
  B <- 5
  set.seed(8)
  fit <- add_bootstrap(fit, n = B, resample = cluster_resample(~ g),
                       refit = "full")
  e <- statmod_bootstrap(fit, resample = cluster_resample(~ g),
                         refit = "full")
  picks <- cluster_picks(8, B, 10)
  brute <- lapply(picks, function(p) {
    db <- stack_picks(d, p)
    statmod(y ~ x + random(~ 1 | newid), g1, db)
  })
  bf <- t(vapply(brute, function(f)
    c(f@coefficients$mu[1:2], f@coefficients$sigma), numeric(3)))
  keep <- c(1:2, ncol(coef(e)))
  expect_equal(max(abs(unname(coef(e)[, keep]) - bf)), 0, tolerance = 1e-4)
  bh <- vapply(brute, function(f) f@hyper$mu[[1L]][["sigma"]], numeric(1))
  expect_equal(unname(e@hyper[, 1L]), bh, tolerance = 1e-3)
  # the variance components are stored, finite and not constant
  expect_true(all(is.finite(e@hyper)))
  expect_gt(stats::sd(e@hyper[, 1L]), 0)
})

test_that("the individual effects have no variance under the cluster scheme", {
  d <- ri_data(G = 8, m = 4, seed = 5)
  fit <- statmod(y ~ x + random(~ 1 | g), g1, d)
  set.seed(1)
  fit <- add_bootstrap(fit, n = 8, resample = cluster_resample(~ g))
  V <- vcov(fit, vcov_type = bootstrap_vcov(resample = cluster_resample(~ g)))
  lab <- coef_labels(fit@spec, statmod_design(fit@spec))
  rnd <- which(lab$term == "random(~1 | g)")
  expect_length(rnd, 8)
  expect_true(all(is.na(V[rnd, ])) && all(is.na(V[, rnd])))
  expect_false(anyNA(V[-rnd, -rnd]))
  expect_true(all(is.finite(diag(V)[-rnd])))
  # the standard errors of the summary leave the effects out
  s <- suppressWarnings(summary(fit, vcov_type = bootstrap_vcov(
    resample = cluster_resample(~ g))))
  expect_true(S7::S7_inherits(s, StatmodSummary))
})

test_that("500 and then 300 are the 800 of one run (small numbers)", {
  d <- ri_data(G = 8, m = 4, seed = 6)
  mk <- function() statmod(y ~ x + random(~ 1 | g), g1, d)
  rs <- cluster_resample(~ g)
  set.seed(9)
  a <- add_bootstrap(mk(), n = 14, resample = rs)
  set.seed(9)
  b <- add_bootstrap(mk(), n = 6, resample = rs)
  b <- add_bootstrap(b, n = 14, resample = rs)
  ca <- coef(statmod_bootstrap(a, resample = rs))
  cb <- coef(statmod_bootstrap(b, resample = rs))
  expect_identical(ca, cb)
  set.seed(9)
  c6 <- add_bootstrap(mk(), n = 6, resample = rs)
  expect_identical(coef(statmod_bootstrap(c6, resample = rs)),
                   ca[1:6, , drop = FALSE])
})

test_that("a smooth beside a random effect keeps the knots of the fit", {
  set.seed(12)
  G <- 10
  d <- data.frame(g = factor(rep(seq_len(G), each = 6)), x = runif(G * 6))
  d$y <- sin(2 * pi * d$x) + rnorm(G)[d$g] + rnorm(G * 6, sd = 0.3)
  fit <- statmod(y ~ s(x, basis7::bspline_smooth(k = 6)) + random(~ 1 | g),
                 g1, d)
  rs <- cluster_resample(~ g)
  plan <- bootstrap_plan(fit, rs, fit@spec, statmod_design(fit@spec),
                         random_modes(fit@spec, "zero"), TRUE)
  build <- attr(plan, "build")
  set.seed(2)
  cb <- build(plan())
  sm <- fit@spec@terms$mu[[2L]]
  db <- cb$spec@newdata
  # the smooth block of the replicate is the fit's basis at the rows drawn
  des <- statmod_design(cb$spec)
  Xs <- as.matrix(des$mu$X[, des$mu$blocks[[2L]]])
  expect_equal(Xs, as.matrix(modelterms7::term_predict(sm, newdata = db)),
               ignore_attr = TRUE)
  # and a basis built again on a few groups would have other knots
  cb2 <- build(list(pick = c(1, 1, 1, 2)))
  Xnew <- as.matrix(modelterms7::term_matrix(modelterms7::term_build(
    sm, cb2$spec@newdata)))
  des2 <- statmod_design(cb2$spec)
  Xs2 <- as.matrix(des2$mu$X[, des2$mu$blocks[[2L]]])
  expect_gt(max(abs(Xs2 - Xnew)), 1e-2)
  expect_equal(Xs2, as.matrix(modelterms7::term_predict(
    sm, newdata = cb2$spec@newdata)), ignore_attr = TRUE)
  set.seed(3)
  fit <- add_bootstrap(fit, n = 4, resample = rs)
  e <- statmod_bootstrap(fit, resample = rs)
  lab <- coef_labels(fit@spec, statmod_design(fit@spec))
  expect_true(all(is.finite(coef(e)[, lab$term != "random(~1 | g)"])))
  expect_true(all(is.na(coef(e)[, lab$term == "random(~1 | g)"])))
})

test_that("a group nested in the cluster is relabeled, a crossed one is not", {
  set.seed(13)
  d <- data.frame(sch = rep(1:4, each = 15), g = rep(1:12, each = 5),
                  x = rnorm(60))
  d$cr <- rep(1:5, 12)
  d$y <- 1 + d$x + rnorm(12)[d$g] + rnorm(5)[d$cr] + rnorm(60)
  fit <- statmod(y ~ x + random(~ 1 | g) + random(~ 1 | cr, label = "rc"), g1, d)
  rs <- cluster_resample(~ sch)
  plan <- bootstrap_plan(fit, rs, fit@spec, statmod_design(fit@spec),
                         random_modes(fit@spec, "zero"), TRUE)
  build <- attr(plan, "build")
  cb <- build(list(pick = c(2, 2, 3, 1)))
  db <- cb$spec@newdata
  expect_equal(nlevels(factor(db$.boot_group1)), 3 * 2 + 3 + 3)
  expect_equal(nlevels(factor(db$.boot_group2)), 5)
  # the coordinates in common are the intercept, the slope and sigma
  expect_equal(sum(!is.na(cb$map)), 3)
})

test_that("a fit that cannot be regrouped is rejected", {
  d <- ri_data(G = 8, m = 4, seed = 7)
  d$u <- d$x
  fl <- statmod(y ~ x + random(~ 1 | u | g), g1, d)
  expect_error(add_bootstrap(fl, n = 3, resample = cluster_resample(~ g)),
               "label")
  fit <- statmod(y ~ x + random(~ 1 | g), g1, d)
  expect_error(add_bootstrap(fit, n = 3, resample = case_resample()),
               "breaks the groups")
  # a variable of the model that is not a column of the data
  d2 <- d[sample(nrow(d)), ]
  zz <- rnorm(nrow(d2))
  fz <- statmod(y ~ zz + random(~ 1 | g), g1, d2)
  expect_error(add_bootstrap(fz, n = 3, resample = cluster_resample(~ g)),
               "not a column of the data")
})

test_that("a replicate that leaves a level at weight zero is counted failed", {
  stream <- function(s, n, B, rare) {
    set.seed(s)
    bad <- logical(0)
    while (sum(!bad) < B) {
      w <- tabulate(sample.int(n, n, TRUE), n)
      bad <- c(bad, sum(w[rare]) == 0)
    }
    bad
  }
  # a logistic fit: the refit returns a value for the rare level, and the
  # engine counts the replicate as unidentified
  set.seed(3)
  n <- 200
  d <- data.frame(f = factor(c(rep("rare", 3), sample(c("b", "c"), n - 3,
                                                        TRUE))))
  d$y <- c(1, 1, 0, rbinom(n - 3, 1, 0.5))
  fb <- statmod(y ~ f, distributions7::bernoulli_distrib(), d)
  B <- 40
  set.seed(15)
  fb <- add_bootstrap(fb, n = B, resample = case_resample())
  e <- statmod_bootstrap(fb, resample = case_resample())
  bad <- stream(15, n, B, 1:3)
  expect_gt(sum(bad), 0)
  expect_identical(e@status, ifelse(bad, 3L, 0L))
  expect_equal(statmod_bootstrap(fb)$unidentified, sum(bad))
  expect_equal(nrow(coef(e)), B)
  expect_true(all(is.finite(coef(e))))
  # a gaussian fit with a dispersion on the marginal criterion: the refit
  # signals an error, and the replicate is failed and replaced as well
  set.seed(14)
  n <- 100
  d <- data.frame(f = factor(c(rep("rare", 3), sample(c("b", "c"), n - 3,
                                                        TRUE))))
  d$y <- 1 + (d$f == "c") + 3 * (d$f == "rare") + rnorm(n)
  fit <- statmod(y ~ f, g1, d)
  set.seed(15)
  fit <- add_bootstrap(fit, n = 30, resample = case_resample())
  e <- statmod_bootstrap(fit, resample = case_resample())
  bad <- stream(15, n, 30, 1:3)
  expect_gt(sum(bad), 0)
  expect_identical(e@status > 0L, bad)
  expect_true(all(is.finite(coef(e))))
})
