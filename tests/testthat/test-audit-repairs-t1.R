# Repairs found by the documentation audit of the fit (statmodels7 0.199.6).
# Every block fails on 0.199.6.

g1 <- distributions7::gaussian1_distrib()

test_that("wcrossprod() rejects blocks of different heights on every route", {
  set.seed(1)
  A <- matrix(stats::rnorm(5000 * 10), 5000)
  B <- matrix(stats::rnorm(4000 * 10), 4000)
  w <- stats::runif(5000)
  expect_error(wcrossprod(A, w, B, 1L), "non-conformable")
  expect_error(wcrossprod(A, w, B, 2L), "non-conformable")
})

test_that("iwls() validates its budget and step_halving counts halvings", {
  expect_error(iwls(step_halving = -1), "step_halving")
  expect_error(iwls(step_halving = NA_real_), "step_halving")
  expect_error(iwls(step_halving = c(1, 2)), "step_halving")
  expect_error(iwls(maxit = Inf), "maxit")
  expect_error(iwls(maxit = 2.5), "maxit")
  expect_error(iwls(tol = NA_real_), "tol")
  expect_error(iwls(solver = "chol_crossprod"), "must be a solver object")
  # step_halving = 0 still tries the full step
  set.seed(5)
  dp <- data.frame(x = stats::runif(200))
  dp$y <- stats::rpois(200, exp(1 + 2 * dp$x))
  f <- statmod(y ~ x, distributions7::poisson_distrib(), dp,
               inner_optimizer = iwls(step_halving = 0))
  ref <- stats::glm(y ~ x, stats::poisson, dp)
  expect_true(f@converged)
  expect_equal(as.numeric(stats::logLik(f)), as.numeric(stats::logLik(ref)),
               tolerance = 1e-8)
})

test_that("the final verdict of iwls_fit() reads the free coordinates", {
  obj <- list(fn = function(b) 0.5 * sum((b - c(1, 5))^2),
              gr = function(b) b - c(1, 5))
  pieces <- function(b) list(R = NULL, C = NULL, A = diag(2))
  r <- iwls_fit(obj, c(0, 0), iwls(solver = solve_chol(), maxit = 1), n = 1,
                pieces_at = pieces, frozen = 2L)
  expect_equal(r$par, c(1, 0))
  expect_true(r$converged)
})

test_that("infinite weights and offsets, and an unnamed start, are rejected", {
  set.seed(2)
  d <- data.frame(x = stats::runif(150))
  d$y <- 1 + d$x + stats::rnorm(150, sd = 0.3)
  expect_error(statmod(y ~ x, g1, d, weights = c(Inf, rep(1, 149))),
               "non-negative, finite and complete")
  expect_error(statmod(y ~ x, g1, d, offsets = list(mu = c(Inf, rep(0, 149)))),
               "offset for 'mu' is not finite")
  expect_error(statmod(y ~ x, g1, d, start = list(c(1, 2))),
               "'start' must be a named list")
  expect_error(statmod(y ~ x, g1, d,
                       start = start_random(fn = function(n) rep(NA_real_, n))),
               "missing or non-numeric values for 'mu'")
})

test_that("criterion_tol() takes the tightest part of crit_all()", {
  expect_equal(criterion_tol(optimizers7::crit_all(optimizers7::crit_grad(1e-3),
                                                   optimizers7::crit_grad(1e-9))),
               1e-9)
  expect_equal(criterion_tol(optimizers7::crit_any(optimizers7::crit_grad(1e-3),
                                                   optimizers7::crit_grad(1e-9))),
               1e-3)
})

test_that("a kinked zero whose reading cannot be computed is not settled", {
  set.seed(9)
  n <- 150
  X <- matrix(stats::rnorm(n * 4), n)
  dd <- data.frame(y = 1 + X[, 1] + stats::rnorm(n))
  dd$X <- X
  fit <- statmod(y ~ lasso(X, lambda = 20), g1, dd)
  sp <- fit@spec
  des <- statmod_design(sp)
  obj <- statmod_objective(sp, fit@hyper, des, FALSE, distributions7::expect_bartlett())
  b <- obj$stack(fit@coefficients)
  expect_true(any(b == 0))
  pen <- statmod_penalty_at
  testthat::local_mocked_bindings(
    statmod_penalty_at = function(spec, coef, hyper, design, what, ...) {
      # the reading probes the penalty's gradient just off each kink
      off <- sqrt(.Machine$double.eps) * .Machine$double.eps
      if (identical(what, "gradient") && any(unlist(coef) == off)) {
        stop("unreadable")
      }
      pen(spec, coef, hyper, design, what, ...)
    })
  r <- alternation_readings(sp, des, obj, b, fit@hyper, FALSE, distributions7::expect_bartlett())
  expect_identical(r[["kkt"]], Inf)
  expect_false(alternation_settled(sp, des, obj, b, fit@hyper, FALSE,
                                   distributions7::expect_bartlett()))
})

test_that("rows of weight zero keep the lasso on coordinate descent", {
  set.seed(9)
  n <- 150
  X <- matrix(stats::rnorm(n * 4), n)
  dd <- data.frame(y = 1 + X[, 1] + stats::rnorm(n))
  dd$X <- X
  fa <- statmod(y ~ lasso(X, lambda = 10), g1, dd,
                weights = c(rep(0, 20), rep(1, 130)))
  fb <- statmod(y ~ lasso(X, lambda = 10), g1, dd[21:n, , drop = FALSE])
  expect_equal(unlist(coef(fa)), unlist(coef(fb)), tolerance = 1e-10)
})

test_that("the path warns only at its dense end", {
  for (s in 1:3) {
    set.seed(s)
    n <- 100
    X <- matrix(stats::rnorm(n * 6), n)
    dd <- data.frame(y = stats::rnorm(n))
    dd$X <- X
    expect_no_warning(statmod(y ~ lasso(X), g1, dd))
  }
})

test_that("the held-out rows of cv() carry their offsets and weights", {
  set.seed(11)
  n <- 120
  x <- matrix(stats::rnorm(n * 4), n, 4)
  o <- 3 * sin(seq_len(n) / 10)
  d <- data.frame(y = 1 + 2 * x[, 1] - 1.5 * x[, 2] + o + stats::rnorm(n),
                  o = o)
  d$x <- x
  fo <- rep(1:5, length.out = n)
  fa <- statmod(y ~ lasso(x, n_lambda = 8), g1, d, offsets = list(mu = o),
                sparse_criterion = cv(folds = fo))
  fb <- statmod(y ~ lasso(x, n_lambda = 8) + offset(o), g1, d,
                sparse_criterion = cv(folds = fo))
  expect_equal(fa@criterion, fb@criterion, tolerance = 1e-8)
  # rows of weight zero do not enter the held-out deviance
  sp <- statmod_spec(y ~ lasso(x, n_lambda = 8), g1, d)
  hy <- statmod_hyper_start(sp)
  hy$mu[[1L]][["lambda"]] <- 2
  extra <- d[1:30, , drop = FALSE]
  extra$y <- extra$y + 50
  dA <- rbind(d, extra)
  w <- c(rep(1, n), rep(0, 30))
  spA <- statmod_spec(y ~ lasso(x, n_lambda = 8), g1, dA, weights = w)
  a <- cv_curve(spA, dA, iwls(), list(hy),
                c(fo, rep(1:5, length.out = 30)))
  b <- cv_curve(sp, d, iwls(), list(hy), fo)
  expect_equal(a$cvm, b$cvm, tolerance = 1e-8)
})

test_that("a restart is compared on the criterion that produced the value", {
  old <- list(criterion = 200, value = 2, criterion_by = bic())
  new <- list(criterion = 100, value = 1, criterion_by = bic())
  expect_true(edge_restart_better(new, old, 200, reml(), bic()))
})

test_that("cv() is a sparse criterion and never an outer one", {
  set.seed(2)
  d3 <- data.frame(x = stats::runif(150), X1 = stats::rnorm(150),
                   X2 = stats::rnorm(150))
  d3$y <- sin(6 * d3$x) + 0.8 * d3$X1 + stats::rnorm(150, sd = 0.3)
  expect_error(statmod(y ~ lasso(cbind(X1, X2)), g1, d3,
                       outer_criterion = cv(nfolds = 3)),
               "sparse_criterion = cv\\(\\)")
  expect_error(statmod(y ~ s(x) + lasso(cbind(X1, X2)), g1, d3,
                       outer_criterion = cv(nfolds = 3),
                       sparse_criterion = NULL),
               "sparse_criterion = cv\\(\\)")
  expect_error(statmod(y ~ s(x), g1, d3, outer_criterion = cv()),
               "nothing to select")
  expect_identical(names(formals(path_rows)), c("spec", "blocks", "hyper"))
})

test_that("a mixed covariance class counts its coordinates", {
  set.seed(77)
  m <- 8L
  ni <- 30L
  g <- factor(rep(seq_len(m), each = ni))
  dd <- data.frame(g = g, y = stats::rnorm(m * ni))
  f <- y ~ random(~ 1 | u | g) +
    gas(p = 1, q = 1, by = g, kappa1 ~ 1 + random(~ 1 | u | g))
  sp <- statmod_spec(f, g1, dd)
  des <- statmod_design(sp)
  u <- Filter(function(z) isTRUE(z$mixed), statmod_penalized(sp, des))[[1L]]
  expect_length(unit_positions(u), 16L)
  expect_identical(unit_positions(u), as.integer(u$joint))
})

test_that("start_search() names a coefficient by its label", {
  set.seed(2)
  d <- data.frame(x = stats::runif(100))
  d$y <- 1 + d$x + stats::rnorm(100, sd = 0.3)
  f <- statmod(y ~ x, g1, d,
               start = start_search(optimizers7::sa(maxit = 20),
                                    over = "sigma:(Intercept)"))
  expect_true(f@converged)
  expect_error(statmod(y ~ x, g1, d,
                       start = start_search(optimizers7::sa(maxit = 20),
                                            over = "(Intercept)")),
               "mu:\\(Intercept\\), sigma:\\(Intercept\\)")
  sp <- statmod_spec(y ~ x, g1, d)
  expect_error(start_at(start_search(), sp, statmod_design(sp), NULL),
               "evaluates the objective")
})

test_that("start_from() projects on the predictor beside a moving block", {
  set.seed(2)
  ds <- data.frame(x = stats::runif(200), z = stats::runif(200))
  ds$y <- 1 + 2 * pmax(ds$x - 0.5, 0) + sin(4 * ds$z) +
    stats::rnorm(200, sd = 0.2)
  r6 <- statmod(y ~ s(z, bspline_smooth(k = 6)) + seg(x), g1, ds)
  eta <- function(sp, des, cf) {
    d <- statmod_design_at(sp, cf, des)$mu
    as.numeric(d$X %*% cf$mu) + if (is.null(d$adj)) 0 else d$adj
  }
  sp <- statmod_spec(y ~ s(z, bspline_smooth(k = 10)) + seg(x), g1, ds)
  de <- statmod_design(sp)
  st <- start_at(start_from(r6), sp, de, NULL)
  expect_true(any(attr(st, "taken")$how == "projected"))
  expect_lt(max(abs(eta(sp, de, st) -
                      eta(r6@spec, statmod_design(r6@spec), r6@coefficients))),
            0.01)
  # a pending block that moves with its coefficients is left to the fallback
  sp2 <- statmod_spec(y ~ s(z, bspline_smooth(k = 6)) + seg(x, npsi = 2),
                      g1, ds)
  st2 <- start_at(start_from(r6), sp2, statmod_design(sp2), NULL)
  expect_false(any(grepl("^seg", attr(st2, "taken")$term)))
})

test_that("the working fits of a frozen block report what the pivot dropped", {
  set.seed(4)
  dj <- data.frame(x = stats::runif(200))
  dj$x2 <- 2 * dj$x
  dj$y <- 1 + 1.5 * (dj$x > 0.5) + dj$x + stats::rnorm(200, sd = 0.2)
  # the pivot's report is marked on the first coordinate of every smooth
  # fit, and each working phase has to hand it on
  seen <- list()
  fw <- fit_working
  fs <- fit_smooth
  testthat::local_mocked_bindings(
    fit_smooth = function(obj, beta, idx, ...) {
      r <- fs(obj, beta, idx, ...)
      r$aliased <- union(r$aliased, idx[[1L]])
      r
    },
    fit_working = function(...) {
      r <- fw(...)
      seen[[length(seen) + 1L]] <<- r$aliased
      r
    })
  f <- statmod(y ~ x + x2 + jump(x), g1, dj, outer_criterion = NULL)
  expect_gt(length(seen), 0L)
  expect_true(all(vapply(seen, function(a) 1L %in% a, logical(1))))
})
