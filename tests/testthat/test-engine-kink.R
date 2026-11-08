# The compiled alternation of a model with a kinked penalty (R/engine_kink.R,
# src/engine_kink.h) against statmod_alternate() in R: a model fitted with
# engine = "auto" and with engine = "r" reaches the same point by the same
# passes, block by block, since the two routes perform the same operations in
# the same order.

kink_data <- function(n = 300, seed = 11) {
  set.seed(seed)
  d <- data.frame(x = runif(n), z = runif(n), w = rnorm(n), v = rnorm(n, 2, 3),
                  u = rnorm(n), g = factor(sample(10, n, replace = TRUE)))
  b <- rnorm(10, 0, 0.4)
  eta <- sin(5 * d$x) + 0.6 * d$z - 0.3 * d$w + b[as.integer(d$g)]
  d$y <- eta + rnorm(n, 0, 0.5)
  d$yg <- stats::rgamma(n, shape = 3, rate = 3 / exp(0.2 + 0.4 * d$z - 0.2 * d$w))
  d$yc <- stats::rpois(n, exp(0.3 + 0.5 * d$z - 0.4 * d$w))
  d
}

# fits the model on both routes, checks that the compiled alternation ran,
# and compares the point and the histories
kink_twin <- function(formula, distrib, data, ...) {
  hit <- 0L
  orig <- statmodels7:::engine_alternate_fit
  # the compiled outer search runs the compiled alternation from C++, the
  # final refit included, so a search with a compiled point counts too
  orig_s <- statmodels7:::engine_outer_search
  local_mocked_bindings(engine_alternate_fit = function(...) {
    r <- orig(...)
    if (!is.null(r)) hit <<- hit + 1L
    r
  }, engine_outer_search = function(cfg, cb) {
    if (!is.null(cfg$compiled)) hit <<- hit + 1L
    orig_s(cfg, cb)
  })
  fa <- statmod(formula, distrib = distrib, data = data, engine = "auto", ...)
  used <- hit
  fr <- statmod(formula, distrib = distrib, data = data, engine = "r", ...)
  expect_gt(used, 0L)
  expect_equal(unlist(coef(fa)), unlist(coef(fr)), tolerance = 1e-10)
  expect_equal(as.numeric(logLik(fa)), as.numeric(logLik(fr)), tolerance = 1e-10)
  ba <- fa@history$blocks
  br <- fr@history$blocks
  expect_identical(nrow(ba), nrow(br))
  expect_identical(ba$pass, br$pass)
  expect_identical(ba$block, br$block)
  expect_identical(as.integer(ba$iterations), as.integer(br$iterations))
  expect_identical(ba$converged, br$converged)
  expect_equal(ba$objective, br$objective, tolerance = 1e-10)
  ia <- fa@history$inner
  ir <- fr@history$inner
  expect_identical(is.null(ia), is.null(ir))
  if (!is.null(ia)) {
    expect_identical(as.integer(ia$iteration), as.integer(ir$iteration))
    expect_identical(ia$route, ir$route)
  }
  invisible(list(auto = fa, r = fr, used = used))
}

test_that("the compiled value and gradient of the kinked units are R's", {
  d <- kink_data()
  forms <- list(
    y ~ s(x) + lasso(~ z + w + v, lambda = 2),
    y ~ s(x) + lasso(~ z + w + v, lambda = 2, standardize = TRUE),
    y ~ enet(~ z + w + v, lambda = 3, alpha = 0.4, standardize = TRUE),
    y ~ scad(~ z + w + v, lambda = 0.3, a = 3.7, standardize = TRUE),
    y ~ mcp(~ z + w + v, lambda = 0.2, gamma = 3))
  for (f in forms) {
    spec <- statmodels7:::statmod_spec(f, distributions7::gaussian1_distrib(), d)
    design <- statmodels7:::statmod_design(spec)
    hyper <- statmodels7:::statmod_hyper_start(spec, design)
    obj <- statmodels7:::statmod_objective(spec, hyper, design)
    beta <- statmodels7:::statmod_start(spec, design, obj, NULL)
    em <- statmodels7:::engine_model(spec, design, obj$split(beta), hyper, "chol")
    expect_false(is.null(em))
    set.seed(5)
    for (r in 1:10) {
      b <- beta + rnorm(length(beta), 0, 0.3)
      if (r %% 3 == 0) b[sample(length(b), 2)] <- 0
      e <- statmodels7:::engine_point_probe(em$model, em$pen, b)
      expect_identical(e$value, obj$fn(b))
      expect_identical(as.numeric(e$gradient), as.numeric(obj$gr(b)))
    }
  }
})

test_that("a lasso beside a smooth, its smoothing parameter estimated, fits alike", {
  d <- kink_data()
  kink_twin(y ~ s(x) + lasso(~ z + w + v + u, lambda = 20),
            distributions7::gaussian1_distrib(), d)
})

test_that("a standardized lasso on a gamma mean fits alike", {
  d <- kink_data()
  kink_twin(yg ~ lasso(~ z + w + v + u, lambda = 15, standardize = TRUE),
            distributions7::gamma1_distrib(), d)
})

test_that("a lasso with no smooth block and no intercept fits alike", {
  d <- kink_data()
  kink_twin(yc ~ 0 + lasso(~ z + w + v, lambda = 5),
            distributions7::poisson_distrib(), d)
})

test_that("two kinked blocks in two equations fit alike", {
  d <- kink_data()
  kink_twin(y ~ enet(~ z + w + v, lambda = 10, alpha = 0.5) | sigma ~ lasso(~ z + w, lambda = 5),
            distributions7::gaussian1_distrib(), d)
})

test_that("a scaled SCAD and an MCP beside a random intercept fit alike", {
  d <- kink_data()
  r <- kink_twin(y ~ scad(~ z + w + v + u, lambda = 20, a = 3.7, standardize = TRUE),
                 distributions7::gaussian1_distrib(), d)
  kink_twin(y ~ mcp(~ z + w + v, lambda = 10, gamma = 3) + random(~ 1 | g),
            distributions7::gaussian1_distrib(), d)
})

test_that("a negbin2 with a lasso on theta, on the observed information, fits alike", {
  # where the observed curvature of the theta block is not positive the
  # working weights fall back to the expected ones (test-coord-fallback.R)
  set.seed(3)
  nf <- 400L
  xf <- matrix(stats::rnorm(nf * 8L), nf, 8L)
  zf <- matrix(stats::rnorm(nf * 4L), nf, 4L)
  df <- data.frame(y = stats::rnbinom(nf, mu = exp(1 + xf[, 1] / 2 - xf[, 2] / 3),
                                       size = exp(1 + 0.5 * zf[, 1])))
  df$X <- xf
  df$Z <- zf
  kink_twin(y ~ lasso(X, lambda = 5) | theta ~ lasso(Z, lambda = 5),
            distributions7::negbin2_distrib(), df,
            inner_information = observed_info())
})

test_that("a lasso path chosen by bic() fits alike", {
  d <- kink_data()
  kink_twin(yg ~ lasso(~ z + w + v + u), distributions7::gamma1_distrib(), d)
})

test_that("a kinked block beside a nonlinear term stays on the R route", {
  d <- kink_data()
  hit <- 0L
  orig <- statmodels7:::engine_alternate_fit
  # the compiled outer search runs the compiled alternation from C++, the
  # final refit included, so a search with a compiled point counts too
  orig_s <- statmodels7:::engine_outer_search
  local_mocked_bindings(engine_alternate_fit = function(...) {
    r <- orig(...)
    if (!is.null(r)) hit <<- hit + 1L
    r
  }, engine_outer_search = function(cfg, cb) {
    if (!is.null(cfg$compiled)) hit <<- hit + 1L
    orig_s(cfg, cb)
  })
  d$yp <- 2 * exp(-1.5 * d$x) * (1 + 0.3 * as.integer(d$g) / 10) + rnorm(nrow(d), 0, 0.1)
  suppressWarnings(statmod(yp ~ 0 + nl(~ a * exp(-r * x), a ~ 1 + lasso(~ g, lambda = 0.5)),
                           distrib = distributions7::gaussian1_distrib(), data = d))
  expect_identical(hit, 0L)
})
