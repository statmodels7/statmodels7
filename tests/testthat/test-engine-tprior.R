# The compiled fit of a model with a univariate Student t prior on a random
# effect (unit kind 9 of engine_penalty_units()), whose Hessian moves with the
# coefficients, against the R route: the value and the gradient of the
# penalty at a point, and the inner fit with the penalty's square root formed
# at every point, which takes the assembled route wherever the prior's
# curvature is negative.

tprior_data <- function(n = 240, G = 20, seed = 7) {
  set.seed(seed)
  d <- data.frame(g = factor(sample(G, n, replace = TRUE)), x = stats::rnorm(n))
  b <- stats::rt(G, df = 2) * 0.5
  d$y <- 1 + 0.5 * d$x + b[as.integer(d$g)] + stats::rnorm(n, 0, 0.4)
  d$yc <- stats::rpois(n, exp(0.2 + 0.3 * d$x + b[as.integer(d$g)] / 2))
  d$yb <- stats::rbinom(n, 1, stats::plogis(0.3 * d$x + b[as.integer(d$g)]))
  d
}
t_free <- function() distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
t_held <- function() distributions7::fixed(distributions7::student_t1_distrib(), mu = 0,
                                            nu = 1e6)

# fits the model on both routes, checks that the compiled inner fit ran, and
# compares the point
tprior_twin <- function(formula, distrib, data, ...) {
  hit <- 0L
  orig <- statmodels7:::engine_inner_fit
  local_mocked_bindings(engine_inner_fit = function(...) {
    r <- orig(...)
    if (!is.null(r)) hit <<- hit + 1L
    r
  })
  fa <- statmod(formula, distrib = distrib, data = data, engine = "auto", ...)
  used <- hit
  fr <- statmod(formula, distrib = distrib, data = data, engine = "r", ...)
  expect_gt(used, 0L)
  expect_equal(unlist(coef(fa)), unlist(coef(fr)), tolerance = 1e-10)
  expect_equal(as.numeric(logLik(fa)), as.numeric(logLik(fr)), tolerance = 1e-10)
  expect_equal(unlist(fa@hyper), unlist(fr@hyper), tolerance = 1e-8)
  invisible(list(auto = fa, r = fr))
}

test_that("the compiled value and gradient of a t prior are R's", {
  d <- tprior_data()
  for (pr in list(t_free(), t_held())) {
    f <- y ~ x + random(~ 1 | g, distrib = pr)
    spec <- statmodels7:::statmod_spec(f, distributions7::gaussian1_distrib(), d)
    design <- statmodels7:::statmod_design(spec)
    hyper <- statmodels7:::statmod_hyper_start(spec, design)
    obj <- statmodels7:::statmod_objective(spec, hyper, design)
    beta <- statmodels7:::statmod_start(spec, design, obj, NULL)
    em <- statmodels7:::engine_model(spec, design, obj$split(beta), hyper, "qr")
    expect_false(is.null(em))
    expect_true(any(vapply(em$pen$units, function(u) identical(u$kind, 9L), logical(1))))
    set.seed(5)
    for (r in 1:10) {
      # coefficients from small to far in the tail, where the t's curvature
      # in the location is negative and its overflow guards are taken
      b <- beta + stats::rnorm(length(beta), 0, 10^(r %% 5 - 1))
      e <- statmodels7:::engine_point_probe(em$model, em$pen, b)
      expect_identical(e$value, obj$fn(b))
      expect_identical(as.numeric(e$gradient), as.numeric(obj$gr(b)))
    }
  }
})

test_that("a t prior with its degrees of freedom estimated fits alike", {
  d <- tprior_data()
  tprior_twin(y ~ x + random(~ 1 | g, distrib = t_free()),
              distributions7::gaussian1_distrib(), d)
})

test_that("a t prior with its degrees of freedom held fits alike", {
  d <- tprior_data()
  tprior_twin(y ~ x + random(~ 1 | g, distrib = t_held()),
              distributions7::gaussian1_distrib(), d)
})

test_that("a t prior under ml() fits alike", {
  d <- tprior_data()
  tprior_twin(y ~ x + random(~ 1 | g, distrib = t_free()),
              distributions7::gaussian1_distrib(), d, outer_criterion = ml())
})

test_that("a t prior on a poisson and a bernoulli response fits alike", {
  d <- tprior_data()
  tprior_twin(yc ~ x + random(~ 1 | g, distrib = t_free()),
              distributions7::poisson_distrib(), d)
  tprior_twin(yb ~ x + random(~ 1 | g, distrib = t_free()),
              distributions7::bernoulli_distrib(), d)
})

test_that("a t prior beside a smooth fits alike", {
  d <- tprior_data()
  d$z <- stats::runif(nrow(d))
  d$y2 <- d$y + sin(5 * d$z)
  tprior_twin(y2 ~ x + s(z) + random(~ 1 | g, distrib = t_free()),
              distributions7::gaussian1_distrib(), d)
})

# a multivariate t prior (unit kind 10): the 2-dimensional prior on the
# default log-Cholesky chart, and an AR(1) prior of dimension 5 on a panel
# observed at five times
mvt_data <- function(seed = 7L, m = 30L, ni = 8L) {
  set.seed(seed)
  g <- factor(rep(seq_len(m), each = ni))
  x <- rep(seq(-1, 1, length.out = ni), m)
  b <- matrix(stats::rnorm(2 * m, 0, 0.4), m, 2)
  b[1:3, ] <- c(2.5, -2.4, 2.8, 1.5, -1.8, 1.2)
  data.frame(g = g, x = x,
             y = 1 + 0.5 * x + b[as.integer(g), 1] + b[as.integer(g), 2] * x +
               stats::rnorm(m * ni, 0, 0.5))
}
ts5_data <- function(seed = 9L, m = 16L) {
  set.seed(seed)
  id <- factor(rep(seq_len(m), each = 10))
  time <- rep(rep(1:5, each = 2), m)
  b <- matrix(stats::rt(5 * m, df = 4) * 0.5, m, 5)
  d <- data.frame(id = id, time_g = factor(time))
  d$y <- 1 + 0.2 * time + b[cbind(as.integer(id), time)] + stats::rnorm(nrow(d), 0, 0.3)
  d
}
mvt2 <- function() distributions7::fixed(distributions7::mvstudent_t1_distrib(2),
                                          mu1 = 0, mu2 = 0)
ts5 <- function() distributions7::fixed(distributions7::mvstudent_t1_distrib(
  n_dim = 5, sigma = parameters7::ar1(5)), mu = 0)

test_that("the compiled value and gradient of a multivariate t prior are R's", {
  cases <- list(list(y ~ x + random(~ x | g, distrib = mvt2()), mvt_data()),
                list(y ~ 1 + random(~ 0 + time_g | id, distrib = ts5()), ts5_data()))
  for (cs in cases) {
    spec <- statmodels7:::statmod_spec(cs[[1]], distributions7::gaussian1_distrib(), cs[[2]])
    design <- statmodels7:::statmod_design(spec)
    hyper <- statmodels7:::statmod_hyper_start(spec, design)
    obj <- statmodels7:::statmod_objective(spec, hyper, design)
    beta <- statmodels7:::statmod_start(spec, design, obj, NULL)
    em <- statmodels7:::engine_model(spec, design, obj$split(beta), hyper, "qr")
    expect_false(is.null(em))
    expect_true(any(vapply(em$pen$units, function(u) identical(u$kind, 10L), logical(1))))
    set.seed(5)
    for (r in 1:10) {
      b <- beta + stats::rnorm(length(beta), 0, 10^(r %% 5 - 1))
      e <- statmodels7:::engine_point_probe(em$model, em$pen, b)
      expect_identical(e$value, obj$fn(b))
      expect_identical(as.numeric(e$gradient), as.numeric(obj$gr(b)))
    }
  }
})

test_that("a multivariate t prior on the log-Cholesky chart fits alike", {
  tprior_twin(y ~ x + random(~ x | g, distrib = mvt2()),
              distributions7::gaussian1_distrib(), mvt_data())
})

test_that("an AR(1) multivariate t prior on a sparse design fits alike", {
  # 24 subjects make the design sparse: S and C move on a sparse S
  d <- ts5_data(m = 24L)
  spec <- statmodels7:::statmod_spec(y ~ 1 + random(~ 0 + time_g | id, distrib = ts5()),
                                     distributions7::gaussian1_distrib(), d)
  expect_true(statmodels7:::design_sparse(statmodels7:::statmod_design(spec)))
  tprior_twin(y ~ 1 + random(~ 0 + time_g | id, distrib = ts5()),
              distributions7::gaussian1_distrib(), d)
})

test_that("an AR(1) multivariate t prior of dimension 5 fits alike", {
  tprior_twin(y ~ 1 + random(~ 0 + time_g | id, distrib = ts5()),
              distributions7::gaussian1_distrib(), ts5_data())
})
