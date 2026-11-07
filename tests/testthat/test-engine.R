# The compiled inner fit (R/engine.R, src/engine_iwls.cpp) against
# iwls_fit(): a model fitted with engine = "auto" and with engine = "r" must
# reach the same point by the same iterations, since the two routes perform
# the same operations in the same order.

engine_twin <- function(formula, distrib, data, ...) {
  fa <- statmod(formula, distrib = distrib, data = data, engine = "auto", ...)
  fr <- statmod(formula, distrib = distrib, data = data, engine = "r", ...)
  expect_equal(unlist(coef(fa)), unlist(coef(fr)), tolerance = 1e-10)
  expect_equal(as.numeric(logLik(fa)), as.numeric(logLik(fr)),
               tolerance = 1e-10)
  invisible(list(auto = fa, r = fr))
}

engine_data <- function(n = 400, seed = 3) {
  set.seed(seed)
  d <- data.frame(x = runif(n), z = runif(n),
                  g = factor(sample(12, n, replace = TRUE)))
  u <- rnorm(12, 0, 0.4)
  eta <- sin(5 * d$x) + 0.5 * d$z + u[as.integer(d$g)]
  d$y <- eta + rnorm(n, 0, 0.5)
  d$yg <- stats::rgamma(n, shape = 3, rate = 3 / exp(eta))
  d$yc <- stats::rpois(n, exp(0.3 + eta))
  d
}

# the inputs of engine_route() for a model, as statmod() builds them
engine_route_of <- function(formula, distrib, data, method = iwls()) {
  spec <- statmodels7:::statmod_spec(formula, distrib, data)
  design <- statmodels7:::statmod_design(spec)
  hyper <- statmodels7:::statmod_hyper_start(spec, design)
  method <- statmodels7:::inner_bind(method, spec@distrib, auto_info())
  obj <- statmodels7:::statmod_objective(spec, hyper, design)
  beta <- statmodels7:::statmod_start(spec, design, obj, NULL)
  statmodels7:::engine_route(spec, design, obj, beta, hyper, method)
}

test_that("a smooth and a random intercept are covered and fit alike", {
  d <- engine_data()
  f <- y ~ s(x, bspline_smooth(k = 10)) + random(~ 1 | g)
  expect_false(is.null(engine_route_of(f, distributions7::gaussian1_distrib(), d)))
  engine_twin(f, distributions7::gaussian1_distrib(), d)
})

test_that("a two-parameter family with a modelled dispersion fits alike", {
  d <- engine_data()
  engine_twin(yg ~ s(x, bspline_smooth(k = 8)) | phi ~ z,
              distributions7::gamma1_distrib(), d)
})

test_that("the three solvers fit alike", {
  d <- engine_data()
  f <- yc ~ x + s(z, bspline_smooth(k = 8))
  for (sv in list(solve_qr(), solve_svd(), solve_chol())) {
    engine_twin(f, distributions7::poisson_distrib(), d,
                inner_optimizer = iwls(solver = sv))
  }
})

test_that("a named information fits alike", {
  d <- engine_data()
  f <- yg ~ x + z
  for (inf in list(observed_info(), expected_info())) {
    engine_twin(f, distributions7::gamma1_distrib(), d,
                inner_information = inf)
  }
})

test_that("a kinked penalty is covered and a user criterion is left to iwls_fit()", {
  d <- engine_data()
  # since phase 3b the kinked units have a compiled route (test-engine-kink.R)
  expect_false(is.null(engine_route_of(y ~ lasso(~ x + z),
                                       distributions7::gaussian1_distrib(), d)))
  expect_null(engine_route_of(
    y ~ x, distributions7::gaussian1_distrib(), d,
    method = iwls(criterion = optimizers7::crit_grad(1e-6))))
})

test_that("a simulated or integrated expected information fits alike", {
  d <- engine_data(n = 150)
  # the simulation draws from R's stream at every evaluation: the same seed
  # before each route, and the same state of the stream after
  twin_seeded <- function(...) {
    set.seed(4)
    fa <- statmod(..., engine = "auto")
    sa <- stats::runif(1)
    set.seed(4)
    fr <- statmod(..., engine = "r")
    sr <- stats::runif(1)
    expect_equal(unlist(coef(fa)), unlist(coef(fr)), tolerance = 1e-10)
    expect_identical(sa, sr)
  }
  for (m in list(distributions7::expect_mc(nsim = 200),
                 distributions7::expect_mc(nsim = 200, identity = "score"),
                 distributions7::expect_integrate(),
                 distributions7::expect_integrate(identity = "score"))) {
    for (inf in list(expected_info(method = m), auto_info(method = m))) {
      twin_seeded(yg ~ x + z, distrib = distributions7::gamma1_distrib(),
                  data = d, inner_information = inf)
      twin_seeded(yc ~ x, distrib = distributions7::poisson_distrib(),
                  data = d, inner_information = inf)
    }
  }
})

test_that("the outer product of the scores fits alike", {
  d <- engine_data()
  f <- yg ~ x + s(z, bspline_smooth(k = 6))
  # named, every iteration on it; and as the expectation of auto_info(),
  # where it is the start of the run or the backup of the observed step
  for (inf in list(expected_info(method = distributions7::expect_opg()),
                   auto_info(method = distributions7::expect_opg()))) {
    engine_twin(f, distributions7::gamma1_distrib(), d,
                inner_information = inf)
    engine_twin(yc ~ x + z, distributions7::poisson_distrib(), d,
                inner_information = inf)
  }
})
