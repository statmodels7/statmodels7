test_that("a static start reads the constant from the model without the term", {
  # share ~ law + gas(): the constant is the predictor of share ~ law at the
  # first month, where law = 0
  sb <- data.frame(month = seq_len(nrow(datasets::Seatbelts)),
                   share = datasets::Seatbelts[, "front"] /
                     (datasets::Seatbelts[, "front"] +
                        datasets::Seatbelts[, "rear"]),
                   law = datasets::Seatbelts[, "law"])
  spec <- statmod_spec(share ~ law + gas(p = 1, q = 1, time = month,
                                         init = "static"),
                       distributions7::beta1_distrib(), sb)
  tm <- Filter(function(t) S7::S7_inherits(t, modelterms7::GasTerm),
               unlist(spec@terms))[[1L]]
  f0 <- statmod(share ~ law, distrib = distributions7::beta1_distrib(),
                data = sb)
  expect_equal(tm@init_value, f0@coefficients$mu[[1L]], tolerance = 1e-8)
  expect_identical(static_equation_rhs(quote(law + gas(p = 1, time = month))),
                   quote(law + 1))
  expect_identical(static_equation_rhs(quote(0 + x + gas(omega ~ 1 + z))),
                   quote(0 + x + (1 + z)))
  nile <- data.frame(year = 1871:1970, flow = as.numeric(Nile))
  # the default start leaves the constant unset
  spec0 <- statmod_spec(flow ~ gas(p = 1, q = 1, time = year),
                        distributions7::gaussian1_distrib(), nile)
  tm0 <- Filter(function(t) S7::S7_inherits(t, modelterms7::GasTerm),
                unlist(spec0@terms))[[1L]]
  expect_length(tm0@init_value, 0L)
})

test_that("the derivatives of a statically started filter are the objective's", {
  # f_1 = c - eta_1 carries the first row's static predictor into every later
  # level without a score in between; the gradient assembled from the
  # adjoint's score sequence missed it until the adjoint returned it
  nile <- data.frame(year = 1871:1970, flow = as.numeric(Nile))
  for (q in 1:2) {
    fit <- statmod(flow ~ gas(p = 1, q = q, time = year, init = "static"),
                   distrib = distributions7::gaussian1_distrib(), data = nile)
    expect_identical(statmod_certificate(fit)$state, "converged")
    spec <- fit@spec
    design <- statmod_design(spec)
    obj <- statmod_objective(spec, fit@hyper, design, FALSE, NULL)
    jp <- statmod_joint_pieces(spec, design, obj, fit@hyper)
    u <- c(obj$stack(fit@coefficients), jp$zeta())
    u <- u + 0.01 * seq_along(u) / length(u)
    expect_equal(jp$gr(u), numDeriv::grad(jp$fn, u), tolerance = 1e-6,
                 ignore_attr = TRUE)
    expect_equal(as.matrix(jp$he(u)), numDeriv::jacobian(jp$gr, u),
                 tolerance = 1e-8, ignore_attr = TRUE)
    # the first fitted predictor is the constant
    tm <- Filter(function(t) S7::S7_inherits(t, modelterms7::GasTerm),
                 unlist(spec@terms))[[1L]]
    expect_equal(fitted(fit, what = "mu")[1L], tm@init_value, tolerance = 1e-10)
  }
})

test_that("a static start removes the second maximum of the DAX log-scale", {
  # with the start at the level the intercept doubles as the starting value
  # at a persistence near one: from kappa 0.01, xi 0.9 the fit reached
  # -2591.37 (xi 0.99958) where the default start reaches -2616.35
  dax <- data.frame(day = 1:1859,
                    r = 100 * diff(log(as.numeric(EuStockMarkets[, "DAX"]))))
  form <- r ~ 1 | sigma ~ gas(p = 1, q = 1, time = day, init = "static")
  a <- statmod(form, distrib = distributions7::gaussian1_distrib(), data = dax)
  expect_identical(statmod_certificate(a)$state, "converged")
  expect_equal(as.numeric(logLik(a)), -2614.4019, tolerance = 1e-7)
  expect_lt(a@structural[[1L]]$parameter[["pacf1"]], 0.99)
  # the static start is not a parameter, so the model keeps four
  expect_equal(attr(logLik(a), "df"), 4)
})

test_that("a joint filter with static starts fits with its exact derivatives", {
  sb <- data.frame(t = seq_len(nrow(datasets::Seatbelts)),
                   y = as.numeric(datasets::Seatbelts[, "DriversKilled"] /
                                    datasets::Seatbelts[, "drivers"]))
  form <- y ~ gas(p = 1, q = 1, time = t, init = "static") |
    phi ~ gas(p = 1, q = 1, time = t, init = "static")
  fit <- statmod(form, distrib = distributions7::beta1_distrib(), data = sb)
  expect_identical(statmod_certificate(fit)$state, "converged")
  spec <- fit@spec
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, fit@hyper, design, FALSE, NULL)
  jp <- statmod_joint_pieces(spec, design, obj, fit@hyper)
  u <- c(obj$stack(fit@coefficients), jp$zeta())
  u <- u + 0.002 * seq_along(u) / length(u)
  expect_equal(jp$gr(u), numDeriv::grad(jp$fn, u), tolerance = 1e-6,
               ignore_attr = TRUE)
  expect_equal(as.matrix(jp$he(u)), numDeriv::jacobian(jp$gr, u),
               tolerance = 1e-8, ignore_attr = TRUE)
})
