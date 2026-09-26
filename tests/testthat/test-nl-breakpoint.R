# A parameter of nl() developed over a break-point term, fitted end to end,
# and the rule that names a coefficient whose own design column vanished.

test_that("a coefficient whose design column vanished is named, and only it", {
  # the rule is pinned on a design built for it rather than on a fit that
  # runs away, where a break-point lands is platform arithmetic
  set.seed(3)
  n <- 80
  d <- data.frame(x = stats::rnorm(n), z0 = 0, g = factor(rep(1:4, 20)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(n, sd = 0.3)
  f <- suppressWarnings(statmod(y ~ x + z0, distributions7::gaussian1_distrib(),
                                d, outer_criterion = NULL))
  spec <- f@spec
  design <- statmod_design(spec)
  cf <- f@coefficients
  K <- statmod_penalized_at(spec, cf, design, f@hyper)
  # the empty column of z0 is named, the live ones are not
  expect_identical(vanished_coords(spec, cf, design, K), 3L)
  expect_true("mu:z0" %in% f@aliased)
  # AT A BOUNDARY the design column is alive and only the information is
  # empty: that is uninformative_coords()' case and is not named here
  Kb <- as.matrix(K)
  Kb[2L, ] <- 0
  Kb[, 2L] <- 0
  Kb[3L, 3L] <- 1
  expect_identical(vanished_coords(spec, cf, design, Kb), integer(0))
  # a penalized coordinate with an empty column is identified by its prior
  fr <- suppressWarnings(statmod(y ~ x + ridge(~ 0 + z0, lambda = 2),
                                 distributions7::gaussian1_distrib(), d,
                                 outer_criterion = NULL))
  dr <- statmod_design(fr@spec)
  Kr <- statmod_penalized_at(fr@spec, fr@coefficients, dr, fr@hyper)
  expect_identical(vanished_coords(fr@spec, fr@coefficients, dr, Kr),
                   integer(0))
})

test_that("a parameter that steps at a break-point is fitted end to end", {
  set.seed(4)
  n <- 600
  d <- data.frame(t = seq(0, 10, length.out = n), x = stats::runif(n, 0, 3))
  lr <- log(0.5) + 0.875 * (d$t > 6)
  d$y <- 5 * exp(-exp(lr) * d$x) + stats::rnorm(n, sd = 0.2)
  f <- statmod(y ~ 0 + nl(~ a * exp(-r * x),
                          r ~ jump(t, smoothed = numericals7::smooth_probit(h = 0.05)),
                          links = list(r = linkfunctions7::log_link()),
                          start = list(a = 4)),
               distributions7::gaussian1_distrib(), d)
  expect_true(f@converged)
  cf <- stats::setNames(f@coefficients$mu,
                        statmod_design(f@spec)$mu$coef_names)
  expect_lt(abs(cf[[grep("psi", names(cf))]] - 6), 0.1)
  expect_lt(abs(cf[["nl.a"]] - 5), 0.1)
  expect_length(f@aliased, 0L)
  # the fitted term carries its sub-term at the fitted position
  tm <- f@spec@terms$mu[[1L]]
  expect_equal(as.numeric(modelterms7::seg_psi(
    modelterms7::term_components(tm)$r$subs[[2L]])),
    cf[[grep("psi", names(cf))]], tolerance = 1e-12)
})
