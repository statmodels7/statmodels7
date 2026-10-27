# The reading of a joint fit over the coefficients and a structural term's
# parameters, and the flag and the certificate at fits that stop away from a
# maximum. Before, newton() ended such a fit on a stop on the change of the
# parameters (converged = TRUE at a log-likelihood 473 below the maximum), the
# certificate read the decrement over the coefficients alone after
# deficient_coords() had dropped the term's parameters from an indefinite
# curvature, and a joint gradient that could not be evaluated was returned as
# zeros, which crit_grad() reads as a stationary point.

test_that("joint_reading() reads the decrement of a positive definite curvature", {
  K <- matrix(c(4, 1, 0, 1, 3, 0.5, 0, 0.5, 2), 3, 3)
  g <- c(1, -2, 0.5)
  r <- joint_reading(g, K, 1:3)
  expect_false(r$ascent)
  expect_equal(r$decrement, 0.5 * sum(g * solve(K, g)))
  # a subset of the coordinates
  r2 <- joint_reading(g, K, c(1L, 3L))
  expect_equal(r2$decrement, 0.5 * sum(g[c(1, 3)] * solve(K[c(1, 3), c(1, 3)], g[c(1, 3)])))
})

test_that("joint_reading() reads a curvature with a direction of ascent as such", {
  # a block with eigenvalues 4.5 and -1.5 beside an eigenvalue of 1e10: the
  # negative one is far beyond rounding, and the decrement over a subset of
  # the coordinates would read the point as nearly a maximum
  K <- diag(c(1e10, 5, 1, 2))
  K[3, 4] <- K[4, 3] <- 3
  g <- c(1, 1, 1, 1)
  r <- joint_reading(g, K, 1:4)
  expect_true(r$ascent)
  expect_true(is.na(r$decrement))
  expect_gt(r$ratio, curvature_ascent_limit())
  # the coordinates of the block alone, without the large one
  r2 <- joint_reading(g, K, 3:4)
  expect_true(r2$ascent)
})

test_that("joint_reading() leaves rounding in a singular curvature to deficient_coords()", {
  # a duplicated column: the information is singular and rounding can make
  # the smallest eigenvalue negative
  X <- cbind(1, c(0.3, 1.2, -0.4, 2.2, 0.9), c(0.3, 1.2, -0.4, 2.2, 0.9))
  K <- crossprod(X)
  K[3, 3] <- K[3, 3] * (1 - 1e-15)
  g <- c(0.2, 0.1, 0.1)
  r <- joint_reading(g, K, 1:3)
  expect_false(r$ascent)
  expect_true(is.finite(r$decrement))
})

test_that("joint_reading() leaves out a coordinate that is flat at an edge", {
  # a coordinate with a diagonal of -4e-13 and no gradient: the rest is read
  K <- diag(c(3, 2, -4e-13))
  g <- c(1, 1, 1e-9)
  r <- joint_reading(g, K, 1:3)
  expect_false(r$ascent)
  expect_equal(r$decrement, 0.5 * (1 / 3 + 1 / 2))
  # nothing finite to read is not a decrement of zero
  expect_true(is.na(joint_reading(rep(NA_real_, 3), K, 1:3)$decrement))
  # no coordinate to read from the start: nothing to improve
  expect_equal(joint_reading(g, K, integer(0))$decrement, 0)
})

test_that("joint_reading() reads a curvature with no positive diagonal as ascent", {
  # an exploded filter had six diagonals of -4e23 to -3e28: leaving every
  # coordinate out as flat read as a decrement of zero
  K <- -diag(c(1e3, 2, 5))
  g <- c(1e5, 1, 1)
  r <- joint_reading(g, K, 1:3)
  expect_true(r$ascent)
  expect_true(is.na(r$decrement))
  expect_gt(r$ratio, curvature_ascent_limit())
})

# the data of the fits below: base datasets
dax <- data.frame(day = 1:1859,
                  r = 100 * diff(log(as.numeric(datasets::EuStockMarkets[, "DAX"]))))
erupt <- data.frame(t = seq_along(datasets::faithful$eruptions),
                    y = datasets::faithful$eruptions)
belt2 <- data.frame(t = 1:192,
                    y = as.numeric(datasets::Seatbelts[, "DriversKilled"] /
                                     datasets::Seatbelts[, "drivers"]))

gas_fit <- function(name) {
  d7 <- distributions7::student_t1_distrib
  switch(name,
    dax_t_nu_d0 = statmod(r ~ 1 | sigma ~ gas(p = 1, q = 1, time = day) |
                            nu ~ gas(p = 1, q = 1, time = day),
                          distrib = d7(), data = dax),
    dax_t_nu_dh = statmod(r ~ 1 | sigma ~ gas(p = 1, q = 1, time = day,
                                              scaling = 0.5) |
                            nu ~ gas(p = 1, q = 1, time = day, scaling = 0.5),
                          distrib = d7(), data = dax),
    erupt_weib = statmod(y ~ gas(p = 1, q = 1, time = t) |
                           sigma ~ gas(p = 1, q = 1, time = t),
                         distrib = distributions7::weibull1_distrib(),
                         data = erupt),
    belt2_beta = statmod(y ~ gas(p = 1, q = 1, time = t) |
                           phi ~ gas(p = 1, q = 1, time = t),
                         distrib = distributions7::beta1_distrib(),
                         data = belt2))
}

# the start of the term's loadings and persistences, on the unconstrained scale
fit_from <- function(name, start) {
  base_start <- structural_zeta_start
  testthat::with_mocked_bindings(
    gas_fit(name),
    structural_zeta_start = function(term, target = NULL) {
      z <- base_start(term, target = target)
      z[names(start)] <- start
      z
    })
}

# the log-likelihood at the best point found for each model over a grid of
# starts of the filter (misure_gas_partenze/c1_census.R)
best_ll <- c(dax_t_nu_d0 = -2480.5442, dax_t_nu_dh = -2479.7165,
             erupt_weib = -413.3641, belt2_beta = 675.0132)

# fits from starts that stop away from a maximum: the model, and the start
stalled <- list(
  # the census of gas() fits of the battery, as numbered there
  k320 = list("dax_t_nu_d0", c(sigma.omega = 0, sigma.kappa1 = -3.5548020394027575,
    sigma.pacf1 = 3.8002011672501994, nu.omega = 0, nu.kappa1 = -1.0717170825601805,
    nu.pacf1 = 3.8002011672501994)),
  k326 = list("dax_t_nu_d0", c(sigma.omega = 0, sigma.kappa1 = -2.3025850929940455,
    sigma.pacf1 = 3.8002011672501994, nu.omega = 0, nu.kappa1 = -2.3025850929940455,
    nu.pacf1 = 3.8002011672501994)),
  k357 = list("dax_t_nu_dh", c(sigma.omega = 0, sigma.kappa1 = -4.6819725148583924,
    sigma.pacf1 = 3.8002011672501994, nu.omega = 0, nu.kappa1 = -3.4404300364371037,
    nu.pacf1 = 3.8002011672501994)),
  k974 = list("erupt_weib", c(mu.omega = 0, mu.kappa1 = -1.2039728043259359,
    mu.pacf1 = 1.4722194895832204, sigma.omega = 0, sigma.kappa1 = -1.2039728043259359,
    sigma.pacf1 = 1.4722194895832204)),
  k975 = list("erupt_weib", c(mu.omega = 0, mu.kappa1 = -1.2039728043259359,
    mu.pacf1 = 2.0922957200349388, sigma.omega = 0, sigma.kappa1 = -1.2039728043259359,
    sigma.pacf1 = 2.0922957200349388)),
  k998 = list("belt2_beta", c(mu.omega = 0, mu.kappa1 = -3.4011973816621555,
    mu.pacf1 = 0.54930614433405489, phi.omega = 0, phi.kappa1 = -3.4011973816621555,
    phi.pacf1 = 0.54930614433405489)),
  k1002 = list("belt2_beta", c(mu.omega = 0, mu.kappa1 = -3.4011973816621555,
    mu.pacf1 = 3.8002011672501994, phi.omega = 0, phi.kappa1 = -3.4011973816621555,
    phi.pacf1 = 3.8002011672501994)),
  k1021 = list("belt2_beta", c(mu.omega = 0, mu.kappa1 = -1.2039728043259359,
    mu.pacf1 = 0, phi.omega = 0, phi.kappa1 = -1.2039728043259359, phi.pacf1 = 0)),
  k1025 = list("belt2_beta", c(mu.omega = 0, mu.kappa1 = -1.2039728043259359,
    mu.pacf1 = 2.6466524123622457, phi.omega = 0, phi.kappa1 = -1.2039728043259359,
    phi.pacf1 = 2.6466524123622457)),
  # exploded filters whose curvature is negative in every coordinate
  k1011 = list("belt2_beta", c(mu.omega = 0, mu.kappa1 = -2.3025850929940455,
    mu.pacf1 = 1.4722194895832204, phi.omega = 0, phi.kappa1 = -2.3025850929940455,
    phi.pacf1 = 1.4722194895832204)),
  k1023 = list("belt2_beta", c(mu.omega = 0, mu.kappa1 = -1.2039728043259359,
    mu.pacf1 = 1.4722194895832204, phi.omega = 0, phi.kappa1 = -1.2039728043259359,
    phi.pacf1 = 1.4722194895832204)),
  k1024 = list("belt2_beta", c(mu.omega = 0, mu.kappa1 = -1.2039728043259359,
    mu.pacf1 = 2.0922957200349388, phi.omega = 0, phi.kappa1 = -1.2039728043259359,
    phi.pacf1 = 2.0922957200349388)),
  # the loadings 0.05 and the persistences 0.9 on the share of drivers killed
  b4 = list("belt2_beta", c(mu.omega = 0, mu.kappa1 = -2.9957322735539909,
    mu.pacf1 = 1.4722194895832204, phi.omega = 0, phi.kappa1 = -2.9957322735539909,
    phi.pacf1 = 1.4722194895832204)))

test_that("a fit that stops away from a maximum is not reported as converged", {
  skip_on_cran()
  for (nm in names(stalled)) {
    model <- stalled[[nm]][[1L]]
    fit <- suppressWarnings(fit_from(model, stalled[[nm]][[2L]]))
    gap <- best_ll[[model]] - as.numeric(logLik(fit))
    # either the fit reached the maximum or it says that it did not
    info <- sprintf("%s: %.4g below the maximum, converged = %s", nm, gap,
                    fit@converged)
    expect_true(!fit@converged || gap < 1, info = info)
    cert <- statmod_certificate(fit)
    expect_true(!(cert$state %in% c("converged", "boundary")) || gap < 1,
                info = paste(nm, cert$state))
  }
})

# the design at the fitted structural state
fitted_pieces_joint <- function(fit) {
  sp <- fit@spec
  des <- statmod_design(sp)
  sst <- statmod_structural_state(des)
  for (k in names(fit@structural)) {
    sst$zeta[[k]] <- fit@structural[[k]]$unconstrained
    sst$held[[k]] <- fit@structural[[k]]$held
  }
  obj <- statmod_objective(sp, fit@hyper, des, FALSE, NULL)
  list(spec = sp, des = des, obj = obj, b = obj$stack(fit@coefficients))
}

test_that("the certificate reads an ascent direction of the joint curvature", {
  skip_on_cran()
  fit <- suppressWarnings(fit_from("belt2_beta", stalled$b4[[2L]]))
  gap <- best_ll[["belt2_beta"]] - as.numeric(logLik(fit))
  skip_if(gap < 1, "the fit reached the maximum on this platform")
  cert <- statmod_certificate(fit)
  expect_identical(cert$state, "not converged")
  p <- fitted_pieces_joint(fit)
  jp <- statmod_joint_pieces(p$spec, p$des, p$obj, fit@hyper,
                             kinds = c("filter", "loglik"))
  u <- c(p$b, jp$zeta())
  r <- joint_reading(jp$gr(u), as.matrix(jp$he(u)), seq_along(u))
  if (r$ascent) {
    expect_match(paste(cert$reason, collapse = " "), "direction of ascent")
    expect_true(is.na(cert$mode_error))
    # the reading that decides an alternation is the same
    expect_identical(
      alternation_readings(p$spec, p$des, p$obj, p$b, fit@hyper, FALSE,
                           NULL)[["mode"]], Inf)
    expect_false(alternation_settled(p$spec, p$des, p$obj, p$b, fit@hyper,
                                     FALSE, NULL))
  } else {
    expect_gt(cert$mode_error, mode_error_limit())
  }
})

test_that("the joint gradient is missing where it cannot be evaluated", {
  skip_on_cran()
  fit <- suppressWarnings(fit_from("belt2_beta", stalled$k1021[[2L]]))
  p <- fitted_pieces_joint(fit)
  jp <- statmod_joint_pieces(p$spec, p$des, p$obj, fit@hyper,
                             kinds = c("filter", "loglik"))
  u <- c(p$b, jp$zeta())
  # at the end of this fit the score may not be evaluable, and then it is
  # missing, never zero
  g0 <- jp$gr(u)
  expect_false(isTRUE(all(g0 == 0)))
  # a score that cannot be computed is missing, where zeros would be a
  # stationary point to crit_grad()
  g <- testthat::with_mocked_bindings(
    jp$gr(u),
    statmod_structural_score = function(...) stop("overflow"))
  expect_length(g, length(u))
  expect_true(all(is.na(g)))
  expect_true(is.na(joint_reading(g, diag(length(u)), seq_along(u))$decrement))
})
