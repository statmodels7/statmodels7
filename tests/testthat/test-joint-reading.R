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

# The joint gradient and curvature, over two coefficients and the four free
# parameters of a joint filter, at three points where the filter had exploded
# (misure_gas_flag2/r10_dput.R, census fits 997, 977 and 976 before the second
# start of a joint filter was repaired). The curvature is numerically of rank
# one, with an eigenvalue of 1e30 beside eigenvalues at rounding level.
exploded <- list(
  k997 = list(
    g = c(-759078761094.63, -283476044500.78375, -282508705919.03528,
          5066353445641.6699, -35469678242.702888, 224022121731.63144),
    K = matrix(c(4.5427607867081516e+28, 1.6964825611320266e+28,
      1.6906934476243765e+28, -3.0319952134725618e+29, 2.1227079924822597e+27,
      -1.3406763178780091e+28, 1.6964825611320266e+28, 6.335471347393225e+27,
      6.313852048972586e+27, -1.1322909672334523e+29, 7.9271994734092752e+26,
      -5.0067219036881945e+27, 1.6906934476243765e+28, 6.313852048972586e+27,
      6.2923065247099495e+27, -1.1284271132315421e+29, 7.9001485277639948e+26,
      -4.9896368583916838e+27, -3.0319952134725618e+29, -1.1322909672334523e+29,
      -1.1284271132315421e+29, 2.0236581687106832e+30, -1.4167685191871943e+28,
      8.9481360993426385e+28, 2.1227079924822597e+27, 7.9271994734092752e+26,
      7.9001485277639948e+26, -1.4167685191871943e+28, 9.9188344553207156e+25,
      -6.2646141166362033e+26, -1.3406763178780091e+28, -5.0067219036881945e+27,
      -4.9896368583916838e+27, 8.9481360993426385e+28, -6.2646141166362033e+26,
      3.9566533958338036e+27), 6, 6)),
  k977 = list(
    g = c(257289676931200.97, -205024693605584.22, -194253498263727.97,
          -1730875229438.3721, 21971661631146.711, -202350533797.46982),
    K = matrix(c(3.5640285886799963e+28, -2.8400434798283226e+28,
      -2.6908387056977846e+28, -2.3976433390733855e+26, 3.0435589615544509e+27,
      -2.8030004778725488e+25, -2.8400434798283226e+28, 2.2631263376882498e+28,
      2.1442305332958008e+28, 1.9105939143465865e+26, -2.4253003501963244e+27,
      2.2336081299746734e+25, -2.6908387056977846e+28, 2.1442305332958008e+28,
      2.0315810493436769e+28, 1.8102187843635298e+26, -2.2978845576142693e+27,
      2.1162630967397068e+25, -2.3976433390733855e+26, 1.9105939143465865e+26,
      1.8102187843635298e+26, 1.6129762818610601e+24, -2.0475057058817834e+25,
      1.8856738261126451e+23, 3.0435589615544509e+27, -2.4253003501963244e+27,
      -2.2978845576142693e+27, -2.0475057058817834e+25, 2.5990956362903193e+26,
      -2.3936668888580274e+24, -2.8030004778725488e+25, 2.2336081299746734e+25,
      2.1162630967397068e+25, 1.8856738261126451e+23, -2.3936668888580274e+24,
      2.2044749315564632e+22), 6, 6)),
  k976 = list(
    g = c(-1252.6142740797723, 605.91789179562204, 350.44668103802769,
          137.77459900563738, -1291.3739258739324, 232.42730347695883),
    K = matrix(c(1.9175484295919468e+18, -1.1672801699937669e+18,
      -1.616662109498297e+17, -5.4333506222585222e+17, 388004625946562.12,
      -26053022490607816, -1.1672801699937669e+18, 7.1056510189451853e+17,
      98411992775547840, 3.3074744502465005e+17, -236192264427771.25,
      15859404671286296, -1.616662109498297e+17, 98411992775547840,
      13629884575297434, 45807928201811728, -32712205199126.492,
      2196499115666937.8, -5.4333506222585222e+17, 3.3074744502465005e+17,
      45807928201811728, 1.5395334234494029e+17, -109940648345771.06,
      7382092873199624, 388004625946562.12, -236192264427771.25,
      -32712205199126.492, -109940648345771.06, 78510425008.165543,
      -5271675587996.6807, -26053022490607816, 15859404671286296,
      2196499115666937.8, 7382092873199624, -5271675587996.6807,
      353972796942626.38), 6, 6)))

test_that("joint_reading() does not read a point from a rank-one curvature that drops the term", {
  e <- exploded$k997
  # the aliasing test names four of the six coordinates, the term's own among
  # them; read over the intercept alone the decrement was 6.3e-6, 3128 below the
  # maximum
  r0 <- joint_reading(e$g, e$K, 1:6)
  expect_false(r0$partial)
  expect_lt(r0$decrement, mode_error_limit())
  r <- joint_reading(e$g, e$K, 1:6, own = 3:6)
  expect_true(r$partial)
  expect_false(r$ascent)
  expect_true(is.na(r$decrement))
})

test_that("a partial reading above the limit is kept as a lower bound", {
  e <- exploded$k977
  r <- joint_reading(e$g, e$K, 1:6, own = 3:6)
  expect_true(r$partial)
  expect_false(r$ascent)
  expect_gt(r$decrement, mode_error_limit())
})

test_that("a negative decrement is a direction of ascent", {
  e <- exploded$k976
  r <- joint_reading(e$g, e$K, 1:6, own = 3:6)
  expect_true(r$ascent)
  expect_true(is.na(r$decrement))
})

# filter_decrement() on a given gradient and curvature, over `npar` coefficients
fd_at <- function(g, K, npar = 2L) {
  testthat::with_mocked_bindings(
    filter_decrement(NULL, NULL, list(obj = 1, par = numeric(npar)), NULL),
    statmod_joint_pieces = function(...)
      list(zeta = function() numeric(length(g) - npar),
           gr = function(u) g, he = function(u) K),
    held_stack = function(...) integer(0))
}

test_that("filter_decrement() reads a positive definite curvature as the decrement", {
  K <- matrix(c(4, 1, 0, 1, 3, 0.5, 0, 0.5, 2), 3, 3)
  g <- c(1, -2, 0.5)
  expect_equal(fd_at(g, K, 1L), 0.5 * sum(g * solve(K, g)))
})

test_that("filter_decrement() does not return Inf for a flat coordinate at an edge", {
  # a loading at its edge has a diagonal of -4e-13 and no gradient: the other
  # coordinates are read, where the singular matrix gave Inf
  K <- diag(c(3, 2, -4e-13))
  g <- c(1, 1, 1e-9)
  expect_equal(fd_at(g, K, 1L), 0.5 * (1 / 3 + 1 / 2))
})

test_that("filter_decrement() reads the exploded filters as not at a maximum", {
  # a negative decrement was returned as its absolute value
  expect_identical(fd_at(exploded$k976$g, exploded$k976$K), Inf)
  # a rank-one curvature that names the term's parameters is no evidence
  expect_identical(fd_at(exploded$k997$g, exploded$k997$K), Inf)
  # above the limit the lower bound is the answer
  expect_gt(fd_at(exploded$k977$g, exploded$k977$K), mode_error_limit())
  expect_true(is.finite(fd_at(exploded$k977$g, exploded$k977$K)))
  # a curvature with a direction of ascent
  K <- diag(c(1e10, 5, 1, 2))
  K[3, 4] <- K[4, 3] <- 3
  expect_identical(fd_at(rep(1, 4), K), Inf)
})

# the data of the fits below: base datasets
dax <-data.frame(day = 1:1859,
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
  # k357 (dax_t_nu_dh) left the list in 0.218.0: with the second start of a
  # joint filter scaled equation by equation it reaches a stationary point
  # 5.58 below the maximum, nu's persistence at the edge of its chart
  # (decrement 1.3e-7, smallest eigenvalue 1e-6), rather than stopping
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

test_that("the certificate does not call a partial reading converged", {
  skip_on_cran()
  fit <- suppressWarnings(gas_fit("erupt_weib"))
  expect_identical(statmod_certificate(fit)$state, "converged")
  # the reading that the aliasing test leaves below the limit, from a curvature
  # that names the term's parameters
  cert <- testthat::with_mocked_bindings(
    statmod_certificate(fit),
    joint_reading = function(g, K, keep, own = integer(0))
      list(decrement = NA_real_, ascent = FALSE, ratio = 0, partial = TRUE))
  expect_identical(cert$state, "unknown")
  expect_true(is.na(cert$mode_error))
  expect_match(paste(cert$reason, collapse = " "), "lower bound")
})

test_that("bfgs() and lbfgs() from the start of b4 end as not converged", {
  skip_on_cran()
  skip_if_not("stop_decrement" %in% S7::prop_names(optimizers7::bfgs()))
  base_start <- structural_zeta_start
  start <- stalled$b4[[2L]]
  spec <- statmod_spec(y ~ gas(p = 1, q = 1, time = t) |
                         phi ~ gas(p = 1, q = 1, time = t),
                       distributions7::beta1_distrib(), belt2)
  design <- testthat::with_mocked_bindings(
    statmod_design(spec),
    structural_zeta_start = function(term, target = NULL) {
      z <- base_start(term, target = target)
      z[names(start)] <- start
      z
    })
  hyper <- statmod_hyper_start(spec, design)
  obj <- statmod_objective(spec, hyper, design, expected = FALSE)
  beta <- statmod_start(spec, design, obj, NULL)
  jp <- statmod_joint_pieces(spec, design, obj, hyper,
                             kinds = c("filter", "loglik"))
  u0 <- c(beta, jp$zeta())
  # both stop at -171.5056 on crit_abs_par(), with a step of 8e-9 and a
  # decrease of 0.02 still predicted at the new point
  for (opt in list(optimizers7::bfgs(), optimizers7::lbfgs())) {
    r <- optimizers7::minimize(opt, jp$fn, u0, gr = jp$gr)
    expect_false(r@converged)
    expect_match(r@message, "stop_decrement")
  }
  r <- optimizers7::minimize(optimizers7::bfgs(stop_decrement = Inf), jp$fn,
                             u0, gr = jp$gr)
  expect_true(r@converged)
})
