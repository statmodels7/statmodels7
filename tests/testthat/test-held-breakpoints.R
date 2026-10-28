# A restricted refit that holds a break-point of a seg() term with two
# break-points: the free one stays on its side of the held one, at a distance
# from it, and the refit converges.

held_fit <- function() {
  set.seed(1)
  dd <- data.frame(x = sort(runif(120, 0, 10)))
  dd$y <- 1 + 0.5 * dd$x + 2 * pmax(dd$x - 4, 0) - 0.25 * pmax(dd$x - 7, 0) +
    rnorm(120, sd = 0.6)
  list(fit = suppressWarnings(statmod(y ~ seg(x, npsi = 2),
                                      distributions7::gaussian1_distrib(), dd)),
       data = dd)
}

hb <- held_fit()
# the nearest position of a free break-point beside a held one at h, with two
# observations between them: the second observation above it, or the second at
# or below it
hb_apart <- function(h) {
  xs <- sort(hb$data$x)
  c(left = sort(xs[xs <= h], decreasing = TRUE)[2L], right = xs[xs > h][2L])
}

# the objective at two positions, the other coefficients at their optimum: a
# gaussian response with its variance at the maximum likelihood estimate
hb_objective <- function(p) {
  x <- hb$data$x
  X <- cbind(1, x, sapply(p, function(q) pmax(x - q, 0)))
  n <- length(x)
  n / 2 * (log(2 * pi * sum(qr.resid(qr(X), hb$data$y)^2) / n) + 1)
}

test_that("a break-point held beyond the free one gives a converged refit with the free one kept apart", {
  calls <- 0L
  fit_smooth_orig <- fit_smooth
  testthat::local_mocked_bindings(
    fit_smooth = function(...) {
      calls <<- calls + 1L
      fit_smooth_orig(...)
    })
  for (b in c(2, 3)) {
    calls <- 0L
    r <- statmod_restrict(hb$fit, "mu", "seg.psi2", b)
    expect_true(r$converged)
    # the alternation used to run its whole budget of 1000 passes here
    expect_lte(calls, 5L)
    s <- statmod_stat_at(hb$fit, "mu", "seg.psi2", b, lr_test())
    # the statistic was 110 and 118 at these two values, off by a factor of 25
    expect_lt(s$statistic, 10)
    cf <- r$coefficients$mu
    pos <- cf[grepl("seg\\.psi", names(cf))]
    expect_equal(unname(pos[grepl("held", names(pos))]), b)
    # two observations at least lie between the free and the held one
    free <- unname(pos[!grepl("held", names(pos))])
    expect_true(free >= hb_apart(b)[["right"]] - 1e-8 ||
                  free <= hb_apart(b)[["left"]] + 1e-8)
  }
})

test_that("the coefficients of the refit are listed in the order of the positions, the held one marked", {
  r <- statmod_restrict(hb$fit, "mu", "seg.psi2", 3)
  cf <- r$coefficients$mu
  nm <- names(cf)
  expect_identical(sub(" [held]", "", nm, fixed = TRUE),
                   c("(Intercept)", "seg.beta", "seg.gamma1", "seg.gamma2",
                     "seg.psi1", "seg.psi2"))
  expect_equal(sum(grepl("held", nm)), 2L)
  expect_identical(grepl("held", nm[3:4]), grepl("held", nm[5:6]))
  # the positions are in order
  expect_false(is.unsorted(unname(cf[5:6])))
  # a break-point held above the free one keeps its own name
  r <- statmod_restrict(hb$fit, "mu", "seg.psi2", 8)
  expect_identical(names(r$coefficients$mu)[5:6],
                   c("seg.psi1", "seg.psi2 [held]"))
  expect_lte(r$coefficients$mu[["seg.psi1"]], hb_apart(8)[["left"]] + 1e-8)
  # a model without a term of this kind is returned as it was
  fit1 <- suppressWarnings(statmod(y ~ x, distributions7::gaussian1_distrib(),
                                   hb$data))
  expect_null(names(statmod_restrict(fit1, "mu", "x", 0.4)$coefficients$mu))
})

test_that("a free break-point that starts too close to the held one is moved out", {
  b <- coef(hb$fit, readable = FALSE)$mu[["seg.psi1"]] + 0.05
  r <- statmod_restrict(hb$fit, "mu", "seg.psi2", b)
  expect_true(r$converged)
  cf <- r$coefficients$mu
  pos <- cf[grepl("seg\\.psi", names(cf))]
  expect_equal(sum(grepl("held", names(pos))), 1L)
  free <- unname(pos[!grepl("held", names(pos))])
  expect_true(free >= hb_apart(b)[["right"]] - 1e-8 ||
                free <= hb_apart(b)[["left"]] + 1e-8)
})

test_that("the restricted objective is the best one on a grid of the allowed positions of the free break-point", {
  # psi1 held at 4.5: the free break-point may be on either side of it
  lim <- as.numeric(stats::quantile(hb$data$x, c(0.05, 0.95), names = FALSE))
  for (b in c(4.5, 8)) {
    r <- statmod_restrict(hb$fit, "mu", "seg.psi1", b)
    a <- hb_apart(b)
    grid <- c(seq(lim[1], a[["left"]], length.out = 400L),
              seq(a[["right"]], lim[2], length.out = 400L))
    ref <- min(vapply(grid, function(f) hb_objective(c(b, f)), numeric(1)))
    expect_equal(r$objective, ref, tolerance = 1e-4)
  }
})

test_that("the neighbors of a break-point are the fitted positions next to it", {
  est <- coef(hb$fit, readable = FALSE)$mu[c("seg.psi1", "seg.psi2")]
  spec <- hb$fit@spec
  des <- statmod_design(spec)
  at <- function(cf) match(cf, des$mu$coef_names)
  nb <- function(cf) neighbor_positions(hb$fit, spec, des, "mu", at(cf))
  expect_equal(nb("seg.psi1"), c(NA, unname(est[2])))
  expect_equal(nb("seg.psi2"), c(unname(est[1]), NA))
  expect_equal(nb("seg.gamma1"), c(NA_real_, NA_real_))
  expect_equal(nb("(Intercept)"), c(NA_real_, NA_real_))
})

test_that("the search for the limits of a break-point stays between the fitted positions of the neighbors", {
  q95 <- stats::qchisq(0.95, 1)
  est <- coef(hb$fit, readable = FALSE)$mu[c("seg.psi1", "seg.psi2")]
  ci2 <- statmod_invert(hb$fit, "mu", "seg.psi2", 0.95, lr_test())
  expect_gt(ci2[["lower"]], est[["seg.psi1"]])
  expect_equal(unname(ci2[["lower"]]), 4.2164, tolerance = 1e-3)
  expect_identical(unname(ci2[["upper"]]), Inf)
  # the limit is a crossing of the cutoff
  s <- statmod_stat_at(hb$fit, "mu", "seg.psi2", ci2[["lower"]], lr_test())
  expect_equal(s$statistic, q95, tolerance = 1e-3)
  ci1 <- statmod_invert(hb$fit, "mu", "seg.psi1", 0.95, lr_test())
  expect_lt(ci1[["upper"]], est[["seg.psi2"]])
  expect_equal(unname(ci1), c(3.618, 4.0195), tolerance = 1e-3)
})

test_that("a request for break-points that are all aliased gives NA limits and a warning", {
  set.seed(8)
  dd <- data.frame(x = sort(runif(100, 0, 10)))
  dd$y <- 1 + 0.5 * dd$x + 2 * pmax(dd$x - 5, 0) + rnorm(100, sd = 0.5)
  fit <- suppressWarnings(statmod(y ~ seg(x, npsi = 2),
                                  distributions7::gaussian1_distrib(), dd))
  skip_if_not("mu:seg.psi1" %in% fit@aliased,
              "the break-point is not aliased on this platform")
  expect_warning(
    ci <- confint(fit, parm = "mu:seg.psi1", test = lr_test(), readable = FALSE),
    "nothing to invert")
  expect_identical(rownames(ci), "mu:seg.psi1")
  expect_true(all(is.na(c(ci$estimate, ci$lower, ci$upper))))
  # the same row among testable ones keeps NA and gives no warning
  expect_no_warning(
    cj <- confint(fit, parm = c("mu:seg.psi1", "mu:seg.beta"), test = lr_test(),
                  readable = FALSE))
  expect_true(all(is.na(c(cj["mu:seg.psi1", "lower"], cj["mu:seg.psi1", "upper"]))))
  expect_true(all(is.finite(c(cj["mu:seg.beta", "lower"], cj["mu:seg.beta", "upper"]))))
})
