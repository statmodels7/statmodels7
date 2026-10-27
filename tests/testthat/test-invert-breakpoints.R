# Intervals by inversion on the position of a break-point: the search stops at
# the range of the variable, and a side the data do not bound is infinite.

two_break_fit <- function(flip = FALSE) {
  set.seed(1)
  dd <- data.frame(x = sort(runif(120, 0, 10)))
  dd$y <- 1 + 0.5 * dd$x + 2 * pmax(dd$x - 4, 0) - 0.25 * pmax(dd$x - 7, 0) +
    rnorm(120, sd = 0.6)
  # the mirror image of the variable swaps the two break-points and the two
  # sides of every interval
  if (flip) {
    dd$x <- 10 - dd$x
    dd <- dd[order(dd$x), ]
  }
  fit <- suppressWarnings(statmod(y ~ seg(x, npsi = 2),
                                  distributions7::gaussian1_distrib(), dd))
  list(fit = fit, data = dd)
}

tb <- two_break_fit()
tbf <- two_break_fit(flip = TRUE)
q95 <- stats::qchisq(0.95, 1)

test_that("a break-point position has the range of its variable, another coefficient has none", {
  des <- statmod_design(tb$fit@spec)
  cn <- des$mu$coef_names
  rng <- range(tb$data$x)
  for (cf in c("seg.psi1", "seg.psi2")) {
    expect_equal(break_edges(tb$fit@spec, des, "mu", match(cf, cn)), rng)
  }
  expect_equal(break_edges(tb$fit@spec, des, "mu", match("seg.beta", cn)),
               c(-Inf, Inf))
  expect_equal(break_edges(tb$fit@spec, des, "mu", match("(Intercept)", cn)),
               c(-Inf, Inf))
  expect_equal(break_edges(tb$fit@spec, des, "sigma", 1L), c(-Inf, Inf))
})

test_that("the restricted refit is the same at and beyond the end of the range", {
  o0 <- joint_objective(tb$fit)
  at <- function(b) {
    statmod_stat_at(tb$fit, "mu", "seg.psi2", b, lr_test(), bayesian_vcov(),
                    objective0 = o0)$statistic
  }
  top <- max(tb$data$x)
  expect_equal(at(top), at(top + 5), tolerance = 1e-8)
  expect_equal(at(top), at(top + 100), tolerance = 1e-8)
  expect_lt(at(top), q95)
})

test_that("a side where the statistic stays under the cutoff is infinite, in few refits", {
  calls <- 0L
  stat_at <- statmod_stat_at
  testthat::local_mocked_bindings(
    statmod_stat_at = function(...) {
      calls <<- calls + 1L
      stat_at(...)
    })
  ci <- statmod_invert(tb$fit, "mu", "seg.psi2", 0.95, lr_test())
  expect_identical(unname(ci[["upper"]]), Inf)
  expect_true(is.finite(ci[["lower"]]))
  # the search used to try 52 positions up to 1.6e8 beyond the data
  expect_lte(calls, 30L)
  # the finite end is a crossing of the cutoff: under it just inside, over it
  # just outside, and equal to it at the root
  o0 <- joint_objective(tb$fit)
  at <- function(b) {
    stat_at(tb$fit, "mu", "seg.psi2", b, lr_test(), bayesian_vcov(),
            objective0 = o0)$statistic
  }
  expect_equal(at(ci[["lower"]]), q95, tolerance = 1e-3)
  expect_lt(at(ci[["lower"]] + 0.05), q95)
  expect_gt(at(ci[["lower"]] - 0.05), q95)
})

test_that("a position bounded on both sides has two finite limits around the estimate", {
  ci <- statmod_invert(tb$fit, "mu", "seg.psi1", 0.95, lr_test())
  psi <- coef(tb$fit, readable = FALSE)$mu[["seg.psi1"]]
  expect_true(all(is.finite(ci)))
  expect_lt(ci[["lower"]], psi)
  expect_gt(ci[["upper"]], psi)
})

test_that("the mirror image of the variable mirrors the intervals", {
  a <- statmod_invert(tb$fit, "mu", "seg.psi2", 0.95, lr_test())
  b <- statmod_invert(tbf$fit, "mu", "seg.psi1", 0.95, lr_test())
  expect_identical(unname(b[["lower"]]), -Inf)
  expect_equal(unname(b[["upper"]]), 10 - unname(a[["lower"]]),
               tolerance = 0.02)
  a1 <- statmod_invert(tb$fit, "mu", "seg.psi1", 0.95, lr_test())
  b2 <- statmod_invert(tbf$fit, "mu", "seg.psi2", 0.95, lr_test())
  expect_equal(unname(b2), sort(10 - unname(a1)), tolerance = 0.02)
})

test_that("confint() reports the infinite side and the graph draws it to the range", {
  ci <- confint(tb$fit, parm = "mu:seg.psi2", test = lr_test(),
                readable = FALSE)
  expect_identical(ci$upper, Inf)
  expect_true(is.finite(ci$lower))
  term <- tb$fit@spec@terms$mu[["seg(x, npsi = 2)"]]
  ctx <- list(param = "mu", name = "seg(x, npsi = 2)", vars = "x",
              cols = NULL, data = tb$data)
  d <- plot_data(plot_breakpoint(profile = FALSE), term, tb$fit, ctx = ctx)
  expect_identical(d$breaks$upper[2L], Inf)
  expect_equal(d$breaks$to[2L], max(tb$data$x))
  expect_equal(d$breaks$from[2L], d$breaks$lower[2L])
  expect_equal(d$breaks$from[1L], d$breaks$lower[1L])
  pdf(NULL)
  on.exit(dev.off())
  p <- autoplot(tb$fit, term = "seg(x, npsi = 2)", data = tb$data,
                plot_type = plot_breakpoint(profile = FALSE))
  expect_no_error(print(p))
})
