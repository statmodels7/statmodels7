test_that("a kinked penalty on a block that moves fits the model it names", {
  # The block of an nl() term is the JACOBIAN at the current coefficients, and
  # what the term contributes is X beta + adj rather than X beta. A coordinate
  # descent reading the block as it was built, against a linear predictor that
  # is not the model's, converges to a point that is not the mode and reports
  # success: measured before the fix, a log-likelihood of -339.74 against the
  # ridge control's +155.45 and a rate of 0.22 against a truth of 0.70.
  #
  # The check is the one that found it. At a HELD lambda small enough that
  # neither shrinks, a lasso and a ridge on the same block must agree.
  skip_on_cran()
  set.seed(9)
  n <- 300L
  m <- 10L
  grp <- factor(rep(seq_len(m), each = n / m))
  x <- stats::runif(n, 0, 3)
  a_true <- rep(3, m) + stats::rnorm(m, sd = 0.4)
  d <- data.frame(x = x, grp = grp,
                  y = a_true[grp] * exp(-0.7 * x) + stats::rnorm(n, sd = 0.15))

  fit <- function(pen) {
    statmod(stats::as.formula(paste("y ~", deparse(pen))),
            distributions7::gaussian1_distrib(), d, outer_criterion = NULL)
  }
  la <- fit(quote(nl(~ a * exp(-r * x), a ~ 0 + lasso(~ grp, lambda = 0.01))))
  ri <- fit(quote(nl(~ a * exp(-r * x), a ~ 0 + ridge(~ grp, lambda = 0.1))))

  expect_equal(as.numeric(logLik(la)), as.numeric(logLik(ri)),
               tolerance = 1e-3)
  # the rate is the coefficient the wrong predictor moved furthest, and it is
  # what says the fit is of this model and not of another
  bl <- unlist(la@coefficients, use.names = FALSE)
  br <- unlist(ri@coefficients, use.names = FALSE)
  expect_equal(bl, br, tolerance = 1e-2)
  # anchored on the truth as well as on the control, so that both agreeing on
  # a wrong answer would still fail: the fitted mean has to be the curve the
  # data were simulated from
  expect_gt(stats::cor(as.numeric(predict(la, type = "response")$mu),
                       a_true[grp] * exp(-0.7 * x)), 0.99)

  # and the selection: with every true amplitude far from zero, a lasso at a
  # small lambda must keep all of them. The path fit that failed reported a
  # rate of 16.4 and is what this guards against.
  design <- statmod_design(la@spec)
  u <- statmod_penalized(la@spec, design)[[1L]]
  expect_identical(sum(abs(bl[u$index]) > 1e-8), length(u$index))
})


test_that("the KKT conditions hold at the point the descent reached", {
  # The reference that shares no arithmetic with the iteration: the fit got
  # there through a proximal table and a compiled sweep, and this differentiates
  # the penalized objective numerically instead. Away from the kink it is
  # smooth, so on coordinates off zero the gradient must vanish; on coordinates
  # held at zero the gradient must sit inside the interval the kink opens.
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  set.seed(9)
  n <- 300L
  m <- 10L
  grp <- factor(rep(seq_len(m), each = n / m))
  x <- stats::runif(n, 0, 3)
  a_true <- rep(3, m) + stats::rnorm(m, sd = 0.4)
  d <- data.frame(x = x, grp = grp,
                  y = a_true[grp] * exp(-0.7 * x) + stats::rnorm(n, sd = 0.15))

  for (lam in c(0.01, 20)) {
    f <- statmod(stats::as.formula(sprintf(
      "y ~ nl(~ a * exp(-r * x), a ~ 0 + lasso(~ grp, lambda = %g))", lam)),
      distributions7::gaussian1_distrib(), d, outer_criterion = NULL)
    design <- statmod_design(f@spec)
    cfg <- inner_settings(iwls())
    obj <- statmod_objective(f@spec, f@hyper, design, cfg$expected, cfg$approx)
    bh <- unlist(f@coefficients, use.names = FALSE)
    g <- numDeriv::grad(function(b) obj$fn(b), bh)

    u <- statmod_penalized(f@spec, design)[[1L]]
    act <- u$index[abs(bh[u$index]) > 1e-8]
    ina <- setdiff(u$index, act)
    if (length(act)) expect_lt(max(abs(g[act])), 1e-6)
    if (length(ina)) {
      s <- kink_scale(u$penalty, as.list(f@hyper[[u$param]][[u$key]]))
      expect_true(all(abs(g[ina]) <= s * (1 + 1e-6)))
    }
  }
})
