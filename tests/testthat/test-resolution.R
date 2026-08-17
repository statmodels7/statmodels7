test_that("the criterion's resolution is read at the fit and tracks the inner rule", {
  # It is what a difference in the criterion has to exceed to mean anything:
  # the criterion is read at the penalized mode, the inner fit stops short of
  # that mode by whatever its rule allows, and the criterion read at the
  # displaced mode differs by exactly what an evaluation from another warm
  # start would differ by. A LOOSER inner rule leaves the mode less well
  # located, so the resolution has to grow with it.
  skip_on_cran()
  set.seed(3)
  n <- 300L
  x <- sort(runif(n, -3, 3))
  d <- data.frame(x = x, y = sin(1.5 * x) + rnorm(n, sd = 0.35))
  meth <- reml(hessian = "observed")

  at_tol <- function(tol) {
    inner <- iwls(tol = tol, maxit = 300L)
    spec <- statmod_spec(y ~ s(x, k = 10),
                         distributions7::gaussian1_distrib(), d, NULL, NULL,
                         linpar = linpar_options())
    design <- statmod_design(spec)
    hyper <- statmod_hyper_start(spec, design)
    blocks <- statmod_blocks(spec, design)
    cfg <- inner_settings(inner)
    obj <- statmod_objective(spec, hyper, design, cfg$expected, cfg$approx)
    beta0 <- statmod_start(spec, design, obj, NULL)
    basis <- integrated_basis(spec, design, meth@kind)
    bud <- method_budget(inner)
    res <- statmod_alternate(spec, design, blocks, hyper, inner, beta0,
                             cfg$expected, cfg$approx, bud$maxit, bud$tol,
                             list(inner = FALSE, blocks = FALSE,
                                  outer = FALSE, path = FALSE))
    cf <- res$obj$split(res$par)
    ctx <- outer_context(spec, design, cf, hyper, cfg$approx)
    m <- statmod_marginal(spec, design, cf, hyper, meth, cfg$approx, basis, ctx)
    st <- list(ok = TRUE, score = res$obj$gr(res$par), split = res$obj$split,
               par = res$par, cf = cf, hy = hyper, ctx = ctx, value = m$value)
    crit_at <- function(cf, hy, par, ctx = NULL) {
      statmod_marginal(spec, design, cf, hy, meth, cfg$approx, basis, ctx)
    }
    list(res = criterion_resolution(st, spec, design, meth, crit_at), st = st,
         crit_at = crit_at, spec = spec, design = design)
  }

  loose <- at_tol(1e-4)
  tight <- at_tol(1e-8)
  expect_true(is.finite(loose$res) && loose$res > 0)
  expect_true(is.finite(tight$res) && tight$res > 0)
  # the structural property: a mode located ten thousand times more loosely
  # leaves a criterion that cannot be read as finely
  expect_gt(loose$res, tight$res)

  # a state with nothing to read from answers NA rather than a number
  bad <- loose$st
  bad$score <- NULL
  expect_true(is.na(criterion_resolution(bad, loose$spec, loose$design, meth,
                                         loose$crit_at)))
  bad <- loose$st
  bad$ok <- FALSE
  expect_true(is.na(criterion_resolution(bad, loose$spec, loose$design, meth,
                                         loose$crit_at)))
})


test_that("the resolution is read off the criterion the SEARCH runs", {
  # ⚠️ The defect this pins: reading the marginal criterion of a fit whose
  # search is aic() answers for a quantity that search never sees. The two
  # differ by orders -- one is a Laplace approximation of order 1e2 here and
  # the other a penalized log-likelihood -- and the number that came back
  # stopped two aic() fits short of their own optimum. The criterion is
  # therefore passed in as a function rather than chosen inside.
  skip_on_cran()
  set.seed(3)
  n <- 300L
  x <- sort(runif(n, -3, 3))
  d <- data.frame(x = x, y = sin(1.5 * x) + rnorm(n, sd = 0.35))

  # the two criteria at one point are not the same number, which is what makes
  # reading the wrong one an error rather than an inaccuracy
  a <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), d,
               outer_criterion = aic())
  r <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), d,
               outer_criterion = reml(hessian = "observed"))
  expect_gt(abs(a@criterion - r@criterion), 1)

  # and an aic() fit reaches its own optimum: the criterion at the reported
  # hyperparameter beats the criterion just either side of it
  spec <- a@spec
  design <- statmod_design(spec)
  nm <- "s(x, k = 10)"
  at <- function(v) {
    f <- statmod(y ~ s(x, k = 10, lambda = v),
                 distributions7::gaussian1_distrib(), d)
    statmod_pe(f@spec, statmod_design(f@spec), f@coefficients, f@hyper,
               aic())$value
  }
  lam <- a@hyper$mu[[nm]][["lambda"]]
  expect_lt(a@criterion, at(lam * 1.3))
  expect_lt(a@criterion, at(lam / 1.3))
})
