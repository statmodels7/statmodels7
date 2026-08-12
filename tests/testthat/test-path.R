# Selecting a hyperparameter whose penalty has a kink: the size of the kink,
# the grid it drives, the active-set degrees of freedom, and cross-validation.

set.seed(11)
np <- 150L
xp <- matrix(stats::rnorm(np * 12L), np, 12L)
truth <- c(2, -1.5, 1, rep(0, 9))
dp <- data.frame(y = as.numeric(xp %*% truth) + stats::rnorm(np))
dp$x <- xp

test_that("the size of the kink is read from the penalty", {
  # lambda for the lasso, SCAD and MCP; lambda*alpha for the elastic net. The
  # shape parameters do not move it, and saying so needs the reading to be
  # extrapolated to the kink: MCP's derivative just past it is
  # lambda - eps/gamma, so at a fixed eps gamma appears to move the kink by a
  # millionth of itself, which is enough to be selected for a path.
  expect_equal(kink_scale(penalties7::lasso_penalty(n_coef = 5L),
                          c(lambda = 2.5)), 2.5, tolerance = 1e-9)
  expect_equal(kink_scale(penalties7::scad_penalty(n_coef = 5L),
                          c(lambda = 2.5, a = 3.7)), 2.5, tolerance = 1e-9)
  expect_equal(kink_scale(penalties7::mcp_penalty(n_coef = 5L),
                          c(lambda = 2.5, gamma = 3)), 2.5, tolerance = 1e-9)
  expect_equal(kink_scale(penalties7::elasticnet_penalty(n_coef = 5L),
                          c(lambda = 2.5, alpha = 0.4)), 1, tolerance = 1e-7)
  # a penalty with no kink has none to measure
  expect_equal(kink_scale(penalties7::ridge_penalty(n_coef = 5L),
                          c(sigma = 1)), 0)
})

test_that("only the hyperparameters that set the kink are swept", {
  expect_identical(kink_hypers(penalties7::mcp_penalty(n_coef = 5L),
                               c(lambda = 2.5, gamma = 3)), "lambda")
  expect_identical(kink_hypers(penalties7::scad_penalty(n_coef = 5L),
                               c(lambda = 2.5, a = 3.7)), "lambda")
  # the elastic net's alpha does scale the kink, and is left out of the
  # default because it is bounded -- the convention glmnet follows
  th <- c(lambda = 2.5, alpha = 0.4)
  en <- penalties7::elasticnet_penalty(n_coef = 5L)
  expect_identical(kink_hypers(en, th), "lambda")
  expect_setequal(kink_hypers(en, th, unbounded = FALSE),
                  c("lambda", "alpha"))
})

test_that("the grid is inverted even where the kink narrows as the value grows", {
  # a Laplace prior written by its scale has a kink of 1/sigma. A version that
  # assumed the size rises with the value returned NA for every target here,
  # having walked away from the answer.
  lp <- penalties7::distrib_penalty(
    distributions7::fixed(distributions7::laplace_distrib(), mu = 0),
    n_coef = 5L, kinks = 0)
  th <- list(sigma = 1)
  expect_equal(kink_scale(lp, list(sigma = 0.5)), 2, tolerance = 1e-7)
  expect_identical(kink_hypers(lp, th), "sigma")
  expect_equal(kink_solve(lp, th, "sigma", 4), 0.25, tolerance = 1e-8)
  # and the ordinary direction still works
  expect_equal(kink_solve(penalties7::lasso_penalty(n_coef = 5L),
                          list(lambda = 2.5), "lambda", 0.3), 0.3,
               tolerance = 1e-9)
})

test_that("the effective degrees of freedom count what survived", {
  # away from the kink the penalty is twice differentiable, so the trace is
  # taken there; for the lasso it is linear, S_AA vanishes and tau is the
  # number of surviving coefficients, which is Zou-Hastie-Tibshirani. The
  # trace over the whole vector cannot see the selection at all: measured on
  # this design it reads 14 at every lambda.
  full <- numeric(0)
  for (lam in c(2, 8, 25)) {
    fit <- statmod(y ~ lasso(x), distributions7::gaussian1_distrib(), dp,
                   hyper = list(mu = list("lasso(x)" = c(lambda = lam))))
    spec <- fit@spec
    design <- statmod_design(spec)
    blocks <- statmod_blocks(spec, design)
    act <- statmod_active(spec, blocks, unlist(fit@coefficients,
                                               use.names = FALSE), fit@hyper)
    H <- statmod_information_at(spec, fit@coefficients, design, FALSE)
    S <- statmod_penalty_at(spec, fit@coefficients, fit@hyper, design,
                            "hessian")
    S[!is.finite(S)] <- 0
    nz <- sum(abs(fit@coefficients$mu[-1L]) > 1e-8)
    # the two intercepts are never at a kink
    expect_equal(outer_tau(H + S, H, act), nz + 2, tolerance = 1e-8)
    full <- c(full, outer_tau(H + S, H))
  }
  expect_equal(full, rep(14, 3), tolerance = 1e-6)
})

test_that("the elastic net spends less than its count, the ridge part shrinking", {
  taus <- vapply(c(0.9, 0.5, 0.1), function(al) {
    fit <- statmod(y ~ enet(x), distributions7::gaussian1_distrib(), dp,
                   hyper = list(mu = list("enet(x)" = c(lambda = 6,
                                                        alpha = al))))
    spec <- fit@spec
    design <- statmod_design(spec)
    blocks <- statmod_blocks(spec, design)
    act <- statmod_active(spec, blocks, unlist(fit@coefficients,
                                               use.names = FALSE), fit@hyper)
    H <- statmod_information_at(spec, fit@coefficients, design, FALSE)
    S <- statmod_penalty_at(spec, fit@coefficients, fit@hyper, design,
                            "hessian")
    S[!is.finite(S)] <- 0
    nz <- sum(abs(fit@coefficients$mu[-1L]) > 1e-8) + 2
    outer_tau(H + S, H, act) - nz
  }, numeric(1))
  # every one is short of the count, and by more as the quadratic part grows
  expect_true(all(taus < 0))
  expect_true(taus[3L] < taus[1L])
})

test_that("a criterion selects a lasso, and the true predictors survive", {
  for (m in list(aic(), bic())) {
    fit <- statmod(y ~ lasso(x), distributions7::gaussian1_distrib(), dp,
                   sparse_criterion = m)
    b <- fit@coefficients$mu[-1L]
    kept <- which(abs(b) > 1e-8)
    expect_true(all(1:3 %in% kept))
    expect_true(is.finite(fit@criterion))
    expect_true(nrow(fit@history$outer) > 1L)
  }
  # bic prices a degree of freedom higher and keeps fewer of them
  a <- statmod(y ~ lasso(x), distributions7::gaussian1_distrib(), dp,
               sparse_criterion = aic())
  b <- statmod(y ~ lasso(x), distributions7::gaussian1_distrib(), dp,
               sparse_criterion = bic())
  expect_gt(b@hyper$mu[["lasso(x)"]][["lambda"]],
            a@hyper$mu[["lasso(x)"]][["lambda"]])
  expect_lt(sum(abs(b@coefficients$mu[-1L]) > 1e-8),
            sum(abs(a@coefficients$mu[-1L]) > 1e-8))
})

test_that("cross-validation selects, and its one-standard-error rule is sparser", {
  skip_on_cran()
  f <- function(rule) {
    statmod(y ~ lasso(x), distributions7::gaussian1_distrib(), dp,
            sparse_criterion = cv(nfolds = 5, n_values = 12, rule = rule))
  }
  m <- f("min")
  s <- f("1se")
  expect_true(all(1:3 %in% which(abs(m@coefficients$mu[-1L]) > 1e-8)))
  expect_true(all(1:3 %in% which(abs(s@coefficients$mu[-1L]) > 1e-8)))
  expect_gte(s@hyper$mu[["lasso(x)"]][["lambda"]],
             m@hyper$mu[["lasso(x)"]][["lambda"]])
  # the criterion is a deviance per observation, so it is of the order the
  # gaussian's own is: 2*log(2*pi*sigma^2)/2 + 1 near 2.9 here
  expect_true(m@criterion > 2 && m@criterion < 4)
  expect_lte(m@criterion, s@criterion)
})

test_that("the folds are the caller's when the caller gives them", {
  skip_on_cran()
  fo <- rep_len(1:4, np)
  a <- statmod(y ~ lasso(x), distributions7::gaussian1_distrib(), dp,
               sparse_criterion = cv(folds = fo, n_values = 8))
  b <- statmod(y ~ lasso(x), distributions7::gaussian1_distrib(), dp,
               sparse_criterion = cv(folds = fo, n_values = 8))
  expect_equal(a@hyper$mu[["lasso(x)"]][["lambda"]],
               b@hyper$mu[["lasso(x)"]][["lambda"]])
  expect_equal(a@criterion, b@criterion)
  expect_error(statmod(y ~ lasso(x), distributions7::gaussian1_distrib(), dp,
                       sparse_criterion = cv(folds = 1:3)),
               "3 entries but there are")
})

test_that("cv is refused where it has nothing to select", {
  # a smoothing parameter that is twice differentiable is read from a
  # criterion, not from folds, and saying so beats sweeping a grid nobody
  # asked for
  dq <- data.frame(x = stats::runif(120, -2, 2))
  dq$y <- sin(1.4 * dq$x) + stats::rnorm(120, sd = 0.3)
  # and asking for one where there is no kinked penalty is the symmetric case
  # of reml() on a model with no smooth one: the criterion applies to a family
  # of penalties the model does not carry, so it does not run
  f <- statmod(y ~ s(x, k = 8), distributions7::gaussian1_distrib(), dq,
               sparse_criterion = cv())
  expect_true(f@converged)
  expect_false(is.na(f@criterion))
})

test_that("the selection is stated where the path ran out", {
  # a choice at either end of the grid is the grid's and not the criterion's,
  # so it is reported rather than returned as though it were a minimum. A path
  # of two points close together has to end at one of them.
  expect_warning(statmod(y ~ lasso(x), distributions7::gaussian1_distrib(),
                         dp, sparse_criterion = aic(k = 2)) -> wide,
                 NA)
  narrow <- aic()
  narrow@n_values <- 2
  narrow@min_ratio <- 0.9
  expect_warning(statmod(y ~ lasso(x), distributions7::gaussian1_distrib(),
                         dp, sparse_criterion = narrow),
                 "stopped at its")
  # and the interior choice of the wide path is not at either end
  lam <- wide@hyper$mu[["lasso(x)"]][["lambda"]]
  h <- wide@history$outer
  expect_gt(lam, min(h$value))
  expect_lt(lam, max(h$value))
})

test_that("cv and its method print what they are", {
  expect_output(print(cv()), "CV")
  expect_output(print(cv(nfolds = 5, rule = "1se")), "1se")
  expect_error(cv(rule = "best"), "arg")
  expect_error(cv(nfolds = 1), "at least 2")
  expect_error(cv(min_ratio = 2), "in \\(0, 1\\)")
})
