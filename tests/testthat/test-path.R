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
                          c(lambda = 1)), 0)
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
    fit <- statmod(y ~ lasso(x, lambda = lam),
                   distributions7::gaussian1_distrib(), dp)
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
    fit <- statmod(y ~ enet(x, lambda = 6, alpha = al),
                   distributions7::gaussian1_distrib(), dp)
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
    statmod(y ~ lasso(x, n_lambda = 12), distributions7::gaussian1_distrib(),
            dp, sparse_criterion = cv(nfolds = 5, rule = rule))
  }
  m <- f("min")
  s <- f("1se")
  expect_true(all(1:3 %in% which(abs(m@coefficients$mu[-1L]) > 1e-8)))
  expect_true(all(1:3 %in% which(abs(s@coefficients$mu[-1L]) > 1e-8)))
  expect_gte(s@hyper$mu[["lasso(x, n_lambda = 12)"]][["lambda"]],
             m@hyper$mu[["lasso(x, n_lambda = 12)"]][["lambda"]])
  # the criterion is a deviance per observation, so it is of the order the
  # gaussian's own is: 2*log(2*pi*sigma^2)/2 + 1 near 2.9 here
  expect_true(m@criterion > 2 && m@criterion < 4)
  expect_lte(m@criterion, s@criterion)
})

test_that("the folds are the caller's when the caller gives them", {
  skip_on_cran()
  fo <- rep_len(1:4, np)
  a <- statmod(y ~ lasso(x, n_lambda = 8), distributions7::gaussian1_distrib(),
               dp, sparse_criterion = cv(folds = fo))
  b <- statmod(y ~ lasso(x, n_lambda = 8), distributions7::gaussian1_distrib(),
               dp, sparse_criterion = cv(folds = fo))
  expect_equal(a@hyper$mu[["lasso(x, n_lambda = 8)"]][["lambda"]],
               b@hyper$mu[["lasso(x, n_lambda = 8)"]][["lambda"]])
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
  expect_warning(statmod(y ~ lasso(x, n_lambda = 2, min_ratio = 0.9),
                         distributions7::gaussian1_distrib(), dp,
                         sparse_criterion = aic()),
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
  # what a PATH does is not the criterion's: the same criterion is put to the
  # smooth hyperparameters of a model, which are read at the mode rather than
  # swept, so it carries neither the length of a grid nor its depth
  expect_error(cv(n_values = 12), "unused argument")
  expect_error(cv(min_ratio = 0.1), "unused argument")
  expect_error(bic(search = "cyclic"), "unused argument")
  expect_false("n_values" %in% names(S7::props(bic())))
  expect_false("min_ratio" %in% names(S7::props(bic())))
})

test_that("a bounded hyperparameter is estimated unless the term holds it", {
  # WHICH hyperparameters are estimated is the term's answer. The default is
  # to estimate every one of them, and the elastic net's alpha is bounded,
  # which is a fact about how it is SWEPT and not about whether it is.
  set.seed(51)
  n2 <- 150
  Z <- matrix(stats::rnorm(n2 * 8), n2, 8)
  colnames(Z) <- paste0("z", 1:8)
  de <- data.frame(y = as.numeric(Z %*% c(2, -1.5, 1, rep(0, 5))) +
                     stats::rnorm(n2, sd = 0.5))
  de$Z <- Z

  free <- statmod(y ~ enet(Z), distributions7::gaussian1_distrib(), de,
                  sparse_criterion = bic())
  tb <- function(f) {
    s <- summary(f)
    k <- vapply(s@tables$mu, `[[`, character(1), "kind")
    s@tables$mu[[which(k == "selection")]]$table
  }
  a <- tb(free)
  expect_identical(a$role[a$name == "alpha"], "estimated")
  expect_identical(a$source[a$name == "alpha"], "bic")
  expect_identical(a$role[a$name == "lambda"], "estimated")

  # and the term holding it is what makes it held, whatever criterion runs
  held <- statmod(y ~ enet(Z, alpha = 0.5),
                  distributions7::gaussian1_distrib(), de,
                  sparse_criterion = bic())
  b <- tb(held)
  expect_identical(b$role[b$name == "alpha"], "fixed")
  expect_identical(b$source[b$name == "alpha"], "fixed")
  expect_equal(b$estimate[b$name == "alpha"], 0.5)
  # while its lambda is still chosen
  expect_identical(b$role[b$name == "lambda"], "estimated")
  expect_output(print(summary(held)), "[fixed]", fixed = TRUE)
})

test_that("a bounded hyperparameter is swept over its own interval", {
  # it used to be accepted and ignored: the sweep walks the SIZE OF THE KINK,
  # which for the elastic net is lambda*alpha, and no admissible alpha empties
  # the block, so every point of that path was dropped and alpha came back at
  # its default marked as though a criterion had chosen it
  set.seed(52)
  n2 <- 150
  Z <- matrix(stats::rnorm(n2 * 8), n2, 8)
  colnames(Z) <- paste0("z", 1:8)
  de <- data.frame(y = as.numeric(Z %*% c(2, -1.5, 1, rep(0, 5))) +
                     stats::rnorm(n2, sd = 0.5))
  de$Z <- Z

  held <- statmod(y ~ enet(Z, alpha = 0.5),
                  distributions7::gaussian1_distrib(), de,
                  sparse_criterion = bic())
  both <- statmod(y ~ enet(Z), distributions7::gaussian1_distrib(), de,
                  sparse_criterion = bic())
  a <- both@hyper$mu[["enet(Z)"]][["alpha"]]
  expect_false(isTRUE(all.equal(a, 0.5)))
  # strictly inside its own bounds: at alpha = 0 there is no kink at all
  expect_gt(a, 0)
  expect_lt(a, 1)
  # a sweep over more values cannot end above the criterion of a sweep over
  # fewer, the held setting being one of the points it visits
  expect_lte(both@criterion, held@criterion)
  # under the product alpha is the outer axis, so it names the combination
  # each point belongs to rather than the axis the path descends
  expect_true(any(grepl("alpha=", both@history$outer$setting, fixed = TRUE)))
  expect_identical(unique(both@history$outer$name), "lambda")
})

test_that("the grid over a bounded hyperparameter excludes its endpoints", {
  pen <- penalties7::elasticnet_penalty(4)
  v <- path_grid(pen, list(lambda = 1, alpha = 0.5), "alpha", 5L)
  expect_length(v, 5L)
  expect_true(all(v > 0 & v < 1))
  expect_true(!is.unsorted(v))
  expect_true(path_bounded(pen, "alpha"))
  expect_false(path_bounded(pen, "lambda"))
  # an unbounded one that SCALES the kink is swept by kink size instead, and
  # the geometric grid is what serves a shape that does neither
  expect_true(path_by_kink(pen, list(lambda = 1, alpha = 0.5), "lambda"))
  expect_false(path_by_kink(pen, list(lambda = 1, alpha = 0.5), "alpha"))
  mc <- penalties7::mcp_penalty(4)
  expect_false(path_by_kink(mc, list(lambda = 1, gamma = 3), "gamma"))
  g <- path_grid(mc, list(lambda = 1, gamma = 3), "gamma", 6L)
  expect_length(g, 6L)
  expect_true(all(g > 1))
})

test_that("the path visits as many values as the term asked for", {
  set.seed(71)
  n2 <- 120
  Z <- matrix(stats::rbinom(n2 * 6, 1, 0.3), n2, 6)
  colnames(Z) <- paste0("z", 1:6)
  dg <- data.frame(y = as.numeric(Z %*% c(0, 0, 1.5, -1, 2, 0)) +
                     stats::rnorm(n2))
  dg$Z <- Z
  n_of <- function(f) nrow(f@history$outer)

  for (k in c(5L, 12L)) {
    f <- suppressWarnings(statmod(y ~ lasso(Z, n_lambda = k),
                                  distributions7::gaussian1_distrib(), dg,
                                  sparse_criterion = bic()))
    expect_identical(n_of(f), k)
  }
  # and the term's own default is a number on its own signature, not a
  # criterion's that nothing prints
  d <- suppressWarnings(statmod(y ~ lasso(Z),
                                distributions7::gaussian1_distrib(), dg,
                                sparse_criterion = bic()))
  expect_identical(n_of(d), as.integer(eval(formals(lasso)$n_lambda)))

  # and it is PER HYPERPARAMETER. The product visits every combination, so
  # the elastic net costs exactly n_lambda * n_alpha points
  e <- suppressWarnings(statmod(y ~ enet(Z, n_lambda = 8, n_alpha = 4),
                                distributions7::gaussian1_distrib(), dg,
                                sparse_criterion = bic()))
  expect_identical(nrow(e@history$outer), 8L * 4L)
  # one row per point: `name` and `value` carry the axis the path descends
  # and `setting` the rest of the combination
  expect_identical(unique(e@history$outer$name), "lambda")
  expect_length(unique(e@history$outer$setting), 4L)

  # the cyclic sweep costs the sum instead, once per pass, and the TERM is
  # what asks for it: the criterion is put to the smooth hyperparameters of
  # the model as well, which are read at the mode rather than swept
  cy <- suppressWarnings(statmod(
    y ~ enet(Z, n_lambda = 8, n_alpha = 4, search = "cyclic"),
    distributions7::gaussian1_distrib(), dg, sparse_criterion = bic()))
  tb <- table(cy@history$outer$name)
  expect_identical(unname(tb[["lambda"]] %% 8L), 0L)
  expect_identical(unname(tb[["alpha"]] %% 4L), 0L)
  expect_lt(nrow(cy@history$outer), 8L * 4L)
})

test_that("the path reaches as far down as the term asked", {
  set.seed(72)
  n2 <- 120
  Z <- matrix(stats::rbinom(n2 * 6, 1, 0.3), n2, 6)
  colnames(Z) <- paste0("z", 1:6)
  dm <- data.frame(y = as.numeric(Z %*% c(0, 0, 1.5, -1, 2, 0)) +
                     stats::rnorm(n2))
  dm$Z <- Z
  ratio <- function(f) {
    v <- f@history$outer$value
    min(v) / max(v)
  }
  for (mr in c(1e-1, 1e-6)) {
    f <- suppressWarnings(statmod(y ~ lasso(Z, n_lambda = 8, min_ratio = mr),
                                  distributions7::gaussian1_distrib(), dm,
                                  sparse_criterion = bic()))
    expect_equal(ratio(f), mr, tolerance = 1e-8)
  }
  # ... and the term's own default where it named none
  d <- suppressWarnings(statmod(y ~ lasso(Z, n_lambda = 8),
                                distributions7::gaussian1_distrib(), dm,
                                sparse_criterion = bic()))
  expect_equal(ratio(d), eval(formals(lasso)$min_ratio), tolerance = 1e-8)
})

test_that("the product visits every combination and each row has its own top", {
  # THE GRID IS NOT A RECTANGLE for the elastic net: the kink is lambda*alpha,
  # so the value emptying the block is kink/alpha and every alpha carries its
  # own lambda axis. The kink itself is a property of the data and does not
  # move, which is what the product of the two must reproduce.
  set.seed(81)
  n2 <- 130
  Z <- matrix(stats::rnorm(n2 * 8), n2, 8)
  colnames(Z) <- paste0("z", 1:8)
  dd <- data.frame(y = as.numeric(Z %*% c(1.6, -1.2, 0.9, rep(0, 5))) +
                     stats::rnorm(n2))
  dd$Z <- Z

  f <- suppressWarnings(statmod(y ~ enet(Z, n_lambda = 6, n_alpha = 3),
                                distributions7::gaussian1_distrib(), dd,
                                sparse_criterion = bic()))
  h <- f@history$outer
  expect_identical(nrow(h), 6L * 3L)
  tops <- tapply(h$value, h$setting, max)
  al <- as.numeric(sub("alpha=", "", names(tops)))
  expect_length(tops, 3L)
  # lambda_max * alpha is the kink that empties the block, one number
  expect_equal(as.numeric(tops * al), rep(as.numeric(tops * al)[[1L]], 3L),
               tolerance = 1e-6)
  # and the three lambda_max really are different, so the test is not
  # satisfied by a rectangle
  expect_gt(max(tops) / min(tops), 1.5)
})

test_that("scad and mcp keep ONE lambda_max whatever the shape", {
  # the controproof against writing the elastic net's relation for them: the
  # shape does NOT scale the kink -- the half-width of the subdifferential at
  # zero is lambda and nothing else -- so lambda_max cannot depend on it
  set.seed(82)
  n2 <- 130
  Z <- matrix(stats::rnorm(n2 * 8), n2, 8)
  colnames(Z) <- paste0("z", 1:8)
  dd <- data.frame(y = as.numeric(Z %*% c(1.6, -1.2, 0.9, rep(0, 5))) +
                     stats::rnorm(n2))
  dd$Z <- Z

  for (call in c("scad(Z, n_lambda = 5, n_a = 3)",
                 "mcp(Z, n_lambda = 5, n_gamma = 3)")) {
    f <- suppressWarnings(statmod(
      stats::as.formula(paste("y ~", call)),
      distributions7::gaussian1_distrib(), dd, sparse_criterion = bic()))
    h <- f@history$outer
    expect_identical(nrow(h), 5L * 3L)
    tops <- tapply(h$value, h$setting, max)
    expect_length(tops, 3L)
    expect_equal(as.numeric(tops), rep(as.numeric(tops)[[1L]], 3L))
  }
})

test_that("a written-out grid is visited as it stands", {
  set.seed(83)
  n2 <- 130
  Z <- matrix(stats::rnorm(n2 * 8), n2, 8)
  colnames(Z) <- paste0("z", 1:8)
  dd <- data.frame(y = as.numeric(Z %*% c(1.6, -1.2, 0.9, rep(0, 5))) +
                     stats::rnorm(n2))
  dd$Z <- Z

  v <- c(0.03, 0.11, 0.4, 1.5, 5.5)
  f <- suppressWarnings(statmod(y ~ lasso(Z, lambda = v),
                                distributions7::gaussian1_distrib(), dd,
                                sparse_criterion = bic()))
  h <- f@history$outer
  expect_equal(sort(h$value), sort(v))
  # walked from the emptiest fit towards the fullest, which for a lasso is
  # downwards -- the order the warm starts need
  expect_true(!is.unsorted(rev(h$value)))
  expect_true(f@hyper$mu[[1L]][["lambda"]] %in% v)
  # ESTIMATED, not held: the caller fixed where to look and not the answer
  s <- summary(f)
  k <- vapply(s@tables$mu, `[[`, character(1), "kind")
  tb <- s@tables$mu[[which(k == "selection")]]$table
  expect_identical(tb$role[tb$name == "lambda"], "estimated")

  # per hyperparameter, and NOT rescaled by the other: the elastic net's
  # lambda_max would divide by alpha, and a written-out grid is not built
  f2 <- suppressWarnings(statmod(y ~ enet(Z, lambda = v, n_alpha = 3),
                                 distributions7::gaussian1_distrib(), dd,
                                 sparse_criterion = bic()))
  h2 <- f2@history$outer
  expect_identical(nrow(h2), length(v) * 3L)
  expect_equal(sort(unique(h2$value)), sort(v))
})

test_that("the shape's floor comes from the step, not from the constant", {
  # SCAD's proximal operator needs t < a - 1, and t = 1/sum(w x^2) under a
  # diagonal map becomes t d^2 -- so which shapes the block can be FITTED at
  # is a property of the data. A grid starting just above the penalty's own
  # bound of 2 names shapes no fit could reach.
  pen <- penalties7::scad_penalty(n_coef = 5L)
  th <- list(lambda = 1, a = 3)
  # short steps: the penalty's own bound binds and nothing is raised
  expect_equal(shape_floor(pen, th, "a", rep(0.005, 5)), 2)
  # long steps: the condition binds and the floor is 1 + t
  expect_equal(shape_floor(pen, th, "a", rep(3, 5)), 4)
  expect_equal(shape_floor(pen, th, "a", rep(10, 5)), 11)
  # and it is asked of the PENALTY, so MCP's own condition, t < gamma, comes
  # out without either constant being written here
  mc <- penalties7::mcp_penalty(n_coef = 5L)
  expect_equal(shape_floor(mc, list(lambda = 1, gamma = 3), "gamma",
                           rep(3, 5)), 3)
  expect_equal(shape_floor(mc, list(lambda = 1, gamma = 3), "gamma",
                           rep(0.005, 5)), 1)
  # with no steps to read there is nothing to raise it above
  expect_equal(shape_floor(pen, th, "a", NULL), 2)

  # the grid starts above whatever the floor is, and every point of it
  # admits the operator -- which the old grid, pinned at the constant, did
  # not: 2.25 with a step of 3 is a shape no fit could use
  g <- path_grid(pen, th, "a", 6L, rep(3, 5))
  expect_true(all(g > 4))
  expect_true(all(vapply(g, function(a) !is.null(penalties7::penalty_prox_spec(
    pen, list(lambda = 1, a = a), rep(3, 5))), logical(1))))
})

test_that("the size of the kink is inverted in closed form", {
  # measured, the size is exactly a power of each hyperparameter -- one for
  # the lasso and for the elastic net in both of its own, minus one for a
  # Laplace prior written by its scale -- so the value giving a target size
  # is closed and uniroot is not run
  lp <- penalties7::distrib_penalty(
    distributions7::fixed(distributions7::laplace_distrib(), mu = 0),
    n_coef = 5L, kinks = 0)
  expect_equal(kink_power(penalties7::lasso_penalty(n_coef = 4L),
                          list(lambda = 3), "lambda")$k, 1)
  expect_equal(kink_power(penalties7::elasticnet_penalty(n_coef = 4L),
                          list(lambda = 3, alpha = 0.5), "alpha")$k, 1)
  expect_equal(kink_power(lp, list(sigma = 0.5), "sigma")$k, -1,
               tolerance = 1e-6)
  # a hyperparameter the kink does not depend on has no exponent to read
  expect_null(kink_power(penalties7::scad_penalty(n_coef = 4L),
                         list(lambda = 3, a = 3.7), "a"))

  # and the closed answer is the one the search would have found
  pen <- penalties7::elasticnet_penalty(n_coef = 4L)
  th <- list(lambda = 3, alpha = 0.5)
  for (target in c(0.25, 2, 17)) {
    expect_equal(kink_solve(pen, th, "lambda", target), target / 0.5,
                 tolerance = 1e-10)
    expect_equal(kink_scale(pen, utils::modifyList(
      th, list(lambda = kink_solve(pen, th, "lambda", target)))), target,
      tolerance = 1e-10)
  }
  # the whole grid through one exponent agrees with solving each target
  v <- path_values(pen, th, "lambda", 20, 7L, 1e-3)
  one <- vapply(exp(seq(log(20), log(0.02), length.out = 7L)),
                function(t) kink_solve(pen, th, "lambda", t), numeric(1))
  expect_equal(v, one, tolerance = 1e-10)
})

test_that("a fold carries a term's matrix input onto its own rows", {
  # data.frame(X = X, y = y) SPLITS the matrix into X.x1 ... X.xp, leaving no
  # column X, so lasso(X) reaches past the data to the matrix in the calling
  # environment: interpret_formula evaluates the call as eval(call, data, env)
  # and the name is looked up in data first, in env after. The fit is right --
  # the matrix is captured once -- but the fold could not rebuild, the name
  # still resolving to all the rows.
  skip_on_cran()
  set.seed(91)
  n2 <- 100L
  Z <- matrix(stats::rbinom(n2 * 8, 1, 0.25), n2, 8)
  colnames(Z) <- paste0("z", 1:8)
  yy <- as.numeric(Z %*% c(0, 0, 1.5, -1, 2, 0, 0, 0)) + stats::rnorm(n2)

  split <- data.frame(Z = Z, y = yy)          # no column 'Z'
  col <- data.frame(y = yy); col$Z <- Z       # the documented spelling
  expect_false("Z" %in% names(split))
  expect_true("Z" %in% names(col))

  fo <- rep_len(1:5, n2)                      # the SAME folds on both sides
  f <- function(dd) {
    suppressWarnings(statmod(y ~ 1 + lasso(Z, n_lambda = 8),
                             distributions7::gaussian1_distrib(), dd,
                             sparse_criterion = cv(folds = fo)))
  }
  a <- f(split)
  b <- f(col)
  # the two spellings are the same model, so the fold must make them the same
  # answer: not merely both finite
  expect_equal(a@hyper$mu[[1L]][["lambda"]], b@hyper$mu[[1L]][["lambda"]])
  expect_equal(a@criterion, b@criterion)
  expect_equal(unname(unlist(a@coefficients)), unname(unlist(b@coefficients)))
})

test_that("carrying a matrix onto a fold keeps it in its own storage", {
  # a sparse input is passed to avoid the memory a dense one costs, so the
  # column the fold is given must not be the densification of it
  set.seed(92)
  n2 <- 60L
  Z <- Matrix::Matrix(matrix(stats::rbinom(n2 * 6, 1, 0.2), n2, 6),
                      sparse = TRUE)
  colnames(Z) <- paste0("z", 1:6)
  dd <- data.frame(y = stats::rnorm(n2))
  spec <- statmod_spec(y ~ lasso(Z, n_lambda = 4),
                       distributions7::gaussian1_distrib(), dd, NULL, NULL)
  keep <- rep(c(TRUE, FALSE), each = n2 / 2)
  sub <- cv_bind_inputs(spec, dd[keep, , drop = FALSE], keep, n2)
  expect_true("Z" %in% names(sub))
  expect_true(methods::is(sub[["Z"]], "Matrix"))
  expect_identical(nrow(sub[["Z"]]), sum(keep))

  # a term whose input the data already carries is left alone, and so is a
  # FORMULA input, which keeps being rebuilt on the fold's own rows
  df <- data.frame(y = stats::rnorm(n2), z = stats::rnorm(n2))
  sp2 <- statmod_spec(y ~ lasso(~z, n_lambda = 4),
                      distributions7::gaussian1_distrib(), df, NULL, NULL)
  expect_identical(names(cv_bind_inputs(sp2, df[keep, , drop = FALSE], keep,
                                        n2)),
                   names(df))
})
