g1 <- distributions7::gaussian1_distrib()

quiet <- function(e) suppressWarnings(suppressMessages(e))

het_fit <- function(n = 150, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = runif(n))
  d$y <- 1 + 2 * d$x + rnorm(n, sd = exp(3 * d$x) / 5)
  list(data = d, fit = statmod(y ~ x, g1, d))
}

test_that("the case bootstrap converges to the sandwich and the parametric to the model", {
  hf <- het_fit()
  fit <- hf$fit
  d <- hf$data
  set.seed(11)
  fit <- add_bootstrap(fit, n = 200, resample = case_resample())
  fit <- add_bootstrap(fit, n = 200, resample = parametric_resample())
  se_case <- sqrt(diag(vcov(fit, vcov_type = bootstrap_vcov(
    resample = case_resample()))))
  se_par <- sqrt(diag(vcov(fit, vcov_type = bootstrap_vcov(
    resample = parametric_resample()))))
  # HC0 written out for the mean's coefficients
  X <- cbind(1, d$x)
  r <- as.numeric(d$y - X %*% fit@coefficients$mu)
  A <- solve(crossprod(X))
  hc0 <- sqrt(diag(A %*% crossprod(X * r) %*% A))
  tol <- 4 / sqrt(2 * 200)
  expect_lt(abs(se_case[["mu:x"]] / hc0[[2]] - 1), tol)
  expect_lt(abs(se_case[["mu:(Intercept)"]] / hc0[[1]] - 1), tol)
  se_f <- sqrt(diag(vcov(fit, vcov_type = frequentist_vcov())))
  expect_lt(abs(se_par[["mu:x"]] / se_f[["mu:x"]] - 1), tol)
  # the two schemes answer different questions: the model's standard error of
  # the slope is well below the sandwich here
  expect_gt(hc0[[2]] / se_f[["mu:x"]], 1.25)
})

test_that("the parametric bootstrap of a Poisson fit is the inverse information", {
  set.seed(2)
  d <- data.frame(x = rnorm(200))
  d$y <- rpois(200, exp(0.5 + 0.4 * d$x))
  fit <- statmod(y ~ x, distributions7::poisson_distrib(), d)
  set.seed(3)
  fit <- add_bootstrap(fit, n = 200)
  se_b <- sqrt(diag(vcov(fit, vcov_type = bootstrap_vcov())))
  se_h <- sqrt(diag(vcov(fit)))
  expect_lt(max(abs(se_b / se_h - 1)), 4 / sqrt(2 * 200))
})

test_that("frequency weights give the coefficients of the duplicated rows", {
  hf <- het_fit(60, 5)
  fit <- hf$fit
  d <- hf$data
  set.seed(1)
  w <- tabulate(sample.int(60, 60, TRUE), 60)
  sp <- S7::set_props(fit@spec, weights = as.numeric(w))
  cfg <- inner_settings(fit@methods$smooth)
  r <- bootstrap_refit(fit, sp, "coefficients", fit@methods$smooth, cfg,
                       unlist(fit@coefficients, use.names = FALSE))
  dup <- statmod(y ~ x, g1, d[rep(seq_len(60), w), ])
  expect_equal(unlist(r$coefficients), unlist(dup@coefficients),
               tolerance = 1e-8)
  # and the weights of a smooth keep the knots of the fit: the replicate's
  # design is the fit's design
  set.seed(2)
  ds <- data.frame(x = runif(80))
  ds$y <- sin(2 * pi * ds$x) + rnorm(80, sd = 0.3)
  fs <- statmod(y ~ s(x, bspline_smooth(k = 8)), g1, ds)
  wi <- tabulate(sample.int(80, 80, TRUE), 80)
  sps <- S7::set_props(fs@spec, weights = as.numeric(wi))
  expect_identical(statmod_design(sps)$mu$coef_names,
                   statmod_design(fs@spec)$mu$coef_names)
  X0 <- modelterms7::term_matrix(fs@spec@terms$mu[[1L]])
  X1 <- modelterms7::term_matrix(sps@terms$mu[[1L]])
  expect_identical(X0, X1)
})

test_that("a cluster bootstrap resamples the groups", {
  set.seed(7)
  dd <- data.frame(g = rep(1:20, each = 5), x = runif(100))
  dd$y <- 1 + 2 * dd$x + rnorm(20)[dd$g] + rnorm(100, sd = 0.4)
  fit <- statmod(y ~ x, g1, dd)
  set.seed(1)
  fit <- add_bootstrap(fit, n = 200, resample = cluster_resample(~ g))
  se_b <- sqrt(diag(vcov(fit, vcov_type = bootstrap_vcov())))
  X <- cbind(1, dd$x)
  r <- as.numeric(dd$y - X %*% fit@coefficients$mu)
  sc <- rowsum(X * r, dd$g)
  A <- solve(crossprod(X))
  cl <- sqrt(diag(A %*% crossprod(sc) %*% A))
  expect_lt(max(abs(se_b[1:2] / cl - 1)), 4 / sqrt(2 * 200))
  expect_error(add_bootstrap(fit, n = 5, resample = cluster_resample(~ nope)),
               "not a column")
  # the data object has gone: it is recovered from `data =`, and wrong data
  # is rejected
  local({
    dd2 <- dd
    fx <- statmod(y ~ x, g1, dd2)
    rm(dd2)
    expect_error(add_bootstrap(fx, n = 5, resample = cluster_resample(~ g)),
                 "Pass it as cluster_resample")
    set.seed(1)
    fx <- add_bootstrap(fx, n = 20, resample = cluster_resample(~ g, data = dd))
    expect_equal(statmod_bootstrap(fx)$replicates, 20L)
    expect_error(add_bootstrap(fx, n = 30, resample = cluster_resample(
      ~ g, data = dd[1:50, ])), "Pass it as cluster_resample")
    # the stored entry does not hold the data
    e <- statmod_bootstrap(fx, resample = cluster_resample(~ g))
    expect_null(e@resample@data)
  })
})

test_that("the draws are reproducible, extendable and independent of the batch", {
  hf <- het_fit(80, 3)
  mk <- function() statmod(y ~ x, g1, hf$data)
  set.seed(4)
  a <- add_bootstrap(mk(), n = 40, resample = case_resample())
  set.seed(4)
  b <- add_bootstrap(mk(), n = 40, resample = case_resample())
  ea <- statmod_bootstrap(a, resample = case_resample())
  eb <- statmod_bootstrap(b, resample = case_resample())
  expect_identical(coef(ea), coef(eb))
  # 25 and then 40 are the 40 of one run
  set.seed(4)
  c1 <- mk()
  c1 <- add_bootstrap(c1, n = 25, resample = case_resample())
  c1 <- add_bootstrap(c1, n = 40, resample = case_resample())
  expect_identical(coef(statmod_bootstrap(c1, resample = case_resample())),
                   coef(ea))
  # n below the stored number computes nothing and keeps the entry
  c1 <- add_bootstrap(c1, n = 10, resample = case_resample())
  expect_equal(statmod_bootstrap(c1)$replicates, 40L)
  # the first rows are those of a fresh run of that size
  set.seed(4)
  d10 <- add_bootstrap(mk(), n = 10, resample = case_resample())
  expect_identical(coef(statmod_bootstrap(d10, resample = case_resample())),
                   coef(ea)[1:10, , drop = FALSE])
  # the parametric scheme too
  set.seed(5)
  p1 <- add_bootstrap(mk(), n = 12)
  set.seed(5)
  p2 <- add_bootstrap(add_bootstrap(mk(), n = 5), n = 12)
  expect_identical(coef(statmod_bootstrap(p1)), coef(statmod_bootstrap(p2)))
})

test_that("the replicates do not depend on the number of workers", {
  skip_if(exists(".__DEVTOOLS__", asNamespace("statmodels7")),
          "worker_map() is sequential under pkgload")
  hf <- het_fit(80, 3)
  run <- function(w) {
    fit <- statmod(y ~ x, g1, hf$data, threads = numericals7::n_threads(workers = w))
    set.seed(8)
    coef(statmod_bootstrap(add_bootstrap(fit, n = 24,
                                         resample = case_resample())))
  }
  expect_identical(run(1L), run(3L))
})

test_that("a replica that fails is replaced and counted by cause", {
  hf <- het_fit(60, 2)
  fit <- hf$fit
  orig <- bootstrap_refit
  calls <- 0L
  local_mocked_bindings(bootstrap_refit = function(...) {
    calls <<- calls + 1L
    if (calls %% 4L == 0L) stop("forced")
    orig(...)
  })
  set.seed(1)
  fit <- add_bootstrap(fit, n = 12, resample = case_resample())
  e <- statmod_bootstrap(fit, resample = case_resample())
  expect_equal(nrow(coef(e)), 12L)
  expect_equal(e@n, 12L)
  tab <- statmod_bootstrap(fit)
  expect_equal(tab$replicates, 12L)
  expect_equal(tab$error, 3)
  expect_equal(tab$mode_error + tab$unidentified, 0)
  expect_equal(tab$replaced, 3)
  expect_equal(sum(e@status == 0L), 12L)
})

test_that("a replica whose mode error exceeds the limit counts as failed", {
  hf <- het_fit(60, 2)
  fit <- hf$fit
  orig <- boot_read_mode
  calls <- 0L
  local_mocked_bindings(boot_read_mode = function(...) {
    r <- orig(...)
    calls <<- calls + 1L
    # the first call reads the fit itself
    if (calls > 1L && calls %% 3L == 1L) r$me <- 1
    r
  })
  set.seed(1)
  fit <- add_bootstrap(fit, n = 10, resample = case_resample())
  tab <- statmod_bootstrap(fit)
  expect_equal(tab$mode_error, 4)
  expect_equal(tab$replicates, 10L)
})

test_that("a coefficient left unidentified by the weights fails the replica", {
  set.seed(12)
  n <- 120
  d <- data.frame(g = factor(c(rep("a", n - 3), rep("rare", 3))),
                  x = rnorm(n))
  d$y <- rbinom(n, 1, plogis(-0.3 + 0.8 * d$x))
  fit <- quiet(statmod(y ~ g + x, distributions7::bernoulli_distrib(), d))
  set.seed(2)
  fit <- quiet(add_bootstrap(fit, n = 60, resample = case_resample()))
  tab <- statmod_bootstrap(fit)
  expect_gt(tab$replaced, 0)
  expect_equal(tab$error, 0)
  # the replicates kept all carry a coefficient for the rare level that the
  # observations identify
  b <- coef(statmod_bootstrap(fit, resample = case_resample()))
  expect_equal(nrow(b), 60L)
  expect_true(all(is.finite(b)))
})

test_that("a bootstrap that fails more than it succeeds stops", {
  hf <- het_fit(60, 2)
  local_mocked_bindings(bootstrap_refit = function(...) stop("forced"))
  expect_error(add_bootstrap(hf$fit, n = 10, resample = case_resample()),
               "More bootstrap replicas failed")
})

test_that("the settings the bootstrap does not cover are rejected", {
  set.seed(6)
  gg <- data.frame(g = factor(rep(1:8, each = 5)), x = rnorm(40))
  gg$y <- 1 + gg$x + rnorm(8)[gg$g] + rnorm(40, sd = 0.5)
  fr <- statmod(y ~ x + random(~ 1 | g), g1, gg)
  expect_error(add_bootstrap(fr, n = 3, resample = case_resample()),
               "cluster_resample")
  expect_equal(statmod_bootstrap(add_bootstrap(fr, n = 3,
    resample = cluster_resample(~ g)))$replicates, 3L)
  hf <- het_fit(60, 2)
  set.seed(3)
  dl <- data.frame(matrix(rnorm(300), 60))
  dl$y <- 1 + dl$X1 + rnorm(60)
  fcv <- quiet(statmod(y ~ lasso(~ X1 + X2 + X3 + X4 + X5), g1, dl,
                       sparse_criterion = cv(nfolds = 3)))
  expect_error(add_bootstrap(fcv, n = 3, refit = "full"), "cross-validated")
  expect_error(add_bootstrap(hf$fit, n = 0), "positive whole number")
  expect_error(add_bootstrap(hf$fit, n = 3, resample = "case"),
               "resample_method")
  expect_error(cluster_resample(y ~ x), "one-sided")
  expect_error(bootstrap_vcov(refit = "partial"), "coefficients")
  # a model that moves a block with its coefficients is bootstrapped, see
  # test-bootstrap-structural.R
  set.seed(8)
  ds <- data.frame(x = runif(120))
  ds$y <- ifelse(ds$x < 0.5, ds$x, 0.5 + 2 * (ds$x - 0.5)) + rnorm(120, sd = 0.1)
  fsg <- quiet(statmod(y ~ seg(x, psi = 0.4), g1, ds))
  expect_no_error(add_bootstrap(fsg, n = 2, resample = case_resample()))
  # a replicate set belongs to the fit that computed it
  fit <- hf$fit
  set.seed(1)
  fit <- add_bootstrap(fit, n = 5, resample = case_resample())
  moved <- S7::set_props(fit, coefficients = lapply(fit@coefficients,
                                                    function(v) v + 0.01))
  expect_equal(nrow(statmod_bootstrap(moved)), 0L)
  expect_equal(nrow(statmod_bootstrap(fit)), 1L)
})

test_that("a variance resolves against the stored entries", {
  hf <- het_fit(60, 2)
  fit <- hf$fit
  # nothing stored: the default is computed at the call, with a message
  expect_message(V0 <- vcov(fit, vcov_type = bootstrap_vcov(n = 6)),
                 "computed 6 replicates")
  expect_equal(statmod_bootstrap(fit)$resample, "parametric")
  # no arguments and one entry: the entry is used and nothing is computed
  expect_no_message(V1 <- vcov(fit, vcov_type = bootstrap_vcov()))
  expect_identical(V0, V1)
  # an argument that differs computes another entry
  expect_message(vcov(fit, vcov_type = bootstrap_vcov(
    n = 5, resample = case_resample())), "case")
  expect_equal(nrow(statmod_bootstrap(fit)), 2L)
  # several entries and no argument: an error that prints the table
  expect_error(vcov(fit, vcov_type = bootstrap_vcov()), "case")
  expect_error(statmod_bootstrap(fit, refit = "coefficients"), "parametric")
  expect_error(statmod_bootstrap(fit, resample = cluster_resample(~ x)),
               "No stored bootstrap")
  # n below the stored number uses the first rows
  e <- statmod_bootstrap(fit, resample = parametric_resample())
  V3 <- vcov(fit, vcov_type = bootstrap_vcov(
    n = 4, resample = parametric_resample()))
  expect_equal(unname(V3), unname(stats::cov(coef(e)[1:4, ])))
})

test_that("summary and confint compute the replicates once", {
  hf <- het_fit(60, 2)
  fit <- hf$fit
  orig <- bootstrap_engine
  calls <- 0L
  local_mocked_bindings(bootstrap_engine = function(...) {
    calls <<- calls + 1L
    orig(...)
  })
  s <- quiet(summary(fit, vcov_type = bootstrap_vcov(
    n = 8, resample = case_resample())))
  expect_equal(calls, 1L)
  expect_equal(s@bootstrap$n, 8L)
  expect_output(print(s), "bootstrap, 8 replicates")
  expect_output(print(s), "MC error about 25.0%")
  ci <- confint(fit, vcov_type = bootstrap_vcov(
    n = 8, resample = case_resample()))
  expect_equal(calls, 1L)
  expect_equal(ci$se, unname(sqrt(diag(vcov(fit, vcov_type = bootstrap_vcov(
    n = 8, resample = case_resample()))))))
  # a second variance setting computes its own
  confint(fit, vcov_type = bootstrap_vcov(n = 5, resample = case_resample(),
                                          refit = "coefficients"))
  expect_equal(calls, 1L)
})

test_that("a coordinate a lasso sets to zero has a bootstrap variance", {
  set.seed(4)
  n <- 100
  dl <- data.frame(matrix(rnorm(n * 6), n))
  dl$y <- 1 + dl$X1 + rnorm(n)
  fl <- quiet(statmod(y ~ lasso(~ X1 + X2 + X3 + X4 + X5 + X6), g1, dl,
                      outer_criterion = bic()))
  zero <- unlist(fl@coefficients$mu) == 0
  expect_true(any(zero))
  set.seed(1)
  fl <- quiet(add_bootstrap(fl, n = 30, resample = case_resample()))
  V <- vcov(fl, vcov_type = bootstrap_vcov())
  expect_true(all(is.finite(diag(V))))
  # the other conventions leave the zero coordinates missing
  expect_true(anyNA(diag(vcov(fl))))
  s <- quiet(summary(fl, vcov_type = bootstrap_vcov()))
  fr <- s@bootstrap$selection
  expect_equal(length(fr), 6L)
  expect_true(all(fr >= 0 & fr <= 1))
  expect_output(print(s), "selection frequency")
  # readable = TRUE is the covariance of the replicates, reordered
  Vr <- vcov(fl, vcov_type = bootstrap_vcov(), readable = FALSE)
  expect_equal(unname(Vr), unname(stats::cov(coef(statmod_bootstrap(
    fl, resample = case_resample())))))
  expect_equal(unname(V), unname(Vr))
})

test_that("a random-effect coordinate stores the prediction error", {
  set.seed(6)
  gg <- data.frame(g = factor(rep(1:10, each = 5)), x = rnorm(50))
  gg$y <- 1 + gg$x + rnorm(10, sd = 1)[gg$g] + rnorm(50, sd = 0.5)
  fr <- statmod(y ~ x + random(~ 1 | g), g1, gg)
  set.seed(1)
  fr <- add_bootstrap(fr, n = 60)
  e <- statmod_bootstrap(fr)
  expect_equal(e$replicates, 60L)
  ent <- statmod_bootstrap(fr, resample = parametric_resample())
  expect_identical(ent@random, "redrawn")
  Vb <- vcov(fr, vcov_type = bootstrap_vcov())
  Vm <- vcov(fr)
  i <- grep("random", rownames(Vb))
  # the effects' standard deviation is the model's prediction error within
  # the Monte Carlo error of 60 replicates
  expect_lt(abs(mean(sqrt(diag(Vb)[i])) / mean(sqrt(diag(Vm)[i])) - 1),
            4 / sqrt(2 * 60))
})

test_that("predictive_bootstrap is the loop it replaced, to the bit", {
  set.seed(5)
  d <- data.frame(x = runif(50))
  d$y <- rpois(50, exp(1 + d$x))
  fit <- statmod(y ~ x, distributions7::poisson_distrib(), d)
  spec0 <- fit@spec
  design0 <- statmod_design(spec0)
  nd <- statmod_respec(spec0, data.frame(x = c(0.3, 0.6)), FALSE)
  design <- statmod_design(nd)
  aside <- random_modes(spec0, "zero")
  aside <- aside[aside$mode != "conditional", , drop = FALSE]
  for (rf in c("coefficients", "full")) {
    # the reference is the loop of the version before the engine
    set.seed(9)
    draw <- bootstrap_simulator(fit, spec0, design0, aside)
    cfg <- inner_settings(fit@methods$smooth)
    beta0 <- unlist(fit@coefficients, use.names = FALSE)
    ref <- list()
    for (b in 1:6) {
      y <- draw()
      r <- bootstrap_refit(fit, S7::set_props(spec0, response = y), rf,
                           fit@methods$smooth, cfg, beta0)
      fit_b <- S7::set_props(fit, coefficients = r$coefficients,
                             hyper = r$hyper)
      cb <- predictive_mixture(fit_b, nd, design, aside,
                               statmod_eta(nd, design, r$coefficients)$eta,
                               NULL, weight = 1,
                               n_draw = max(1L, ceiling(5000 / 6)))
      ref <- c(ref, cb)
    }
    ref <- lapply(ref, function(cc) {
      cc$weight <- cc$weight / 6
      cc
    })
    set.seed(9)
    got <- predictive_bootstrap(fit, nd, design, aside, 6L, rf)
    attr(got, "infinite_variance") <- NULL
    expect_identical(got, ref, label = rf)
  }
})
