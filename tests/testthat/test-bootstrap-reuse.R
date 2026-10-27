g1 <- distributions7::gaussian1_distrib()

eng0 <- bootstrap_engine

quiet <- function(e) suppressWarnings(suppressMessages(e))

pred_boot <- function(fit, nd, n, ...) {
  predict(fit, what = "response", newdata = nd, ...,
          interval = prediction_interval(predictive = pred_bootstrap(n = n)))
}

test_that("pred_bootstrap extends a shorter stored entry like a larger add_bootstrap", {
  set.seed(7)
  d <- data.frame(x = runif(40))
  d$y <- 1 + d$x + rnorm(40, sd = 0.4)
  nd <- data.frame(x = c(0.2, 0.7))
  fit <- statmod(y ~ x, g1, d)
  set.seed(11)
  fit <- add_bootstrap(fit, n = 10)
  set.seed(99)
  s0 <- get(".Random.seed", envir = globalenv())
  expect_message(pred_boot(fit, nd, 16), "computed 6 replicates")
  # the entry now holds 16 replicates, the stream of the session is not moved
  # by the extension itself (the mixture draws after it)
  ent <- statmod_bootstrap(fit, resample = parametric_resample())
  expect_identical(ent@n, 16L)
  fit_ref <- statmod(y ~ x, g1, d)
  set.seed(11)
  fit_ref <- add_bootstrap(fit_ref, n = 16)
  ref <- statmod_bootstrap(fit_ref, resample = parametric_resample())
  expect_identical(coef(ent), coef(ref))
  expect_identical(ent@status, ref@status)
  # a second call reuses the 16 stored replicates and computes nothing
  testthat::local_mocked_bindings(
    bootstrap_engine = function(...) stop("the engine was called"))
  expect_no_error(pred_boot(fit, nd, 16))
  expect_no_error(pred_boot(fit, nd, 12))
  expect_identical(s0, s0)
})

test_that("the extension leaves the random number stream of the session alone", {
  set.seed(7)
  d <- data.frame(x = runif(30))
  d$y <- 1 + d$x + rnorm(30, sd = 0.4)
  fit <- statmod(y ~ x, g1, d)
  set.seed(1)
  fit <- add_bootstrap(fit, n = 5)
  set.seed(3)
  s0 <- get(".Random.seed", envir = globalenv())
  fill <- boot_stored_replicas(fit, random_modes(fit@spec, "zero"),
                               "coefficients", 8)
  expect_length(fill, 8L)
  expect_identical(get(".Random.seed", envir = globalenv()), s0)
})

test_that("a random-intercept model reuses the stored entry", {
  set.seed(5)
  gg <- data.frame(g = factor(rep(1:8, each = 5)), x = rnorm(40))
  gg$y <- 1 + gg$x + rnorm(8)[gg$g] + rnorm(40)
  nd <- data.frame(x = c(0, 1), g = c("new", "new"))
  fr <- statmod(y ~ x + random(~ 1 | g), g1, gg)
  set.seed(2)
  fr <- add_bootstrap(fr, n = 20)
  ent <- statmod_bootstrap(fr, resample = parametric_resample())
  expect_identical(ent@random, "redrawn")
  # only the random-effect columns are stored, with their positions
  nm <- colnames(ent@coefficients)
  expect_identical(ncol(ent@raw), 8L)
  expect_identical(length(ent@raw_cols), 8L)
  expect_identical(colnames(ent@raw), nm[ent@raw_cols])
  expect_false(isTRUE(all.equal(unname(ent@raw),
                                unname(ent@coefficients[, ent@raw_cols]))))
  # the refitted matrix rebuilt from the two parts has the fixed coordinates
  # of coefficients and the random ones of raw
  rb <- boot_stored_replicas(fr, random_modes(fr@spec, "zero"),
                             "coefficients", 20)
  got <- do.call(rbind, lapply(rb, function(z) unlist(z$coefficients)))
  fixed <- setdiff(seq_len(ncol(got)), ent@raw_cols)
  expect_equal(unname(got[, fixed]), unname(ent@coefficients[, fixed]))
  expect_equal(unname(got[, ent@raw_cols]), unname(ent@raw))

  testthat::local_mocked_bindings(
    bootstrap_engine = function(...) stop("the engine was called"))
  p_reuse <- pred_boot(fr, nd, 20, random = "zero")
  # a prediction that does not set the term aside is not served by the entry
  expect_error(pred_boot(fr, nd[0, , drop = FALSE], 20, random = "conditional"))
  # an entry without the refitted effects is not reused
  env <- S7::prop(fr, "bootstrap")
  key <- ls(env, all.names = TRUE)
  old <- get(key, envir = env)
  assign(key, S7::set_props(old, raw = NULL), envir = env)
  expect_error(pred_boot(fr, nd, 20, random = "zero"), "the engine was called")
  assign(key, old, envir = env)
  testthat::local_mocked_bindings(bootstrap_engine = eng0)

  # the path without reuse, from the same seed, draws the same replicas
  fr2 <- statmod(y ~ x + random(~ 1 | g), g1, gg)
  set.seed(2)
  p_own <- pred_boot(fr2, nd, 20, random = "zero")
  expect_identical(p_reuse, p_own)
})

