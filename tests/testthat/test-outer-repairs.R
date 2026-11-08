# Repairs of the R route of outer_fit() made before its port to C++
# (piano_motore.md, D44): each block fails with the code before the repair.
# The repairs are of the R route, and the injections reach its functions only
# there, so each block runs the search in R (outer_search_settings() NULL).

local_outer_r <- function(env = parent.frame()) {
  local_mocked_bindings(outer_search_settings = function(...) NULL,
                        .package = "statmodels7", .env = env)
}

outer_repair_data <- function(n = 300, seed = 11) {
  set.seed(seed)
  x <- runif(n)
  data.frame(x = x, z = runif(n), y = sin(2 * pi * x) + rnorm(n, sd = 0.3))
}

outer_crit <- function(fit) {
  h <- fit@history$outer
  h$criterion[nrow(h)]
}

test_that("a NULL exact gradient makes the point unavailable (D5)", {
  local_outer_r()
  d <- outer_repair_data()
  ref <- statmod(y ~ s(x), distrib = gaussian1_distrib(), data = d)
  orig <- statmod_marginal_grad
  k <- 0L
  local_mocked_bindings(statmod_marginal_grad = function(...) {
    k <<- k + 1L
    if (k == 3L) return(NULL)
    orig(...)
  }, .package = "statmodels7")
  fit <- statmod(y ~ s(x), distrib = gaussian1_distrib(), data = d)
  expect_gte(k, 3L)
  # read as zeros, the third gradient stopped the search short of the
  # optimum with the fit reported as converged
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-6)
  expect_equal(outer_crit(fit), outer_crit(ref), tolerance = 1e-8)
})

test_that("an error in the inner score at a usable point makes it unavailable (D6)", {
  local_outer_r()
  d <- outer_repair_data()
  ref <- statmod(y ~ s(x), distrib = gaussian1_distrib(), data = d)
  orig <- statmod_alternate
  k <- 0L
  local_mocked_bindings(statmod_alternate = function(...) {
    r <- orig(...)
    k <<- k + 1L
    if (k == 4L) r$obj$gr <- function(p) stop("injected gradient error")
    r
  }, .package = "statmodels7")
  fit <- expect_no_error(statmod(y ~ s(x), distrib = gaussian1_distrib(), data = d))
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-6)
})

test_that("the pinned branch of the first evaluation leaves the caller's m alone (D8)", {
  local_outer_r()
  d <- outer_repair_data()
  assign("m", "the user's m", envir = globalenv())
  on.exit(rm("m", envir = globalenv()), add = TRUE)
  orig <- statmod_marginal
  k <- 0L
  local_mocked_bindings(
    statmod_marginal = function(...) {
      k <<- k + 1L
      if (k == 1L) return(NULL)
      orig(...)
    },
    deficient_coords = function(K, held = integer(0)) 1L,
    .package = "statmodels7")
  fit <- tryCatch(statmod(y ~ s(x), distrib = gaussian1_distrib(), data = d),
                  error = function(e) NULL)
  expect_identical(get("m", envir = globalenv()), "the user's m")
})

test_that("the trace prints the gradient of the point on its own line (D1)", {
  local_outer_r()
  d <- outer_repair_data()
  out <- capture.output(statmod(y ~ s(x), distrib = gaussian1_distrib(), data = d,
                                verbose = c(outer = TRUE)))
  # the evaluation lines: number, criterion, |grad|max, the hyperparameters
  ev <- grep("^ +[0-9]+ +-?[0-9]+[.][0-9]+ +[-0-9.e+]+ +", out, value = TRUE)
  expect_gt(length(ev), 1L)
  first <- strsplit(trimws(ev[1L]), " +")[[1L]]
  # the first line printed the previous evaluation's gradient, of which
  # there is none, as 0
  expect_false(identical(as.numeric(first[3L]), 0))
})

test_that("an unavailable line shows the estimated coefficients (D3)", {
  local_outer_r()
  d <- outer_repair_data()
  # the second evaluation's inner fit fails, so that point is unavailable;
  # the criterion also estimates sigma2's coefficients
  orig <- statmod_alternate
  k <- 0L
  local_mocked_bindings(statmod_alternate = function(...) {
    k <<- k + 1L
    if (k == 2L) stop("injected inner failure")
    orig(...)
  }, .package = "statmodels7")
  out <- capture.output(statmod(y ~ s(x) | sigma2 ~ z, distrib = gaussian2_distrib(),
                                data = d, verbose = c(outer = TRUE)))
  un <- grep("unavailable", out, value = TRUE)
  expect_gte(length(un), 1L)
  expect_false(any(grepl("\\bNA\\b", un)))
})

test_that("the resolution does not move a coordinate pinned as deficient (D4)", {
  local_outer_r()
  d <- outer_repair_data()
  orig_m <- statmod_marginal
  orig_r <- criterion_resolution
  k <- 0L
  scores <- list()
  # the first criterion is unavailable and coordinate 1 is named deficient,
  # so the first evaluation pins it
  local_mocked_bindings(
    statmod_marginal = function(...) {
      k <<- k + 1L
      if (k == 1L) return(NULL)
      orig_m(...)
    },
    deficient_coords = function(K, held = integer(0)) 1L,
    criterion_resolution = function(st, ...) {
      scores[[length(scores) + 1L]] <<- st$score
      orig_r(st, ...)
    },
    .package = "statmodels7")
  fit <- statmod(y ~ s(x), distrib = gaussian1_distrib(), data = d)
  expect_identical(as.integer(fit@methods$pinned), 1L)
  expect_gt(length(scores), 0L)
  expect_true(all(vapply(scores, function(s) identical(s[[1L]], 0), logical(1))))
})

test_that("the edge probe's readings leave the resolution as it was (D7)", {
  local_outer_r()
  d <- outer_repair_data()
  orig_r <- criterion_resolution
  orig_p <- edge_probe_hyper
  probing <- FALSE
  # every coordinate is probed, and every reading taken during the probe is
  # made large
  local_mocked_bindings(
    edge_probe_hyper = function(par, start, nh, value, fn, edge = 8,
                                tol = mode_error_limit()) {
      probing <<- TRUE
      on.exit(probing <<- FALSE)
      orig_p(par, start, nh, value, fn, edge = -1, tol = Inf)
    },
    criterion_resolution = function(st, ...) {
      r <- orig_r(st, ...)
      if (probing) 1e3 else r
    },
    .package = "statmodels7")
  fit <- statmod(y ~ s(x), distrib = gaussian1_distrib(), data = d)
  expect_lt(fit@methods$search@line_search@resolution, 1)
})
