# nl(curvature =) on the inner joint Newton beside a filter: the three
# settings reach the same mode, "auto" falls back to the Gauss-Newton matrix
# after a failed line search, and the outer gradient does not depend on the
# setting. Then the same on the route without a filter (iwls_fit()).

nl_curvature_data <- function() {
  set.seed(31)
  S <- 4L
  Tn <- 25L
  n <- S * Tn
  d <- data.frame(t = rep(seq_len(Tn), S), id = factor(rep(seq_len(S), each = Tn)),
                  xs = runif(n, 0, 3))
  e <- numeric(n)
  for (i in 2:n) e[i] <- if (d$t[i] == 1) 0 else 0.7 * e[i - 1] + rnorm(1, sd = 0.4)
  d$y <- 1.5 * exp(-0.8 * d$xs) + e + rnorm(n, sd = 0.5)
  d
}

nl_curvature_cache <- new.env()

nl_curvature_fit <- function(cv, d) {
  if (!is.null(nl_curvature_cache[[cv]])) return(nl_curvature_cache[[cv]])
  fm <- eval(bquote(y ~ nl(~ b * exp(-r * xs), start = list(b = 1, r = 0.5),
                           curvature = .(cv)) +
                      gas(p = 1, q = 1, omega ~ random(~ 1 | id), time = t,
                          by = id)))
  f <- statmod(fm, distributions7::gaussian1_distrib(), data = d,
               outer_criterion = NULL)
  assign(cv, f, envir = nl_curvature_cache)
  f
}

test_that("the three curvature settings reach the same mode", {
  d <- nl_curvature_data()
  fits <- lapply(c(auto = "auto", full = "full", gn = "gauss_newton"),
                 nl_curvature_fit, d = d)
  ll <- vapply(fits, function(f) as.numeric(logLik(f)), numeric(1))
  expect_lt(max(abs(ll - ll[["auto"]])), 1e-6)
  cf <- lapply(fits, function(f) unlist(f@coefficients))
  expect_lt(max(abs(cf$full - cf$auto)), 1e-4)
  expect_lt(max(abs(cf$gn - cf$auto)), 1e-4)
})

test_that("the curvature units follow the setting", {
  d <- nl_curvature_data()
  for (cv in c("auto", "full", "gauss_newton")) {
    fm <- eval(bquote(y ~ nl(~ b * exp(-r * xs), start = list(b = 1, r = 0.5),
                             curvature = .(cv)) +
                        gas(p = 1, q = 1, time = t, by = id)))
    sp <- statmod_spec(fm, distributions7::gaussian1_distrib(), data = d)
    u <- nl_curvature_units(statmod_design(sp))
    if (cv == "gauss_newton") {
      expect_length(u, 0L)
    } else {
      expect_length(u, 1L)
      expect_identical(u[[1L]]$curvature, cv)
    }
  }
})

test_that("auto falls back to Gauss-Newton after a failed line search", {
  # No measured fit has had its line search fail on the full step, so the
  # failure is built deliberately: an Armijo search that accepts only a
  # decrease of 0.999999 of the linear prediction and stops after one
  # halving rejects the full step from the start of the run. The run is
  # then repeated on the Gauss-Newton matrix, which the counter records.
  d <- nl_curvature_data()
  f <- nl_curvature_fit("auto", d)
  sp <- f@spec
  des <- statmod_design(sp)
  obj <- statmod_objective(sp, f@hyper, des, FALSE)
  beta <- unlist(lapply(sp@distrib@params, function(p) {
    b <- f@coefficients[[p]]
    b + 0.3
  }), use.names = FALSE)
  strict <- optimizers7::newton(
    line_search = optimizers7::armijo(c1 = 0.999999, max_step = 1))
  r <- statmod_fit_joint(sp, des, obj, beta, f@hyper, optimizer = strict)
  expect_identical(r$curvature_fallbacks[["search"]], 1L)
  # the Gauss-Newton matrix is the full one without the block's term
  jp <- statmod_joint_pieces(sp, des, obj, f@hyper)
  u0 <- c(unlist(f@coefficients[sp@distrib@params], use.names = FALSE),
          jp$zeta())
  Hf <- jp$he(u0)
  jp$curvature$gauss_newton <- TRUE
  Hg <- jp$he(u0)
  expect_gt(max(abs(Hf - Hg)), 1e-6)
  expect_equal(Hf[-seq_len(jp$nb), ], Hg[-seq_len(jp$nb), ], tolerance = 1e-12)
})

test_that("the outer gradient does not depend on the setting", {
  # each at its own fitted mode: the setting changes the path of the inner
  # fit and not its optimum
  d <- nl_curvature_data()
  m <- reml(information = distributions7::observed_info())
  grad_at <- function(f) {
    sp <- f@spec
    des <- statmod_design(sp)
    idx <- outer_hyper_index(sp, statmod_blocks(sp, des))
    statmod_marginal_grad(sp, des, f@coefficients, f@hyper, m, idx)
  }
  ga <- grad_at(nl_curvature_fit("auto", d))
  gg <- grad_at(nl_curvature_fit("gauss_newton", d))
  expect_true(all(is.finite(ga)))
  expect_lt(max(abs(ga - gg)), 1e-6)
})


# the route without a filter: every iwls_fit() call of one statmod() fit
nl_iwls_fit <- function(cv, fail_full = FALSE) {
  puro <- Puromycin[Puromycin$state == "treated", ]
  runs <- list()
  orig_fit <- iwls_fit
  orig_step <- iwls_newton_step
  testthat::local_mocked_bindings(iwls_fit = function(...) {
    r <- orig_fit(...)
    runs[[length(runs) + 1L]] <<- r
    r
  })
  if (fail_full) {
    # the full step's line search fails at every iteration: the Newton step
    # taken first is reported as finding no acceptable length
    testthat::local_mocked_bindings(iwls_newton_step = function(obj, beta,
                                                                value, g, H,
                                                                frozen,
                                                                halvings,
                                                                exact = FALSE) {
      if (exact) {
        return(list(ok = FALSE, cand = beta, vnew = value, step_used = 0))
      }
      orig_step(obj, beta, value, g, H, frozen, halvings)
    })
  }
  fm <- eval(bquote(rate ~ 0 + nl(~ Vm * conc / (K + conc),
                                  links = list(Vm = linkfunctions7::log_link(),
                                               K = linkfunctions7::log_link()),
                                  start = list(Vm = 150, K = 0.1),
                                  curvature = .(cv))))
  f <- statmod(fm, distributions7::gaussian1_distrib(), data = puro,
               outer_criterion = NULL)
  cnt <- Reduce(`+`, lapply(runs, function(r) r$curvature))
  list(fit = f, curvature = cnt, runs = runs)
}

test_that("without a filter the three settings reach the same mode", {
  r <- lapply(c(auto = "auto", full = "full", gn = "gauss_newton"),
              nl_iwls_fit)
  cf <- lapply(r, function(x) unlist(x$fit@coefficients))
  expect_lt(max(abs(cf$full - cf$gn)), 1e-6)
  expect_lt(max(abs(cf$auto - cf$gn)), 1e-6)
  # the full step is taken first under "full" and "auto", never under
  # "gauss_newton", and no line search failed
  expect_gt(r$full$curvature[["full"]], 0L)
  expect_identical(r$auto$curvature, r$full$curvature)
  expect_identical(r$gn$curvature, c(full = 0L, search = 0L))
})

test_that("without a filter a failed full step falls back to scoring", {
  ra <- nl_iwls_fit("auto", fail_full = TRUE)
  rf <- nl_iwls_fit("full", fail_full = TRUE)
  rg <- nl_iwls_fit("gauss_newton")
  # "auto" leaves the full step after its first failure in each run;
  # "full" tries it at every iteration and takes the scoring step instead
  fails_a <- vapply(ra$runs, function(r) r$curvature[["search"]], integer(1))
  expect_true(all(fails_a <= 1L))
  expect_gt(sum(fails_a), 0L)
  expect_identical(ra$curvature[["full"]], 0L)
  fails_f <- vapply(rf$runs, function(r) r$curvature[["search"]], integer(1))
  expect_gt(max(fails_f), 1L)
  expect_identical(rf$curvature[["full"]], 0L)
  # the scoring steps reach the mode of "gauss_newton"
  cg <- unlist(rg$fit@coefficients)
  expect_lt(max(abs(unlist(ra$fit@coefficients) - cg)), 1e-6)
  expect_lt(max(abs(unlist(rf$fit@coefficients) - cg)), 1e-6)
})
