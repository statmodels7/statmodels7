# nl(curvature =) on the inner joint Newton beside a filter: the three
# settings reach the same mode, "auto" falls back to the Gauss-Newton matrix
# where the full one has no Cholesky factor, and the outer gradient does not
# depend on the setting.

nl_curvature_data <- function() {
  set.seed(31)
  S <- 4L
  Tn <- 30L
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

test_that("auto falls back where the full matrix has no factor", {
  d <- nl_curvature_data()
  f <- nl_curvature_fit("auto", d)
  sp <- f@spec
  des <- statmod_design(sp)
  obj <- statmod_objective(sp, f@hyper, des, FALSE)
  jp <- statmod_joint_pieces(sp, des, obj, f@hyper)
  nb <- jp$nb
  u0 <- c(unlist(f@coefficients[sp@distrib@params], use.names = FALSE),
          jp$zeta())
  # a point far from the mode in the nonlinear parameters, where the score
  # times the block's second derivative makes the full matrix indefinite
  bn <- grep("^mu\\.nl", unlist(lapply(sp@distrib@params, function(p)
    paste0(p, ".", des[[p]]$coef_names))))
  found <- FALSE
  for (s in c(2, 4, 8, -2, -4)) {
    u <- u0
    u[bn] <- u[bn] + s
    jp$curvature$fallbacks[] <- 0L
    H <- jp$he(u)
    if (jp$curvature$fallbacks[["cholesky"]] > 0L) {
      found <- TRUE
      break
    }
  }
  expect_true(found)
  # the fallback is the Gauss-Newton matrix
  jp$curvature$gauss_newton <- TRUE
  expect_equal(H, jp$he(u), tolerance = 1e-12)
  jp$curvature$gauss_newton <- FALSE
  # and at the mode the full matrix is taken
  jp$curvature$fallbacks[] <- 0L
  Hm <- jp$he(u0)
  expect_identical(jp$curvature$fallbacks[["cholesky"]], 0L)
  jp$curvature$gauss_newton <- TRUE
  expect_gt(max(abs(Hm - jp$he(u0))), 1e-6)
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
