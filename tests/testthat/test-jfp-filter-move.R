# A Jeffreys block outside a filter, beside a term whose block moves in the
# equation that the filter drives: the block's derivatives reach the filtered
# predictors through the filter (jfp_static_rows(), jfp_static_curv()), and
# the penalty's Hessian and its derivatives along directions are exact

jfp_filter_move_case <- function(fm) {
  set.seed(9)
  S <- 4L
  Tn <- 30L
  n <- S * Tn
  d <- data.frame(t = rep(1:Tn, S), id = factor(rep(1:S, each = Tn)),
                  xs = stats::runif(n, 0, 3), z = stats::rnorm(n))
  e <- numeric(n)
  for (i in 2:n) {
    e[i] <- if (d$t[i] == 1) 0 else 0.7 * e[i - 1] + stats::rnorm(1, sd = 0.5)
  }
  lev <- stats::rnorm(S, sd = 0.5)
  d$yn <- stats::rnbinom(n, mu = exp(0.5 + 1.5 * exp(-0.8 * d$xs) +
                                       0.4 * lev[d$id] + 0.5 * e),
                         size = exp(1 + 0.4 * d$z))
  d$ys <- stats::rnbinom(n, mu = exp(0.5 + 0.3 * d$xs +
                                       0.6 * pmax(d$xs - 1.5, 0) +
                                       0.4 * lev[d$id] + 0.5 * e),
                         size = exp(1 + 0.4 * d$z))
  # the hyperparameters at their start and three Newton steps: the
  # derivatives are read at a point away from the mode
  f <- suppressWarnings(statmod(fm, distrib = distributions7::negbin2_distrib(),
                                data = d, outer_criterion = NULL,
                                inner_optimizer = optimizers7::newton(
                                  maxit = 3)))
  sp <- f@spec
  des <- statmod_design(sp)
  sst <- statmod_structural_state(des)
  key <- names(sst$zeta)[1]
  free <- setdiff(names(sst$zeta[[key]]), sst$held[[key]])
  nb <- sum(vapply(des, function(z) z$npar, integer(1)))
  obj <- statmod_objective(sp, f@hyper, des)
  jp <- statmod_joint_pieces(sp, des, obj, f@hyper)
  u <- c(unlist(f@coefficients[sp@distrib@params], use.names = FALSE),
         as.numeric(sst$zeta[[key]][free]))
  cf_at <- function(uu) {
    jp$setz(uu[nb + seq_along(free)])
    obj$split(uu[seq_len(nb)])
  }
  list(sp = sp, des = des, u = u, cf_at = cf_at)
}

test_that("a Jeffreys penalty beside a filter and a moving block is exact", {
  gasr <- quote(gas(p = 1, q = 1, omega ~ 1 + random(~ 1 | id), time = t,
                    by = id))
  fms <- list(
    bquote(yn ~ nl(~ b * exp(-r * xs), start = list(b = 1, r = 0.5)) +
             .(gasr) | theta ~ jfp(~ z)),
    bquote(ys ~ seg(xs, smoothed = numericals7::smooth_probit()) + .(gasr) |
             theta ~ jfp(~ z)))
  rel <- function(a, b) max(abs(a - b)) / max(abs(b))
  cd <- function(f, x, v, h) (f(x + h * v) - f(x - h * v)) / (2 * h)
  rich <- function(f, x, v, h) (4 * cd(f, x, v, h / 2) - cd(f, x, v, h)) / 3
  for (fm in fms) {
    k <- jfp_filter_move_case(eval(fm))
    sp <- k$sp
    des <- k$des
    expect_true(jfp_outer_ok(sp, des, 1L))
    expect_true(jfp_outer_ok(sp, des, 2L))
    set.seed(2)
    u <- k$u + 0.05 * stats::rnorm(length(k$u))
    v <- stats::rnorm(length(u))
    w <- stats::rnorm(length(u))
    pen <- function(x, what) jfp_penalty(sp, des, k$cf_at(x), what,
                                         joint = TRUE)
    D3 <- function(x) jfp_penalty_move(sp, des, k$cf_at(x), v)
    errs <- function() {
      H <- pen(u, "hessian")
      c(rel(as.numeric(H %*% v), rich(function(x) pen(x, "gradient"), u, v,
                                       1e-4)),
        rel(D3(u), rich(function(x) pen(x, "hessian"), u, v, 1e-4)),
        rel(jfp_penalty_move(sp, des, k$cf_at(u), v, w = w),
            rich(D3, u, w, 1e-4)))
    }
    expect_lt(max(errs()), 1e-8)
    # negative control: without the block's parts through the filter
    ctl <- testthat::with_mocked_bindings(errs(),
      jfp_static_rows = function(st, dirs) lapply(st$params, function(p)
        matrix(0, st$spec@n_obs, st$joint$mf)),
      jfp_static_curv = function(st, g, dirs = list())
        matrix(0, st$joint$mf, st$joint$mf))
    expect_gt(min(ctl), 1e-6)
  }
})

test_that("a Jeffreys block inside the filter beside a moving block keeps the refusal", {
  set.seed(3)
  S <- 3L
  Tn <- 25L
  n <- S * Tn
  d <- data.frame(t = rep(1:Tn, S), id = factor(rep(1:S, each = Tn)),
                  xs = stats::runif(n, 0, 3), z = stats::rnorm(n))
  d$yn <- stats::rnbinom(n, mu = exp(0.5 + exp(-0.8 * d$xs)), size = 3)
  sp <- statmod_spec(yn ~ nl(~ b * exp(-r * xs), start = list(b = 1, r = 0.5)) +
                       gas(p = 1, q = 1, kappa1 ~ 1 + jfp(~ z), time = t,
                           by = id),
                     distributions7::negbin2_distrib(), d)
  des <- statmod_design(sp)
  expect_false(jfp_filter_outer_ok(sp, des, 1L))
  expect_false(jfp_filter_outer_ok(sp, des, 2L))
  # the same model without the moving block is admitted
  sp0 <- statmod_spec(yn ~ gas(p = 1, q = 1, kappa1 ~ 1 + jfp(~ z), time = t,
                               by = id),
                      distributions7::negbin2_distrib(), d)
  expect_true(jfp_filter_outer_ok(sp0, statmod_design(sp0), 1L))
})
