# Repairs found by the documentation audit of the outer criterion and the
# structural terms (statmodels7 0.199.6 with the tranche 1 repairs, b58f1b62).
# Every block fails on b58f1b62.

g1 <- distributions7::gaussian1_distrib()

# a Poisson smooth under ml(), and its ML criterion written out by hand: the
# coefficients at the penalized mode, the determinant over the penalty's range
ml_poisson <- function() {
  set.seed(2)
  d <- data.frame(x = stats::runif(150))
  d$cnt <- stats::rpois(150, exp(0.5 + sin(4 * d$x)))
  fit <- statmod(cnt ~ s(x, bspline_smooth(k = 8)),
                 distributions7::poisson_distrib(), d, outer_criterion = ml())
  X <- statmod_design(fit@spec)$mu$X
  D <- diag(c(0, 0, rep(1, 6)))
  y <- d$cnt
  V <- function(rho) {
    lam <- exp(rho)
    b <- fit@coefficients$mu
    for (it in 1:100) {
      m <- exp(drop(X %*% b))
      st <- solve(crossprod(X * m, X) + lam * D,
                  crossprod(X, y - m) - lam * D %*% b)
      b <- b + drop(st)
      if (max(abs(st)) < 1e-13) break
    }
    eta <- drop(X %*% b)
    H <- crossprod(X * exp(eta), X) + lam * D
    sum(y * eta - exp(eta)) - lam * sum(b * (D %*% b)) / 2 + 3 * rho -
      as.numeric(determinant(H[3:8, 3:8])$modulus) / 2
  }
  r0 <- log(fit@hyper$mu[[1L]][["lambda"]])
  h <- 2e-3
  list(fit = fit, ref = (V(r0 + h) - 2 * V(r0) + V(r0 - h)) / h^2)
}

test_that("the edf correction and the unconditional variance of an ml() fit read the ML criterion's curvature", {
  m <- ml_poisson()
  seen <- list()
  orig <- statmod_marginal_hess
  local_mocked_bindings(statmod_marginal_hess = function(...) {
    h <- orig(...)
    seen[[length(seen) + 1L]] <<- h
    h
  })
  summary(m$fit, correct = TRUE)
  suppressWarnings(vcov(m$fit, type = unconditional_vcov()))
  expect_gte(length(seen), 2L)
  for (h in seen) expect_equal(as.numeric(h), m$ref, tolerance = 1e-5)
})

test_that("the structural table reports standard errors beside a mixed covariance class", {
  set.seed(77)
  S <- matrix(c(0.36, 0, 0, 0.25), 2, 2)
  z <- matrix(stats::rnorm(16), 8, 2) %*% chol(S)
  g <- factor(rep(1:8, each = 30))
  y <- unlist(lapply(1:8, function(i) {
    f <- 0
    out <- numeric(30)
    for (t in 1:30) {
      mu <- z[i, 1L] + f
      out[t] <- stats::rnorm(1, mu, 1)
      f <- 0.25 * exp(z[i, 2L]) * (out[t] - mu) + 0.6 * f
    }
    out
  }))
  d <- data.frame(g = g, y = y)
  fit <- statmod(y ~ random(~ 1 | u | g) +
                   gas(p = 1, q = 1, by = g, kappa1 ~ 1 + random(~ 1 | u | g)),
                 g1, d, outer_criterion = NULL)
  tab <- statmod_structural_table(fit)
  free <- !grepl("held", tab$name) & tab$name != "omega"
  expect_true(all(is.finite(tab$se[free])))
  expect_true(all(tab$se[free] > 0))
})

test_that("a context without an approximation reads the default one", {
  expect_identical(ctx_approx(NULL), "opg")
})

test_that("statmod_breakpoint_test() rejects k and n_boot that are not whole", {
  set.seed(2)
  d <- data.frame(x = stats::runif(150))
  d$y <- 1 + 2 * pmax(d$x - 0.5, 0) + stats::rnorm(150, sd = 0.2)
  fb <- statmod(y ~ seg(x, psi = 0.4), g1, d)
  expect_error(statmod_breakpoint_test(fb, k = 3.7), "whole number of at least 3")
  expect_error(statmod_breakpoint_test(fb, n_boot = 2.5),
               "positive whole number")
})

test_that("marginal = \"all\" beside a structural term names what was written", {
  set.seed(3)
  d <- data.frame(t = 1:120, x = stats::rnorm(120))
  d$y <- stats::rnorm(120)
  err <- tryCatch(statmod(y ~ x + gas(p = 1, q = 1, time = t), g1, d,
                          outer_criterion = reml(marginal = "all")),
                  error = function(e) conditionMessage(e))
  expect_match(err, "marginal = \"all\" includes 'mu'", fixed = TRUE)
  expect_match(err, "marginal = \"none\"", fixed = TRUE)
  expect_no_match(err, "leave 'mu' out")
})

test_that("a label rejects two links of one class with different parameters", {
  mem <- data.frame(parameter = c("mu", "mu"), term = c("a", "b"),
                    name = c("g", "g"))
  expect_error(index_group(mem, list(linkfunctions7::power_link(lambda = 0.5),
                                     linkfunctions7::power_link(lambda = 2)),
                           c("L", "L")),
               "links differ")
  same <- index_group(mem, list(linkfunctions7::log_link(),
                                linkfunctions7::log_link()), c("L", "L"))
  expect_identical(nrow(same), 1L)
})

test_that("a frozen block on new data is read at its committed coefficients", {
  set.seed(3)
  x <- sort(stats::runif(120, 0, 10))
  d <- data.frame(x = x, y = ifelse(x < 6, 1, 3) + stats::rnorm(120, sd = 0.2))
  dnew <- data.frame(x = seq(0, 10, length.out = 25), y = 0)
  spec <- statmod_spec(y ~ jump(x, npsi = 1, psi = 6), g1, d)
  spec_nd <- statmod_respec(spec, dnew, need_response = FALSE)
  des <- statmod_design(spec_nd)
  tm <- attr(des, "state")$terms$mu[["jump(x, npsi = 1, psi = 6)"]]
  bc <- as.numeric(tm@blueprint$coef)
  c0 <- list(mu = c(0, bc))
  c1 <- list(mu = c(0, bc + c(0.5, 0)))
  a0 <- statmod_design_at(spec_nd, c0, des)$mu$adj
  a1 <- statmod_design_at(spec_nd, c1, des)$mu$adj
  expect_equal(a1, a0, tolerance = 1e-12)
})

test_that("a test whose mode error is missing says that it could not be read", {
  set.seed(1)
  d <- data.frame(x = stats::runif(80))
  d$y <- 1 + d$x + stats::rnorm(80)
  tt <- statmod_test(statmod(y ~ x, g1, d), param = "mu", coefname = "x",
                     test = lr_test())
  tt@mode_error <- NA_real_
  out <- paste(utils::capture.output(print(tt)), collapse = "\n")
  expect_match(out, "could not be read")
  expect_no_match(out, "stopped above its own mode")
})

test_that("a restart budget of zero proposes nothing", {
  set.seed(1)
  d <- data.frame(x = stats::runif(120, 0, 10))
  d$y <- 1 + 0.5 * d$x - 0.8 * pmax(d$x - 6, 0) + stats::rnorm(120, sd = 0.5)
  fit <- statmod(y ~ seg(x), g1, d)
  spec <- fit@spec
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  cfg <- inner_settings(fit@methods$smooth)
  res <- statmod_alternate(spec, design, blocks, fit@hyper, fit@methods$smooth,
                           unlist(fit@coefficients, use.names = FALSE),
                           cfg$expected, cfg$approx, cfg$maxit, cfg$tol,
                           verbosity(0))
  calls <- 0L
  local_mocked_bindings(seg_polish = function(...) {
    calls <<- calls + 1L
    stop("not expected")
  }, .package = "modelterms7")
  out <- statmod_boot_restart(spec, design, blocks, fit@hyper,
                              fit@methods$smooth, res, cfg$expected,
                              cfg$approx, cfg$maxit, cfg$tol, verbosity(0), 0L)
  expect_identical(calls, 0L)
  expect_identical(out$par, res$par)
})

test_that("the trace names the optimizer that the joint pass runs", {
  set.seed(4)
  d <- data.frame(t = 1:150, x = stats::rnorm(150))
  d$y <- 0.5 * d$x + stats::rnorm(150)
  out <- utils::capture.output(
    statmod(y ~ x + gas(p = 1, q = 1, time = t), g1, d,
            verbose = c(outer = FALSE, blocks = TRUE, inner = FALSE,
                        optimizer = FALSE)))
  hd <- grep("joint system", out, value = TRUE)
  expect_true(length(hd) > 0L)
  expect_true(all(grepl("[newton]", hd, fixed = TRUE)))
})

test_that("a failed edf correction is reported as not computed", {
  set.seed(1)
  d <- data.frame(x = stats::runif(100))
  d$y <- sin(3 * d$x) + stats::rnorm(100, sd = 0.3)
  fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), g1, d)
  local_mocked_bindings(statmod_edf_correction = function(...) stop("boom"))
  s <- summary(fit, correct = TRUE)
  nt <- paste(s@notes, collapse = "\n")
  expect_match(nt, "could not be")
  expect_no_match(nt, "estimated no hyperparameter here")
})

test_that("the edf correction reads the information that the fit used", {
  set.seed(3)
  d <- data.frame(x = stats::runif(300))
  d$y <- sin(4 * d$x) + stats::rt(300, df = 4) * 0.3
  fit <- statmod(y ~ s(x, bspline_smooth(k = 8)),
                 distributions7::pseudohuber_distrib(), d)
  expect_false(fit_expected(fit))
  ref <- statmod_edf_correction(fit@spec, fit@coefficients, fit@hyper,
                                statmod_design(fit@spec), fit@methods$outer,
                                expected = FALSE,
                                pinned = fit@methods$pinned)$total
  expect_equal(summary(fit, correct = TRUE)@df - summary(fit)@df, ref,
               tolerance = 1e-10)
})

test_that("frequency weights count as observations in bic() and nobs()", {
  set.seed(5)
  d <- data.frame(x = stats::runif(200))
  d$y <- sin(4 * d$x) + stats::rnorm(200, sd = 0.4)
  fw <- statmod(y ~ s(x, bspline_smooth(k = 10)), g1, d,
                weights = rep(2, 200), outer_criterion = bic())
  fs <- statmod(y ~ s(x, bspline_smooth(k = 10)), g1, rbind(d, d),
                outer_criterion = bic())
  expect_equal(unlist(fw@hyper), unlist(fs@hyper), tolerance = 1e-5)
  expect_equal(nobs(fw), 400)
  expect_equal(BIC(fw), BIC(fs), tolerance = 1e-6)
  expect_equal(summary(fw)@bic, summary(fs)@bic, tolerance = 1e-6)
})

# a smooth and a break-point in mu, sigma ~ w, and the default marginal, which
# estimates sigma's and nu's coefficients on the criterion; at the mode the
# break-point sits on an observation, where the objective has a kink
kink_fit <- function(hessian) {
  set.seed(5)
  d <- data.frame(x = stats::runif(250), z = stats::runif(250),
                  w = stats::runif(250))
  d$y <- sin(4 * d$x) + pmax(d$z - 0.5, 0) * 2 + 0.4 * stats::rt(250, 5)
  g <- distributions7::student_t1_distrib()
  fit <- statmod(y ~ s(x, bspline_smooth(k = 8)) + seg(z, psi = 0.5) |
                   sigma ~ w, g, d, outer_criterion = reml(hessian = hessian))
  sp <- fit@spec
  ds <- statmod_design(sp)
  bl <- statmod_blocks(sp, ds)
  idx <- outer_hyper_index(sp, bl)
  m <- fit@methods$outer
  mc <- fit_marginal_context(sp, ds, fit@coefficients, m, fit@methods$pinned)
  basis <- integrated_basis(mc$spec, mc$design, m@kind, gamma = TRUE)
  inner <- iwls_resolve(iwls(), g)
  inner@tol <- 1e-11
  inner@criterion <- NULL
  cfg <- inner_settings(inner)
  beta0 <- unlist(fit@coefficients[g@params], use.names = FALSE)
  nh <- nrow(idx)
  ng <- length(mc$gam$where)
  v0 <- c(hyper_to_eta(fit@hyper, idx), beta0[mc$gam$where])
  mode <- function(v) {
    hy <- eta_to_hyper(v[seq_len(nh)], idx, fit@hyper)
    spv <- statmod_hold(mc$spec, mc$gam, v[nh + seq_len(ng)])
    b0 <- beta0
    b0[mc$gam$where] <- v[nh + seq_len(ng)]
    r <- statmod_alternate(spv, mc$design, bl, hy, inner, b0, cfg$expected,
                           cfg$approx, cfg$maxit, cfg$tol, verbosity(0),
                           hold_refresh = TRUE)
    list(sp = spv, cf = r$obj$split(r$par), hy = hy)
  }
  value <- function(v) {
    r <- mode(v)
    statmod_marginal(r$sp, mc$design, r$cf, r$hy, m, "opg", basis)$value
  }
  grad <- function(v) {
    r <- mode(v)
    statmod_marginal_grad(r$sp, mc$design, r$cf, r$hy, m, idx, basis,
                          gam = mc$gam)
  }
  step <- function(j, h) h * if (j <= nh) 1 else max(1, abs(v0[[j]]))
  list(fit = fit, mc = mc, idx = idx, basis = basis, v0 = v0, nv = nh + ng,
       value = value, grad = grad, step = step, m = m)
}

test_that("the outer gradient holds a break-point that a kink holds at the mode", {
  k <- kink_fit("observed")
  expect_true(k$fit@converged)
  v <- k$v0
  v[1L] <- v[1L] + 0.3
  g <- k$grad(v)
  fd <- vapply(seq_len(k$nv), function(j) {
    s <- k$step(j, 1e-4)
    e <- replace(numeric(k$nv), j, s)
    (k$value(v + e) - k$value(v - e)) / (2 * s)
  }, numeric(1))
  expect_equal(g, fd, tolerance = 5e-4)
})

test_that("the expected route's Hessian spans the coefficients that the criterion estimates", {
  k <- kink_fit("expected")
  H <- statmod_marginal_hess(k$mc$spec, k$mc$design, k$fit@coefficients,
                             k$fit@hyper, k$m, k$idx, k$basis, gam = k$mc$gam)
  expect_identical(dim(H), c(k$nv, k$nv))
  f0 <- k$value(k$v0)
  d2 <- vapply(seq_len(k$nv), function(j) {
    s <- k$step(j, 4e-3)
    e <- replace(numeric(k$nv), j, s)
    (k$value(k$v0 + e) - 2 * f0 + k$value(k$v0 - e)) / s^2
  }, numeric(1))
  expect_equal(diag(H), d2, tolerance = 1e-4)
})

test_that("a break-point developed over a covariate slides along the hyperplane of its kink", {
  set.seed(5)
  n <- 150
  d <- data.frame(x = stats::runif(n), z = stats::runif(n, 0, 10),
                  w = stats::runif(n), v = stats::runif(n))
  d$y <- sin(4 * d$x) + 0.5 * d$z + 1.5 * pmax(d$z - (3 + 4 * d$w), 0) +
    stats::rnorm(n, sd = exp(-0.5 + 0.5 * d$v))
  fit <- statmod(y ~ s(x, bspline_smooth(k = 8)) + seg(z, psi ~ w) |
                   sigma ~ v, g1, d)
  ds <- statmod_design(fit@spec)
  k <- grep("psi", ds$mu$coef_names)
  kap <- fit@coefficients$mu[k]
  j <- which.min(abs(d$z - (kap[1L] + kap[2L] * d$w)))
  # the mode is on the kink of observation j
  expect_lt(abs(d$z[j] - (kap[1L] + kap[2L] * d$w[j])), 1e-6)
  # the inner fit started on that kink but away from the constrained optimum
  # along its hyperplane returns to it, rather than holding the two
  # coordinates where it started
  sp <- fit@spec
  bl <- statmod_blocks(sp, ds)
  b0 <- unlist(fit@coefficients[g1@params], use.names = FALSE)
  dirn <- c(-d$w[j], 1)
  start <- b0
  start[k] <- start[k] + 0.05 * dirn
  cfg <- inner_settings(fit@methods$smooth)
  r <- statmod_alternate(sp, ds, bl, fit@hyper, fit@methods$smooth, start,
                         cfg$expected, cfg$approx, cfg$maxit, cfg$tol,
                         verbosity(0))
  obj <- statmod_objective(sp, fit@hyper, ds, FALSE, "opg")
  f <- function(t) {
    b <- r$par
    b[k] <- b[k] + t * dirn
    obj$fn(b)
  }
  h <- 1e-4
  slope <- (f(h) - f(-h)) / (2 * h)
  curv <- (f(h) - 2 * f(0) + f(-h)) / h^2
  expect_lt(abs(slope / curv), 1e-6)
  expect_lte(obj$fn(r$par), obj$fn(b0) + 1e-8 * abs(obj$fn(b0)))
})
