# Several gas() terms in different equations, filtered as one recursion: the
# joining, the fit, every derivative the fit and the outer criterion read,
# prediction and simulation.

# a gaussian whose mean and log-scale are both driven by their scores
jt_sim <- function(n = 120, seed = 5) {
  set.seed(seed)
  x <- stats::runif(n)
  f <- c(0, 0)
  y <- numeric(n)
  for (t in seq_len(n)) {
    mu <- 1 + 0.5 * x[t] + f[1]
    s <- exp(-0.2 + f[2])
    y[t] <- stats::rnorm(1, mu, s)
    r <- y[t] - mu
    f <- c(0.3 * r / s^2 + 0.6 * f[1], 0.1 * (r^2 / s^2 - 1) + 0.7 * f[2])
  }
  data.frame(t = seq_len(n), x = x, y = y)
}

# the exact gradient over the coefficients and the filter's parameters, and
# the joint information, against numDeriv on the log-likelihood at a point
# away from the optimum
jt_derivs <- function(form, distrib, data, u0) {
  spec <- statmod_spec(form, distrib, data)
  des <- statmod_design(spec)
  sst <- statmod_structural_state(des)
  nmz <- names(sst$zeta[[1L]])
  params <- spec@distrib@params
  npar <- vapply(des, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  nb <- sum(npar)
  np <- length(nmz)
  split_u <- function(u) list(
    coef = stats::setNames(lapply(seq_along(params), function(a)
      u[offs[a] + seq_len(npar[a])]), params),
    zeta = u[nb + seq_len(np)])
  setz <- function(zz) {
    sst$zeta[[1L]] <- stats::setNames(zz, nmz)
    sst$key <- NULL
  }
  L <- function(u) {
    s <- split_u(u)
    setz(s$zeta)
    statmod_loglik_at(spec, s$coef, des)
  }
  G <- function(u) {
    s <- split_u(u)
    setz(s$zeta)
    c(unlist(statmod_score_at(spec, s$coef, des), use.names = FALSE),
      statmod_structural_score(spec, s$coef, des)[[1L]])
  }
  g_an <- G(u0)
  g_num <- numDeriv::grad(L, u0)
  setz(split_u(u0)$zeta)
  I_an <- statmod_full_information(spec, split_u(u0)$coef, des)
  I_num <- -numDeriv::jacobian(G, u0)
  held <- sst$held[[1L]]
  keep <- setdiff(seq_along(u0), nb + match(held, nmz))
  c(gradient = max(abs(g_an[keep] - g_num[keep])) /
      max(1, max(abs(g_num[keep]))),
    information = max(abs(I_an - I_num[keep, keep])) /
      max(abs(I_num[keep, keep])))
}

test_that("gas() terms in several equations are joined into one structural term", {
  dd <- jt_sim()
  spec <- statmod_spec(y ~ x + gas(p = 1, q = 1, time = t) |
                         sigma ~ gas(p = 2, q = 1, time = t),
                       distributions7::gaussian1_distrib(), dd)
  su <- statmod_structural(spec)
  expect_length(su, 1L)
  expect_identical(su[[1L]]$term, "gas_joint")
  tm <- spec@terms$mu$gas_joint
  expect_true(S7::S7_inherits(tm, modelterms7::GasJoint))
  expect_identical(tm@params, c("mu", "sigma"))
  expect_identical(unname(tm@keys), c("gas(p = 1, q = 1, time = t)",
                                      "gas(p = 2, q = 1, time = t)"))
  expect_null(spec@terms$sigma[["gas(p = 2, q = 1, time = t)"]])
  # one gas() term is left as it is
  s1 <- statmod_spec(y ~ x + gas(p = 1, q = 1, time = t),
                     distributions7::gaussian1_distrib(), dd)
  expect_true(S7::S7_inherits(s1@terms$mu[["gas(p = 1, q = 1, time = t)"]],
                              modelterms7::GasTerm))
})

test_that("a joint fit maximizes the log-likelihood of the recursion written out", {
  dd <- jt_sim()
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t, p_cross = c(sigma = 1)) |
                   sigma ~ gas(p = 1, q = 1, time = t, p_cross = c(mu = 1)),
                 distributions7::gaussian1_distrib(), dd)
  cf <- coef(fit)
  ll_of <- function(v) {
    f <- c(0, 0)
    ll <- 0
    for (t in seq_len(nrow(dd))) {
      mu <- v[1] + v[2] * dd$x[t] + f[1]
      ls <- v[3] + f[2]
      s2 <- exp(2 * ls)
      ll <- ll + stats::dnorm(dd$y[t], mu, sqrt(s2), log = TRUE)
      r <- dd$y[t] - mu
      u <- c(r / s2, r^2 / s2 - 1)
      f <- c(v[4] * u[1] + v[5] * u[2] + v[6] * f[1],
             v[7] * u[2] + v[8] * u[1] + v[9] * f[2])
    }
    ll
  }
  v0 <- c(cf$mu[["(Intercept)"]], cf$mu[["x"]], cf$sigma[["(Intercept)"]],
          cf$mu[["gas.kappa1"]], cf$mu[["gas.kappa1_sigma"]],
          cf$mu[["gas.xi1"]], cf$sigma[["gas.kappa1"]],
          cf$sigma[["gas.kappa1_mu"]], cf$sigma[["gas.xi1"]])
  expect_equal(ll_of(v0), as.numeric(logLik(fit)), tolerance = 1e-10)
  # the fit is a maximum of that function: its gradient there is nil, the
  # loadings on their log scale and the persistences on atanh
  w <- c(v0[1:3], log(v0[4]), v0[5], atanh(v0[6]), log(v0[7]), v0[8],
         atanh(v0[9]))
  back <- function(w) c(w[1:3], exp(w[4]), w[5], tanh(w[6]), exp(w[7]), w[8],
                        tanh(w[9]))
  g <- numDeriv::grad(function(w) ll_of(back(w)), w)
  expect_lt(max(abs(g)), 1e-4)
})

test_that("the gradient and the information of a joint filter are exact", {
  skip_if_not_installed("numDeriv")
  dd <- jt_sim()
  g <- jt_derivs(y ~ x + gas(p = 1, q = 1, time = t) |
                   sigma ~ gas(p = 1, q = 1, time = t),
                 distributions7::gaussian1_distrib(), dd,
                 c(0.9, 0.4, -0.1, 0, log(0.2), 0.8, 0, log(0.08), 1.1))
  expect_lt(max(g), 1e-7)
  g <- jt_derivs(y ~ x + gas(p = 1, q = 2, time = t, p_cross = c(sigma = 1)) |
                   sigma ~ gas(p = 2, q = 1, time = t, p_cross = c(mu = 2)),
                 distributions7::gaussian1_distrib(), dd,
                 c(0.9, 0.4, -0.1, 0, log(0.2), 0.8, 0.2, 0.05,
                   0, log(0.08), log(0.03), 1.1, -0.1, 0.07))
  expect_lt(max(g), 1e-7)
  # one way only, at lags beyond the term's own p: mu loads sigma's score at
  # the lags 1..3, sigma loads its own score only
  g <- jt_derivs(y ~ x + gas(p = 1, q = 1, time = t,
                             p_cross = c(sigma = 3)) |
                   sigma ~ gas(p = 1, q = 1, time = t),
                 distributions7::gaussian1_distrib(), dd,
                 c(0.9, 0.4, -0.1, 0, log(0.2), 0.8, 0.05, -0.03, 0.02,
                   0, log(0.08), 1.1))
  expect_lt(max(g), 1e-7)
  # scaled jointly: the information of a Student t couples sigma and nu
  set.seed(3)
  dt <- data.frame(t = 1:70, y = stats::rt(70, df = 5) * 1.2)
  for (d in c(0.5, 1)) {
    g <- jt_derivs(y ~ 1 | sigma ~ gas(p = 1, q = 1, time = t, scaling = d) |
                     nu ~ gas(p = 1, q = 1, time = t, scaling = d),
                   distributions7::student_t1_distrib(), dt,
                   c(0.1, 0.2, 1.6, 0, log(0.1), 0.7, 0, log(0.1), 0.6))
    expect_lt(max(g), 1e-6, label = sprintf("scaling %g", d))
  }
})

test_that("the derivative rows of a joint filter and of its continuation are exact", {
  skip_if_not_installed("numDeriv")
  dd <- jt_sim(80)
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t, p_cross = c(sigma = 1)) |
                   sigma ~ gas(p = 1, q = 1, time = t, p_cross = c(mu = 1)),
                 distributions7::gaussian1_distrib(), dd)
  spec <- fit@spec
  des <- statmod_design(spec)
  sst <- statmod_structural_state(des)
  cf <- fit@coefficients
  z0 <- sst$zeta$gas_joint
  free <- setdiff(names(z0), sst$held$gas_joint)
  u0 <- c(unlist(cf[c("mu", "sigma")], use.names = FALSE), z0[free])
  nb <- length(unlist(cf))
  at <- function(u) {
    co <- list(mu = u[1:2], sigma = u[3])
    z <- z0
    z[free] <- u[nb + seq_along(free)]
    sst$zeta$gas_joint <- z
    sst$key <- NULL
    co
  }
  fj <- joint_filter_jacobian(spec, des, cf)
  for (p in c("mu", "sigma")) {
    num <- numDeriv::jacobian(function(u) statmod_eta(spec, des, at(u))$eta[[p]],
                              u0)
    at(u0)
    expect_equal(unname(fj$V[[p]]), num, tolerance = 1e-6, info = p)
  }
  # past the series: the continued predictor's rows against a difference of
  # the continuation itself
  nd <- data.frame(t = 81:84, x = c(0.2, 0.4, 0.6, 0.8))
  ospec <- statmod_respec(spec, nd, need_response = FALSE)
  odes <- statmod_design(ospec)
  cont <- statmod_eta_continued(fit, ospec, odes, deriv = TRUE)
  # the continuation reads the filter's parameters from the fitted
  # specification, so a perturbed copy of the fit carries them there
  fit_at <- function(u) {
    co <- list(mu = u[1:2], sigma = u[3])
    z <- z0
    z[free] <- u[nb + seq_along(free)]
    sp <- fit@spec
    stt <- sp@structural
    stt$gas_joint <- z
    sp@structural <- stt
    fit2 <- fit
    fit2@coefficients <- co
    fit2@spec <- sp
    statmod_eta_continued(fit2, ospec, odes)$eta
  }
  for (p in c("mu", "sigma")) {
    num <- numDeriv::jacobian(function(u) fit_at(u)[[p]], u0)
    expect_equal(unname(cont$cont_cols[[p]]$X), num, tolerance = 1e-6,
                 info = paste("continued", p))
  }
})

test_that("the outer gradient and Hessian are exact beside a joint filter", {
  skip_on_cran()
  # a curved effect, so that the smoothing parameter stays inside its range:
  # at the edge, where a linear effect sends it, the criterion is flat to
  # 1e-6 and a difference of the gradient cannot judge the Hessian
  dd <- jt_sim(100)
  dd$y <- dd$y + sin(2 * pi * dd$x)
  fit <- statmod(y ~ s(x, pspline_smooth(k = 6)) + gas(p = 1, q = 1, time = t) |
                   sigma ~ gas(p = 1, q = 1, time = t),
                 distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  spec <- fit@spec
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  method <- fit@methods$outer
  idx <- outer_hyper_index(spec, blocks)
  basis <- integrated_basis(spec, design, method@kind)
  inner <- iwls()
  cfg <- inner_settings(inner, spec@distrib)
  beta0 <- unlist(fit@coefficients[spec@distrib@params], use.names = FALSE)
  sst <- statmod_structural_state(design)
  z0 <- sst$zeta
  refit <- function(eta) {
    sst$zeta <- z0
    hy <- eta_to_hyper(eta, idx, fit@hyper)
    r <- statmod_alternate(spec, design, blocks, hy, inner, beta0,
                           cfg$expected, cfg$approx, cfg$maxit, cfg$tol,
                           verbosity(0), hold_refresh = TRUE)
    list(cf = r$obj$split(r$par), hy = hy)
  }
  fn <- function(eta) {
    r <- refit(eta)
    statmod_marginal(spec, design, r$cf, r$hy, method, basis = basis)$value
  }
  gr <- function(eta) {
    r <- refit(eta)
    statmod_marginal_grad(spec, design, r$cf, r$hy, method, idx, basis)
  }
  eta0 <- hyper_to_eta(fit@hyper, idx)
  # no penalty covers the filter's own parameters, so the criterion's
  # determinant is over the coefficients and its Hessian is the one the fit
  # reads, statmod_marginal_hess()
  A <- statmod_marginal_hess(spec, design, fit@coefficients, fit@hyper,
                             method, idx, basis)
  expect_false(is.null(A))
  eta1 <- eta0 - 0.5
  cd <- function(f, eta, h) (f(eta + h) - f(eta - h)) / (2 * h)
  expect_equal(gr(eta1), cd(fn, eta1, 1e-3), tolerance = 1e-4)
  expect_equal(A[1, 1], cd(gr, eta0, 1e-3), tolerance = 1e-4)

  # the pieces of the assembly on the joint vector, which the criterion
  # reads where a penalty covers the filter's parameters, against
  # differences of tr(M K(u)) with M held
  coef <- fit@coefficients
  hyper <- fit@hyper
  jd <- joint_design_rows(spec, design, coef)
  K <- statmod_marginal_full(spec, design, coef, hyper, NULL)
  M <- chol2inv(chol(K))
  key <- jd$f$term
  z0 <- sst$zeta[[key]]
  free <- setdiff(names(z0), sst$held[[key]])
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  nb <- sum(npar)
  u0 <- c(unlist(coef[params], use.names = FALSE), z0[free])
  setu <- function(u) {
    co <- stats::setNames(lapply(seq_along(params), function(a)
      u[offs[a] + seq_len(npar[a])]), params)
    z <- z0
    z[free] <- u[nb + seq_along(free)]
    sst$zeta[[key]] <- z
    sst$key <- NULL
    co
  }
  Kat <- function(u) as.matrix(statmod_marginal_full(spec, design, setu(u),
                                                     hyper, NULL))
  setu(u0)
  st <- joint_grad_parts(spec, design, coef, jd, M, order = 4L)
  set.seed(1)
  v <- stats::rnorm(length(u0)) * 0.1
  w <- stats::rnorm(length(u0)) * 0.1
  dk <- joint_dk_matrix(spec, design, jd, st, v)
  h <- 1e-5
  dKfd <- (Kat(u0 + h * v) - Kat(u0 - h * v)) / (2 * h)
  expect_lt(max(abs(dk - dKfd)) / max(abs(dKfd)), 1e-7)
  setu(u0)
  first <- sum(st$u * v) + joint_chain_extra(spec, design, jd, M, st, v)
  expect_equal(first, sum(M * dk), tolerance = 1e-10)
  second <- joint_chain_extra2(spec, design, jd, M, st, v, w)
  h2 <- 1e-3
  trK <- function(u) sum(M * Kat(u))
  fd2 <- (trK(u0 + h2 * v + h2 * w) - trK(u0 + h2 * v - h2 * w) -
            trK(u0 - h2 * v + h2 * w) + trK(u0 - h2 * v - h2 * w)) / (4 * h2^2)
  setu(u0)
  expect_equal(second, fd2, tolerance = 1e-5)
})

test_that("rstatmod draws from the joint recursion", {
  dd <- jt_sim(60)
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t) |
                   sigma ~ gas(p = 1, q = 1, time = t),
                 distributions7::gaussian1_distrib(), dd)
  set.seed(9)
  sim <- simulate(fit, nsim = 1)
  expect_length(sim[[1L]], 60L)
  # the drawn response, filtered at the fitted parameters, reproduces the
  # levels it was drawn at
  set.seed(9)
  spec <- fit@spec
  des <- statmod_design(spec)
  ev <- statmod_eta(spec, des, fit@coefficients)
  r <- rstatmod_eta(spec, des, fit@coefficients)
  d2 <- dd
  d2$y <- r$y
  s2 <- statmod_spec(y ~ x + gas(p = 1, q = 1, time = t) |
                       sigma ~ gas(p = 1, q = 1, time = t),
                     distributions7::gaussian1_distrib(), d2)
  ds2 <- statmod_design(s2)
  sst2 <- statmod_structural_state(ds2)
  sst2$zeta <- statmod_structural_state(des)$zeta
  e2 <- statmod_eta(s2, ds2, fit@coefficients)
  expect_equal(e2$eta$mu, r$eta$mu, tolerance = 1e-12)
  expect_equal(e2$eta$sigma, r$eta$sigma, tolerance = 1e-12)
})

test_that("joint filters that cannot be fitted are rejected by name", {
  dd <- jt_sim(40)
  g <- distributions7::gaussian1_distrib()
  expect_error(statmod_spec(y ~ gas(p = 1, q = 1, time = t) |
                              sigma ~ gas(p = 1, q = 1, time = t, scaling = 1),
                            g, dd), "same 'scaling'")
  # a score from an equation without a gas() term, with one term or two
  expect_error(statmod_spec(y ~ gas(p = 1, q = 1, time = t,
                                    p_cross = c(sigma = 1)), g, dd),
               "'sigma' has no gas\\(\\) term")
  expect_error(statmod_spec(y ~ gas(p = 1, q = 1, time = t,
                                    p_cross = c(nu = 1)) |
                              sigma ~ gas(p = 1, q = 1, time = t), g, dd),
               "'nu' has no gas\\(\\) term")
  expect_error(statmod_spec(y ~ gas(p = 1, q = 1, time = t) +
                              gas(p = 2, q = 1, time = t, label = "g2"), g, dd),
               "carries two gas\\(\\) terms")
})

test_that("a joint filter carries subformulas on its terms' parameters", {
  skip_on_cran()
  set.seed(31)
  m <- 6L
  ni <- 30L
  u <- stats::rnorm(m, sd = 0.6)
  g <- factor(rep(seq_len(m), each = ni))
  y <- unlist(lapply(seq_len(m), function(i) {
    f <- 0
    out <- numeric(ni)
    for (t in seq_len(ni)) {
      mu <- 0.4 + f
      out[t] <- stats::rnorm(1, mu, 1)
      f <- 0.25 * exp(u[i]) * (out[t] - mu) + 0.6 * f
    }
    out
  }))
  dp <- data.frame(y = y, g = g, t = rep(seq_len(ni), m))
  fit <- function(f, ...) statmod(f, distributions7::gaussian1_distrib(), dp, ...)
  plain <- fit(y ~ gas(p = 1, q = 1, by = g, time = t) |
                 sigma ~ gas(p = 1, q = 1, by = g, time = t))
  one <- fit(y ~ gas(p = 1, q = 1, by = g, time = t, kappa1 ~ 1) |
               sigma ~ gas(p = 1, q = 1, by = g, time = t))
  # a subformula of one is the scalar loading on its log chart
  expect_equal(as.numeric(logLik(one)), as.numeric(logLik(plain)),
               tolerance = 1e-8)
  nd <- data.frame(g = factor(rep(1:2, each = 2), levels = levels(dp$g)),
                   t = rep(ni + 1:2, 2))
  expect_equal(predict(one, what = "mu", newdata = nd),
               predict(plain, what = "mu", newdata = nd), tolerance = 1e-6)
  dev <- fit(y ~ gas(p = 1, q = 1, by = g, time = t,
                     kappa1 ~ 1 + random(~ 1 | g)) |
               sigma ~ gas(p = 1, q = 1, by = g, time = t),
             outer_criterion = reml())
  expect_true(is.finite(as.numeric(logLik(dev))))
  # each equation reports its own term's components, named and written as in
  # that term's formula: the development of mu's loading appears under mu and
  # not under sigma
  out <- utils::capture.output(print(summary(dev)))
  at_sigma <- grep("^=== sigma", out)
  hdr <- "^  kappa1 ~ 1 \\+ random\\(~1 \\| g\\)$"
  expect_length(grep(hdr, out[seq_len(at_sigma)]), 1L)
  expect_length(grep(hdr, out[-seq_len(at_sigma)]), 0L)
  expect_length(grep("mu\\.kappa1", out), 0L)
  # a coordinate of a development that the fit holds is marked as held: the
  # intercept of omega ~ x is aliased with sigma's own intercept
  dx <- dp
  dx$x <- stats::rnorm(nrow(dx))
  dsx <- statmod(y ~ gas(p = 1, q = 1, by = g, time = t) |
                   sigma ~ gas(p = 1, q = 1, by = g, time = t, omega ~ x),
                 distributions7::gaussian1_distrib(), dx)
  ox <- utils::capture.output(print(summary(dsx)))
  expect_length(grep("^  omega ~ x$", ox), 1L)
  expect_length(grep("^    \\(Intercept\\) \\(held\\)", ox), 1L)
  expect_true(all(is.finite(predict(dev, what = "mu", newdata = nd))))
})

# a panel whose leverage loading (the signed residual's loading in log sigma)
# varies by group
lev_panel <- function(m, n, seed = 5) {
  set.seed(seed)
  kg <- stats::rnorm(m, -0.06, 0.06)
  y <- numeric(m * n)
  for (i in seq_len(m)) {
    f1 <- 0; f2 <- 0
    for (t in seq_len(n)) {
      r <- (i - 1L) * n + t
      mu <- 0.2 + f1
      ls <- f2
      y[r] <- stats::rnorm(1, mu, exp(ls))
      s1 <- (y[r] - mu) / exp(2 * ls)
      s2 <- (y[r] - mu)^2 / exp(2 * ls) - 1
      f1 <- 0.2 * s1 * exp(2 * ls) + 0.5 * f1
      f2 <- 0.08 * s2 + kg[i] * s1 * exp(2 * ls) + 0.85 * f2
    }
  }
  data.frame(y = y, g = factor(rep(seq_len(m), each = n)),
             t = rep(seq_len(n), m))
}

test_that("a cross loading carries a subformula, reported under its equation", {
  skip_on_cran()
  dl <- lev_panel(6L, 150L)
  fit <- function(f, ...) statmod(f, distributions7::gaussian1_distrib(), dl, ...)
  plain <- fit(y ~ gas(p = 1, q = 1, by = g, time = t, p_cross = c(sigma = 1)) |
                 sigma ~ gas(p = 1, q = 1, by = g, time = t, p_cross = c(mu = 1)))
  one <- fit(y ~ gas(p = 1, q = 1, by = g, time = t, p_cross = c(sigma = 1)) |
               sigma ~ gas(p = 1, q = 1, by = g, time = t, p_cross = c(mu = 1),
                           kappa1_mu ~ 1))
  expect_equal(as.numeric(logLik(one)), as.numeric(logLik(plain)),
               tolerance = 1e-8)
  dev <- fit(y ~ gas(p = 1, q = 1, by = g, time = t, p_cross = c(sigma = 1)) |
               sigma ~ gas(p = 1, q = 1, by = g, time = t, p_cross = c(mu = 1),
                           kappa1_mu ~ 1 + random(~ 1 | g)),
             outer_criterion = reml())
  out <- utils::capture.output(print(summary(dev)))
  at_sigma <- grep("^=== sigma", out)
  hdr <- which(out == "  kappa1_mu ~ 1 + random(~1 | g)")
  expect_length(hdr, 1L)
  expect_gt(hdr, at_sigma)
  # the random effect's standard deviation is shown with it, although the
  # joint filter keeps its hyperparameters under mu
  sd_row <- startsWith(out, "      sigma [reml]")
  expect_identical(sum(sd_row[-seq_len(at_sigma)]), 1L)
  expect_identical(sum(sd_row[seq_len(at_sigma)]), 0L)
})

test_that("a random effect in the second term of a joint filter shows its scale", {
  skip_on_cran()
  dl <- lev_panel(6L, 100L, seed = 7)
  fit <- statmod(y ~ gas(p = 1, q = 1, by = g, time = t) |
                   sigma ~ gas(p = 1, q = 1, by = g, time = t,
                               kappa1 ~ 1 + random(~ 1 | g)),
                 distributions7::gaussian1_distrib(), dl,
                 outer_criterion = reml())
  out <- utils::capture.output(print(summary(fit)))
  at_sigma <- grep("^=== sigma", out)
  sd_row <- startsWith(out, "      sigma [reml]")
  expect_identical(sum(sd_row[-seq_len(at_sigma)]), 1L)
})
