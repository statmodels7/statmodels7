# gas(scaling = d): the filter driven by the score times the expected
# information to the power -d. What has to hold is that every derivative the
# fit and the outer criterion read is the derivative of THAT recursion -- the
# score's derivatives replaced by those of u = s I^-d wherever the recursion
# reads them, and left the log-likelihood's wherever they weigh it -- and that
# d = 0 is the recursion the package ran before.

sc_pois <- function(n = 150, seed = 4, d = 0.5) {
  set.seed(seed)
  x <- stats::rnorm(n)
  z <- stats::runif(n)
  f <- 0
  y <- integer(n)
  for (t in seq_len(n)) {
    lam <- exp(0.7 + 0.3 * x[t] + f)
    y[t] <- stats::rpois(1, lam)
    f <- 0.15 * (y[t] - lam) * lam^(-d) + 0.8 * f
  }
  data.frame(t = seq_len(n), x = x, z = z, y = y)
}

# the exact gradient over the coefficients and the term's parameters, and the
# joint information, against numDeriv on the log-likelihood at a point AWAY
# from the optimum
sc_derivs <- function(form, distrib, data, u0) {
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
  rows <- if (length(g_an) == length(u0)) keep else seq_along(g_an)
  gk <- if (length(g_an) == length(u0)) g_an[keep] else g_an
  c(gradient = max(abs(gk - g_num[keep])) / max(1, max(abs(g_num[keep]))),
    information = max(abs(I_an - I_num[rows, keep])) /
      max(abs(I_num[rows, keep])))
}

test_that("gas() takes a scaling and rejects anything but one finite number", {
  expect_identical(modelterms7::gas()@scaling, 0)
  expect_identical(modelterms7::gas(scaling = 0.5)@scaling, 0.5)
  expect_error(modelterms7::gas(scaling = "a"), "single finite number")
  expect_error(modelterms7::gas(scaling = NA), "single finite number")
  expect_error(modelterms7::gas(scaling = c(0, 1)), "single finite number")
})

test_that("at scaling 0 the driving arrays are the family's own", {
  dd <- sc_pois()
  spec <- statmod_spec(y ~ x + gas(p = 1, q = 1, time = t), poisson_distrib(), dd)
  th <- list(mu = rep(2, nrow(dd)))
  gl <- distributions7::distrib_gradient(spec@distrib, dd$y, th, scale = "link")
  H <- distributions7::distrib_hessian(spec@distrib, dd$y, th, scale = "link")
  dr <- filter_driving(spec, th, 1L, 0, gl, H)
  expect_identical(dr$gl, gl)
  expect_identical(dr$H, H)
})

test_that("the gradient and the information are exact for a scaled score", {
  skip_if_not_installed("numDeriv")
  dd <- sc_pois()
  # the information of a Poisson mean moves with its own predictor
  for (d in c(0.5, 1, -0.7)) {
    g <- sc_derivs(y ~ x + gas(p = 1, q = 1, time = t, scaling = d),
                   poisson_distrib(), dd, c(0.6, 0.25, 0, log(0.12), 1.2))
    expect_lt(max(g), 1e-8, label = sprintf("poisson, d = %g", d))
  }
  # and with the size of a negative binomial, which carries covariates of
  # its own, so the scaling couples the two equations
  g <- sc_derivs(y ~ x + gas(p = 1, q = 1, time = t, scaling = 0.5) | theta ~ z,
                 negbin2_distrib(), dd, c(0.6, 0.25, 1.0, 0.3, 0, log(0.12), 1.2))
  expect_lt(max(g), 1e-8)
  # a gaussian mean: the information is 1/sigma^2, so the scaled score reads
  # the scale's own predictor
  set.seed(8)
  dg <- dd
  dg$y <- dd$y + stats::rnorm(nrow(dd))
  g <- sc_derivs(y ~ x + gas(p = 1, q = 1, time = t, scaling = 1) - 1 | sigma ~ z,
                 gaussian1_distrib(), dg, c(0.3, 0.1, -0.2, 0.2, log(0.3), 1.0))
  expect_lt(max(g), 1e-8)
})

test_that("a scaled Poisson filter on the identity link is INGARCH(1, 1)", {
  disc <- data.frame(year = 1860:1959, n = as.numeric(discoveries))
  fit <- statmod(n ~ gas(p = 1, q = 1, time = year, scaling = 1),
                 distrib = poisson_distrib(link_mu = linkfunctions7::identity_link()),
                 data = disc)
  cf <- coef(fit)$mu
  lam <- as.numeric(predict(fit, what = "mu"))
  # lambda_{t+1} = c (1 - b) + a y_t + (b - a) lambda_t, the intercept c
  # being the level the filter fluctuates around
  c0 <- cf[["(Intercept)"]]
  a <- cf[["gas.alpha1"]]
  b <- cf[["gas.phi1"]]
  rec <- numeric(nrow(disc))
  rec[1] <- c0
  for (t in 2:nrow(disc)) {
    rec[t] <- c0 * (1 - b) + a * disc$n[t - 1] + (b - a) * rec[t - 1]
  }
  expect_lt(max(abs(rec - lam)), 1e-8)
})

test_that("a constant information only rescales the loading", {
  nile <- data.frame(year = 1871:1970, flow = as.numeric(Nile))
  f0 <- statmod(flow ~ gas(p = 1, q = 1, time = year),
                distrib = gaussian1_distrib(), data = nile)
  f1 <- statmod(flow ~ gas(p = 1, q = 1, time = year, scaling = 1),
                distrib = gaussian1_distrib(), data = nile)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f0)), tolerance = 1e-6)
  s2 <- exp(2 * coef(f0)$sigma[["(Intercept)"]])
  expect_equal(coef(f1)$mu[["gas.alpha1"]], coef(f0)$mu[["gas.alpha1"]] / s2,
               tolerance = 1e-4)
  expect_equal(coef(f1)$mu[["gas.phi1"]], coef(f0)$mu[["gas.phi1"]],
               tolerance = 1e-5)
})

test_that("a family without the information's second derivative is rejected", {
  set.seed(2)
  dl <- data.frame(t = 1:60, y = stats::rnorm(60))
  expect_error(statmod(y ~ gas(p = 1, q = 1, time = t, scaling = 1),
                       distrib = enet_distrib(), data = dl),
               "analytic second derivative")
})

# the outer criterion of a panel whose loading is developed with a random
# effect: its gradient reads the recursion to the third order and its Hessian
# to the fourth, so the third and fourth derivatives of the information enter
sc_panel <- function(seed = 51L, m = 6L, ni = 25L) {
  set.seed(seed)
  g <- factor(rep(seq_len(m), each = ni))
  n <- m * ni
  dev <- stats::rnorm(m, 0, 0.5)
  y <- integer(n)
  for (l in seq_len(m)) {
    rows <- which(as.integer(g) == l)
    f <- 0
    a <- exp(log(0.2) + dev[l])
    for (t in seq_along(rows)) {
      lam <- exp(1 + f)
      y[rows[t]] <- stats::rpois(1, lam)
      f <- a * (y[rows[t]] - lam) / lam + 0.5 * f
    }
  }
  data.frame(y = y, g = g, t = rep(seq_len(ni), m))
}

sc_parts <- function(fit) {
  spec <- fit@spec
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  method <- fit@methods$outer
  list(spec = spec, design = design, blocks = blocks,
       idx = outer_hyper_index(spec, blocks), method = method,
       basis = integrated_basis(spec, design, method@kind),
       coef = fit@coefficients, hyper = fit@hyper)
}

# the criterion and its exact gradient at a hyperparameter, the mode refitted
# from the same start with the structural state restored before every probe
sc_outer <- function(p) {
  inner <- iwls()
  cfg <- inner_settings(inner, p$spec@distrib)
  beta0 <- unlist(p$coef[p$spec@distrib@params], use.names = FALSE)
  sst <- statmod_structural_state(p$design)
  z0 <- sst$zeta
  refit <- function(eta) {
    sst$zeta <- z0
    hy <- eta_to_hyper(eta, p$idx, p$hyper)
    r <- statmod_alternate(p$spec, p$design, p$blocks, hy, inner, beta0,
                           cfg$expected, cfg$approx, cfg$maxit, cfg$tol,
                           verbosity(0), hold_refresh = TRUE)
    list(cf = r$obj$split(r$par), hy = hy)
  }
  list(fn = function(eta) {
         r <- refit(eta)
         statmod_marginal(p$spec, p$design, r$cf, r$hy, p$method,
                          basis = p$basis)$value
       },
       gr = function(eta) {
         r <- refit(eta)
         statmod_marginal_grad(p$spec, p$design, r$cf, r$hy, p$method,
                               p$idx, p$basis)
       })
}

test_that("the outer gradient and Hessian are exact for a scaled score", {
  skip_on_cran()
  d <- sc_panel()
  fit <- statmod(y ~ gas(p = 1, q = 1, by = g, time = t, scaling = 1,
                         alpha1 ~ 1 + random(~ 1 | g)),
                 poisson_distrib(), d, outer_criterion = reml())
  p <- sc_parts(fit)
  o <- sc_outer(p)
  eta0 <- hyper_to_eta(p$hyper, p$idx)
  nh <- length(eta0)
  cdiff <- function(fn, eta, h) vapply(seq_len(nh), function(m) {
    ep <- eta; em <- eta
    ep[[m]] <- ep[[m]] + h
    em[[m]] <- em[[m]] - h
    (fn(ep) - fn(em)) / (2 * h)
  }, numeric(1))
  # the Hessian is read FIRST, at the fitted state: every refit below leaves
  # the design's structural state where that refit ended
  A <- statmod_structural_hess(p$spec, p$design, p$coef, p$hyper, p$method,
                               p$idx, p$basis)
  expect_false(is.null(A))

  # the gradient, away from the optimum where it is not nearly zero
  eta1 <- eta0 - 0.5
  expect_equal(o$gr(eta1), cdiff(o$fn, eta1, 1e-3), tolerance = 1e-4)

  # the Hessian against central differences of the exact gradient, read by
  # the rate at which the gap closes
  fd <- function(h) {
    Hn <- matrix(0, nh, nh)
    for (m in seq_len(nh)) {
      ep <- eta0; ep[[m]] <- ep[[m]] + h
      em <- eta0; em[[m]] <- em[[m]] - h
      Hn[, m] <- (o$gr(ep) - o$gr(em)) / (2 * h)
    }
    (Hn + t(Hn)) / 2
  }
  sc <- max(abs(A))
  gaps <- vapply(c(1e-2, 3e-3, 1e-3), function(h) max(abs(A - fd(h))) / sc,
                 numeric(1))
  expect_gt(gaps[[1L]] / gaps[[2L]], 4)
  expect_lt(gaps[[3L]], 1e-4)
})

test_that("a forecast of a scaled filter and its standard error are exact", {
  skip_if_not_installed("numDeriv")
  set.seed(66)
  m <- 200
  y <- numeric(m)
  f <- 1
  for (i in seq_len(m)) {
    y[i] <- f + rnorm(1, 0, 1.5)
    f <- 0.4 + 0.15 * (y[i] - f) / 1.5 + 0.6 * f
  }
  d <- data.frame(y = y, t = seq_len(m), x = rnorm(m))
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t, scaling = 0.5) - 1,
                 gaussian1_distrib(), d)
  tn <- names(fit@spec@terms$mu)[[which(vapply(fit@spec@terms$mu,
    function(z) S7::S7_inherits(z, modelterms7::structural_term),
    logical(1)))]]
  nd <- data.frame(y = NA_real_, t = m + 1:5, x = c(1, -1, 0.5, 0, 2))
  expect_warning(pr <- predict(fit, "mu", nd, se = TRUE), "parameters alone")

  # THE REFERENCE runs the filter over the extended series with the driving
  # quantity written out, u = (y - eta)/sigma^2 * sigma = (y - eta)/sigma, zero
  # past the data, and differentiates the continued rows numerically in the
  # mean's coefficient, the scale's intercept and the term's parameters
  ext <- rbind(d, nd)
  tmx <- modelterms7::term_build(fit@spec@terms$mu[[tn]], ext)
  ost <- statmod_structural_state(statmod_design(fit@spec))
  zt <- ost$zeta[[tn]]
  s2 <- statmod_respec(fit@spec, ext, need_response = FALSE)
  d2 <- statmod_design_at(s2, fit@coefficients, statmod_design(s2))
  yv <- ext$y
  fc_ref <- function(u) {
    z <- zt
    z[] <- u[-(1:2)]
    ps <- structural_psi(tmx, z)
    sgu <- exp(u[[2L]])
    e <- as.numeric(d2$mu$X %*% u[[1L]])
    as.numeric(modelterms7::term_filter(
      tmx, e, yv,
      function(e, i) if (is.na(yv[[i]])) 0 else (yv[[i]] - e) / sgu,
      function(e, i) if (is.na(yv[[i]])) 0 else -1 / sgu, ps)$eta[m + 1:5])
  }
  u0 <- c(coef(fit)$mu[["x"]], coef(fit)$sigma[["(Intercept)"]], zt)
  expect_equal(pr$fit, fc_ref(u0), tolerance = 1e-10)
  g <- numDeriv::jacobian(fc_ref, u0)
  keys <- c("mu:x", "sigma:(Intercept)", paste0("mu:gas.", names(zt)))
  V <- vcov(fit, readable = FALSE)[keys, keys]
  expect_equal(pr$se, sqrt(rowSums((g %*% V) * g)), tolerance = 1e-6)
  # the scale's column is not zero: the scaled score reads sigma
  expect_gt(max(abs(g[, 2L])), 1e-3)
})

# the term's own parameters developed by covariates, by a smooth and by a
# random effect: the derivatives of the developed coordinates go through the
# subformula route of modelterms7, which reads the same driving arrays
sc_sub_panel <- function(seed = 7, m = 6L, ni = 30L) {
  set.seed(seed)
  id <- factor(rep(seq_len(m), each = ni))
  n <- m * ni
  z <- stats::runif(n, -1, 1)
  x <- stats::rnorm(n)
  y <- integer(n)
  for (g in seq_len(m)) {
    rows <- which(as.integer(id) == g)
    f <- 0
    for (k in seq_along(rows)) {
      i <- rows[k]
      lam <- exp(0.8 + 0.2 * x[i] + f)
      y[i] <- stats::rpois(1, lam)
      a <- exp(log(0.15) + 0.4 * z[i])
      f <- a * (y[i] - lam) / sqrt(lam) + 0.6 * f
    }
  }
  data.frame(id = id, t = rep(seq_len(ni), m), x = x, z = z, y = y)
}

# at the start, every free coordinate jittered away from any optimum
sc_sub_derivs <- function(form, distrib, data, jitter = 0.07) {
  spec <- statmod_spec(form, distrib, data)
  des <- statmod_design(spec)
  sst <- statmod_structural_state(des)
  key <- names(sst$zeta)[[1L]]
  hyper <- statmod_hyper_start(spec, des)
  obj <- statmod_objective(spec, hyper, des, FALSE, "bartlett")
  beta <- statmod_start(spec, des, obj, NULL)
  nmz <- names(sst$zeta[[key]])
  free <- setdiff(nmz, sst$held[[key]])
  set.seed(1)
  z0 <- sst$zeta[[key]]
  z0[free] <- z0[free] + stats::runif(length(free), -jitter, jitter)
  b0 <- beta + stats::runif(length(beta), -jitter, jitter)
  nb <- length(b0)
  setz <- function(zz) {
    v <- z0
    v[free] <- zz
    sst$zeta[[key]] <- v
    sst$key <- NULL
  }
  u0 <- c(b0, as.numeric(z0[free]))
  L <- function(u) {
    setz(u[nb + seq_along(free)])
    statmod_loglik_at(spec, obj$split(u[seq_len(nb)]), des)
  }
  G <- function(u) {
    setz(u[nb + seq_along(free)])
    cf <- obj$split(u[seq_len(nb)])
    gz <- statmod_structural_score(spec, cf, des)[[1L]]
    c(unlist(statmod_score_at(spec, cf, des), use.names = FALSE),
      if (length(gz) == length(nmz)) gz[match(free, nmz)] else gz)
  }
  g_an <- G(u0)
  g_num <- numDeriv::grad(L, u0)
  setz(u0[nb + seq_along(free)])
  I_an <- statmod_full_information(spec, obj$split(b0), des)
  I_num <- -numDeriv::jacobian(G, u0)
  c(gradient = max(abs(g_an - g_num)) / max(1, max(abs(g_num))),
    information = max(abs(I_an - I_num)) / max(abs(I_num)))
}

test_that("a scaled filter is exact with its parameters developed", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  dp <- sc_sub_panel()
  forms <- list(
    y ~ 0 + x + gas(p = 1, q = 1, by = id, time = t, scaling = 1,
                    omega ~ z, alpha1 ~ z, pacf1 ~ z),
    y ~ x + gas(p = 1, q = 1, by = id, time = t, scaling = 0.5,
                alpha1 ~ s(z, bspline_smooth(k = 6))),
    y ~ x + gas(p = 1, q = 1, by = id, time = t, scaling = 1,
                pacf1 ~ random(~ 1 | id)))
  for (fm in forms) {
    g <- sc_sub_derivs(fm, poisson_distrib(), dp)
    expect_lt(max(g), 1e-7, label = deparse(fm[[3L]])[1L])
  }
})
