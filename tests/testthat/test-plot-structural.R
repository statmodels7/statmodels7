seg_fit <- function() {
  set.seed(1)
  dd <- data.frame(x = sort(runif(120, 0, 10)))
  dd$y <- 1 + 0.5 * dd$x + 2 * pmax(dd$x - 6, 0) + rnorm(120, sd = 0.4)
  list(fit = statmod(y ~ seg(x, npsi = 1), distributions7::gaussian1_distrib(),
                     dd), data = dd, name = "seg(x, npsi = 1)")
}

gas_fit <- function() {
  set.seed(4)
  n <- 90
  x <- rnorm(n)
  f <- 0
  y <- integer(n)
  for (t in seq_len(n)) {
    lam <- exp(0.7 + 0.3 * x[t] + f)
    y[t] <- rpois(1, lam)
    f <- 0.15 * (y[t] - lam) * lam^(-0.5) + 0.8 * f
  }
  dd <- data.frame(t = seq_len(n), x = x, y = y)
  nm <- "gas(p = 1, q = 1, time = t, scaling = 0.5)"
  list(fit = statmod(y ~ x + gas(p = 1, q = 1, time = t, scaling = 0.5),
                     distributions7::poisson_distrib(), dd),
       data = dd, name = nm)
}

reg_fit <- function() {
  set.seed(4)
  n <- 80
  s <- integer(n)
  s[1] <- 1L
  for (t in 2:n) s[t] <- if (runif(1) < 0.9) s[t - 1] else 3L - s[t - 1]
  dd <- data.frame(t = seq_len(n), x = rnorm(n))
  dd$y <- c(0, 3)[s] + 0.5 * dd$x + rnorm(n, sd = 0.7)
  list(fit = statmod(y ~ x + regime(k = 2, time = t),
                     distributions7::gaussian1_distrib(), dd),
       data = dd, name = "regime(k = 2, time = t)")
}

sg <- seg_fit()

test_that("the effect of a break-point term is the broken line", {
  tm <- sg$fit@spec@terms$mu[[sg$name]]
  d <- plot_data(plot_effect(), tm, sg$fit,
                 ctx = list(param = "mu", name = sg$name, vars = "x",
                            cols = NULL, data = sg$data))
  b <- coef(sg$fit, readable = FALSE)$mu
  beta <- b[["seg.beta"]]
  gam <- b[["seg.gamma1"]]
  psi <- b[["seg.psi1"]]
  manual <- beta * d$curve$x + gam * pmax(d$curve$x - psi, 0)
  expect_equal(d$curve$fit, manual, tolerance = 1e-10)
  expect_equal(d$curve$fit, modelterms7::term_value(tm, newdata =
    data.frame(x = d$curve$x)), tolerance = 1e-12)
  # the band is the delta method with the derivative of the broken line in
  # (beta, gamma, psi)
  V <- vcov(sg$fit, readable = FALSE)[paste0("mu:seg.", c("beta", "gamma1",
                                                         "psi1")),
                                       paste0("mu:seg.", c("beta", "gamma1",
                                                           "psi1"))]
  xg <- d$curve$x
  J <- cbind(xg, pmax(xg - psi, 0), -gam * (xg > psi))
  se <- sqrt(rowSums((J %*% V) * J))
  expect_equal(d$curve$upper - d$curve$fit, stats::qnorm(0.975) * se,
               tolerance = 1e-8)
  # partial residuals: the contribution at the data plus the residual
  ofit <- beta * sg$data$x + gam * pmax(sg$data$x - psi, 0)
  expect_equal(d$points$y, ofit + sg$data$y - fitted(sg$fit, "mu"),
               tolerance = 1e-8)
})

test_that("the effect of a nonlinear term is its function", {
  set.seed(4)
  dn <- data.frame(x = runif(100, 0, 3))
  dn$y <- 5 * exp(-0.8 * dn$x) + rnorm(100, sd = 0.2)
  fn <- statmod(y ~ 0 + nl(~ a * exp(-r * x), start = list(a = 4, r = 1),
                           links = list(r = linkfunctions7::log_link())),
                distributions7::gaussian1_distrib(), dn)
  nm <- names(fn@spec@terms$mu)[1L]
  p <- autoplot(fn, term = nm, data = dn)
  d <- p$layers[[2]]$data
  b <- coef(fn, readable = FALSE)$mu
  expect_equal(d$fit, b[["nl.a"]] * exp(-exp(b[["nl.r"]]) * d$x),
               tolerance = 1e-10)
  expect_true(all(is.finite(d$lower) & d$lower < d$fit))
})

test_that("the break-point graph carries the position, the interval and the profile", {
  d <- plot_data(plot_breakpoint(), sg$fit@spec@terms$mu[[sg$name]], sg$fit,
                 ctx = list(param = "mu", name = sg$name, vars = "x",
                            cols = NULL, data = sg$data))
  psi <- coef(sg$fit, readable = FALSE)$mu[["seg.psi1"]]
  expect_equal(d$breaks$psi, psi, tolerance = 1e-10)
  ci <- confint(sg$fit, parm = "mu:seg.psi1", test = lr_test(),
                readable = FALSE)
  expect_equal(c(d$breaks$lower, d$breaks$upper),
               c(ci$lower, ci$upper), tolerance = 1e-8)
  # the profile is zero at its minimum, which is the fitted position up to
  # the grid, and crosses the cutoff on both sides of the interval
  expect_equal(min(d$profile$delta), 0)
  expect_equal(d$profile$x[which.min(d$profile$delta)], psi, tolerance = 0.2)
  expect_equal(d$cutoff$cutoff, stats::qchisq(0.95, 1))
  inside <- d$profile$x > d$breaks$lower & d$profile$x < d$breaks$upper
  expect_true(mean(d$profile$delta[inside]) < d$cutoff$cutoff)
  expect_true(all(d$profile$delta[d$profile$x < psi - 1.5 |
                                    d$profile$x > psi + 1.5] >
                    d$cutoff$cutoff))
  pdf(NULL)
  on.exit(dev.off())
  p <- autoplot(sg$fit, term = sg$name, data = sg$data)
  expect_s3_class(p, "ggplot")
  expect_no_error(print(p))
  q <- autoplot(sg$fit, term = sg$name, data = sg$data,
                plot_type = plot_breakpoint(profile = FALSE, test = wald_test()))
  wi <- confint(sg$fit, parm = "mu:seg.psi1", readable = FALSE)
  d2 <- q$layers[[5]]$data
  expect_equal(c(d2$lower, d2$upper), c(wi$lower, wi$upper), tolerance = 1e-8)
  expect_error(plot_breakpoint(profile = NA), "TRUE or FALSE")
  expect_error(autoplot(sg$fit, term = sg$name, plot_type = plot_breakpoint()),
               "does not keep its data")
})

gs <- gas_fit()

test_that("the filter graph draws the fitted path and its band", {
  p <- autoplot(gs$fit, term = gs$name)
  d <- plot_data(plot_filter(), gs$fit@spec@terms$mu[[gs$name]], gs$fit,
                 ctx = list(param = "mu", name = gs$name, vars = NULL,
                            cols = NULL, data = NULL))
  expect_equal(d$path$fit, gs$fit@fitted$mu, tolerance = 1e-12)
  pr <- predict(gs$fit, what = "mu", se = TRUE)
  expect_equal(d$path$lower, pr$lower, tolerance = 1e-12)
  expect_true(all(d$path$lower < d$path$fit & d$path$fit < d$path$upper))
  expect_equal(d$path$x, gs$data$t)
  dl <- plot_data(plot_filter(scale = "link", level = 0.9),
                  gs$fit@spec@terms$mu[[gs$name]], gs$fit,
                  ctx = list(param = "mu", name = gs$name))
  expect_equal(dl$path$fit, log(gs$fit@fitted$mu), tolerance = 1e-12)
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(print(autoplot(gs$fit, term = gs$name,
                                 plot_type = plot_filter(response = TRUE))))
  expect_error(plot_filter(scale = "link", response = TRUE), "scale")
})

test_that("a panel is drawn by group and too many groups are named", {
  set.seed(2)
  m <- 4
  ni <- 40
  dp <- data.frame(id = factor(rep(1:m, each = ni)), t = rep(seq_len(ni), m))
  dp$y <- rpois(m * ni, 4)
  fp <- statmod(y ~ gas(p = 1, q = 1, by = id, time = t),
                distributions7::poisson_distrib(), dp)
  nm <- names(fp@spec@terms$mu)[2L]
  d <- plot_data(plot_filter(), fp@spec@terms$mu[[nm]], fp,
                 ctx = list(param = "mu", name = nm))
  expect_equal(levels(d$path$group), as.character(1:m))
  expect_equal(d$path$fit, fp@fitted$mu, tolerance = 1e-12)
  d2 <- plot_data(plot_filter(groups = "2"), fp@spec@terms$mu[[nm]], fp,
                  ctx = list(param = "mu", name = nm))
  expect_equal(nrow(d2$path), ni)
  expect_error(plot_data(plot_filter(groups = "9"), fp@spec@terms$mu[[nm]],
                         fp, ctx = list(param = "mu", name = nm)),
               "does not have")
})

test_that("the news impact is the loading times the scaled score", {
  tm <- gs$fit@spec@terms$mu[[gs$name]]
  d <- plot_data(plot_news_impact(), tm, gs$fit,
                 ctx = list(param = "mu", name = gs$name))
  kappa <- gs$fit@structural[[gs$name]]$parameter[["kappa1"]]
  mu <- gs$fit@fitted$mu
  ref <- d$reference$row
  expect_equal(mu[ref], stats::median(mu), tolerance = 0.05)
  zero <- d$curve[d$curve$z == 0, ]
  # a Poisson with the log link: score y - mu, information mu
  y0 <- stats::qpois(0.5, mu[ref])
  expect_equal(zero$y, y0)
  expect_equal(zero$impact, kappa * (y0 - mu[ref]) * mu[ref]^(-0.5),
               tolerance = 1e-10)
  # a count family gives a step function of z
  expect_true(attr(d, "discrete"))
  expect_true(all(diff(d$curve$impact) >= 0))
  d3 <- plot_data(plot_news_impact(at = 3), tm, gs$fit,
                  ctx = list(param = "mu", name = gs$name))
  expect_equal(d3$reference$row, 3)
  expect_error(plot_news_impact(range = c(1, 0)), "increasing")
  expect_error(plot_data(plot_news_impact(at = 1000), tm, gs$fit,
                         ctx = list(param = "mu", name = gs$name)),
               "1000")
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(print(autoplot(gs$fit, term = gs$name,
                                 plot_type = plot_news_impact())))
})

rg <- reg_fit()

test_that("the state probabilities are distributions that agree with statmod_latent()", {
  tm <- rg$fit@spec@terms$mu[[rg$name]]
  ctx <- list(param = "mu", name = rg$name, data = rg$data)
  ds <- plot_data(plot_states(), tm, rg$fit, ctx = ctx)
  df <- plot_data(plot_states(smoothed = FALSE), tm, rg$fit, ctx = ctx)
  n <- nrow(rg$data)
  sm <- matrix(ds$states$prob, n, 2)
  fl <- matrix(df$states$prob, n, 2)
  expect_equal(rowSums(sm), rep(1, n), tolerance = 1e-12)
  expect_equal(rowSums(fl), rep(1, n), tolerance = 1e-12)
  expect_equal(sm, as.matrix(statmod_latent(rg$fit)), ignore_attr = TRUE,
               tolerance = 1e-12)
  # the two kinds of probability agree at the last observation
  expect_equal(fl[n, ], sm[n, ], tolerance = 1e-10)
  expect_false(isTRUE(all.equal(fl, sm, tolerance = 1e-6)))
  # the forward recursion reproduces the log-likelihood of the term
  r <- regime_entry(rg$fit, rg$name)
  lk <- modelterms7::term_loglik(r$tm, r$eta_static, rg$fit@spec@response,
                                 r$cb$logdens, r$cb$score, r$psi)$loglik
  expect_equal(regime_filtered(r)$loglik, lk, tolerance = 1e-10)
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(print(autoplot(rg$fit, term = rg$name, data = rg$data)))
})

test_that("the regime densities are the family's density at each level", {
  tm <- rg$fit@spec@terms$mu[[rg$name]]
  d <- plot_data(plot_regime_density(at = 5), tm, rg$fit,
                 ctx = list(param = "mu", name = rg$name))
  r <- regime_entry(rg$fit, rg$name)
  rp <- regime_parameters(r)
  sg1 <- rg$fit@fitted$sigma[5]
  for (j in 1:2) {
    dj <- d$density[d$density$regime == paste("regime", j), ]
    expect_equal(dj$density, stats::dnorm(dj$y, r$eta_static[5] + rp$mu[j],
                                          sg1), tolerance = 1e-10)
  }
  expect_equal(d$reference$level, rp$mu)
  expect_equal(rp$mu[2] - rp$mu[1], unlist(r$psi)[["gap2"]])
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(print(autoplot(rg$fit, term = rg$name,
                                 plot_type = plot_regime_density())))
})

test_that("the default graph of each term is its own", {
  expect_identical(default_plot(sg$fit@spec@terms$mu[[sg$name]])@label,
                   "plot_breakpoint()")
  expect_identical(default_plot(gs$fit@spec@terms$mu[[gs$name]])@label,
                   "plot_filter()")
  expect_identical(default_plot(rg$fit@spec@terms$mu[[rg$name]])@label,
                   "plot_states()")
})

test_that("the interval of a jump() position is read from its profile", {
  set.seed(2)
  dj <- data.frame(x = sort(runif(150, 0, 10)))
  dj$y <- 1 + 0.3 * dj$x + 2 * (dj$x > 6) + rnorm(150, sd = 0.4)
  fj <- statmod(y ~ jump(x, npsi = 1), distributions7::gaussian1_distrib(), dj)
  nm <- grep("jump", names(fj@spec@terms$mu), value = TRUE)
  tm <- fj@spec@terms$mu[[nm]]
  ctx <- list(param = "mu", name = nm, vars = "x", data = dj)
  d <- plot_data(plot_breakpoint(profile = FALSE), tm, fj, ctx = ctx)
  expect_equal(nrow(d$profile), 0)
  expect_true(d$breaks$lower <= d$breaks$psi && d$breaks$psi <= d$breaks$upper)
  dp <- plot_data(plot_breakpoint(), tm, fj, ctx = ctx)
  inside <- dp$profile$x >= dp$breaks$lower & dp$profile$x <= dp$breaks$upper
  expect_true(all(dp$profile$delta[inside] <= dp$cutoff$cutoff))
  expect_equal(c(dp$breaks$lower, dp$breaks$upper),
               c(d$breaks$lower, d$breaks$upper))
})

test_that("profile_interval keeps the run connected to the estimate", {
  x <- 1:9
  delta <- c(0, 5, 1, 0.5, 0, 0.5, 1, 5, 0)
  expect_equal(profile_interval(x, delta, psi = 5, cutoff = 3.84), c(3, 7))
})
