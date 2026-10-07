test_that("a panel whose level is developed is continued past its series", {
  # The span check of a developed level compared the development, built on
  # the rows of the fit, with the design at the new rows, and stopped every
  # forecast with "'qr' and 'y' must have the same number of rows".
  set.seed(5)
  m <- 6L
  per <- 40L
  dd <- data.frame(id = factor(rep(seq_len(m), each = per)),
                   t = rep(seq_len(per), m), x = stats::rnorm(m * per))
  om <- stats::rnorm(m, sd = 0.5)
  dd$y <- 0
  for (g in seq_len(m)) {
    r <- which(dd$id == g)
    f <- om[g] / (1 - 0.6)
    for (k in seq_along(r)) {
      if (k > 1L) f <- om[g] + 0.3 * e + 0.6 * f
      e <- stats::rnorm(1)
      dd$y[r[k]] <- 0.5 * dd$x[r[k]] + f + e
    }
  }
  fit <- statmod(y ~ 0 + x + gas(p = 1, q = 1, omega ~ 1 + random(~ 1 | id),
                                by = id, time = t),
                 distributions7::gaussian1_distrib(), dd)
  nd <- data.frame(id = factor(2L, levels = seq_len(m)), t = per + 1:3,
                   x = c(0.2, -0.1, 0.4))
  p <- predict(fit, what = "mu", newdata = nd)
  # the continuation written out: the filter over group 2 at the estimates,
  # then the recursion with the score at its mean of zero
  cf <- coef(fit)$mu
  s2 <- exp(2 * coef(fit)$sigma[["(Intercept)"]])
  a <- cf[["gas.alpha1"]]
  b <- cf[["gas.phi1"]]
  w <- cf[["gas.omega.(Intercept)"]] + cf[["gas.omega.random.2"]]
  d2 <- dd[dd$id == 2L, ]
  f <- w / (1 - b)
  for (k in seq_len(per)) {
    if (k > 1L) f <- w + a * (d2$y[k - 1L] - cf[["x"]] * d2$x[k - 1L] - fp) / s2 +
      b * fp
    fp <- f
  }
  fc <- numeric(3)
  fc[1] <- w + a * (d2$y[per] - cf[["x"]] * d2$x[per] - f) / s2 + b * f
  for (h in 2:3) fc[h] <- w + b * fc[h - 1L]
  expect_equal(p, cf[["x"]] * nd$x + fc, tolerance = 1e-10)
  # and one row is the first of them
  expect_equal(predict(fit, what = "mu", newdata = nd[1, ]), p[1],
               tolerance = 1e-12)
})
