skip_on_cran()

orange_fit <- function() {
  orange <- as.data.frame(datasets::Orange)
  orange$Tree <- factor(as.character(orange$Tree))
  list(data = orange,
       fit = statmod(circumference ~ 0 +
                       nl(~ Asym / (1 + exp((xmid - age) / scal)),
                          Asym ~ 1 + random(~ 1 | Tree)),
                     distrib = distributions7::gaussian1_distrib(),
                     data = orange))
}

test_that("random = 'zero' reaches a random effect inside nl()", {
  o <- orange_fit()
  f <- o$fit
  cf <- coef(f)$mu
  pop <- function(p, age) p[1] / (1 + exp((p[2] - age) / p[3]))
  p0 <- c(cf[["nl.Asym.(Intercept)"]], cf[["nl.xmid"]], cf[["nl.scal"]])
  # a tree the fit saw and one it did not, both read as the typical tree
  nd <- data.frame(age = c(500, 1000, 1500), Tree = factor(c("1", "1", "9")))
  expect_equal(predict(f, what = "mu", newdata = nd, random = "zero"),
               pop(p0, nd$age), tolerance = 1e-12)
  # the conditional prediction of tree 1 is not the typical tree's
  cond <- predict(f, what = "mu", newdata = nd[1:2, ])
  expect_gt(max(abs(cond - pop(p0, nd$age[1:2]))), 1)
  # at the fitting rows
  expect_equal(predict(f, what = "mu", random = "zero"),
               pop(p0, o$data$age), tolerance = 1e-12)
  # the delta method on the three fixed coefficients, by numDeriv
  V <- vcov(f)
  k <- c("mu:nl.Asym.(Intercept)", "mu:nl.xmid", "mu:nl.scal")
  se_h <- vapply(nd$age, function(a) {
    g <- numDeriv::grad(function(p) pop(p, a), p0)
    sqrt(drop(t(g) %*% V[k, k] %*% g))
  }, numeric(1))
  zs <- predict(f, what = "mu", newdata = nd, random = "zero", se = TRUE)
  expect_equal(zs$se, se_h, tolerance = 1e-6)
  # the curve is linear in Asym, so the average over the prior is the
  # typical tree's curve
  expect_equal(predict(f, what = "mu", newdata = nd, random = "marginal"),
               pop(p0, nd$age), tolerance = 1e-10)
})

test_that("random = 'marginal' averages a nonlinear term over its prior", {
  th <- as.data.frame(datasets::Theoph)
  th <- th[th$Time > 0, ]
  th$Subject <- factor(as.character(th$Subject))
  f <- statmod(conc ~ 0 + nl(~ Dose * exp(lKe + lKa - lCl) *
                               (exp(-exp(lKe) * Time) - exp(-exp(lKa) * Time)) /
                               (exp(lKa) - exp(lKe)),
                             lKa ~ 1 + random(~ 1 | Subject),
                             lCl ~ 1 + random(~ 1 | Subject),
                             start = list(lKe = -2.5, lKa = 0.5, lCl = -3)),
               distrib = distributions7::gaussian1_distrib(), data = th)
  cf <- coef(f)$mu
  tau <- hyper(f)$estimate
  fol <- function(t, p, ua, uc) {
    4.5 * exp(p[1] + p[2] + ua - p[3] - uc) *
      (exp(-exp(p[1]) * t) - exp(-exp(p[2] + ua) * t)) /
      (exp(p[2] + ua) - exp(p[1]))
  }
  p0 <- c(cf[["nl.lKe"]], cf[["nl.lKa.(Intercept)"]], cf[["nl.lCl.(Intercept)"]])
  nd <- data.frame(Time = c(0.5, 2, 8, 20), Dose = 4.5, Subject = factor("99"))
  # the typical subject
  expect_equal(predict(f, what = "mu", newdata = nd, random = "zero"),
               fol(nd$Time, p0, 0, 0), tolerance = 1e-12)
  # the population average, by a Gauss-Hermite product grid written here
  k <- 40L
  e <- eigen(local({
    J <- matrix(0, k, k)
    i <- seq_len(k - 1L)
    J[cbind(i, i + 1L)] <- J[cbind(i + 1L, i)] <- sqrt(i / 2)
    J
  }), symmetric = TRUE)
  x <- e$values * sqrt(2)
  w <- e$vectors[1L, ]^2
  avg <- function(p) {
    vapply(nd$Time, function(t) {
      sum(outer(w, w) * outer(x * tau[1], x * tau[2],
                              function(a, b) fol(t, p, a, b)))
    }, numeric(1))
  }
  m <- predict(f, what = "mu", newdata = nd, random = "marginal", se = TRUE)
  expect_equal(m$fit, avg(p0), tolerance = 1e-10)
  # its delta method, the hyperparameters held
  V <- vcov(f)
  kk <- c("mu:nl.lKe", "mu:nl.lKa.(Intercept)", "mu:nl.lCl.(Intercept)")
  G <- numDeriv::jacobian(avg, p0)
  expect_equal(m$se, sqrt(rowSums((G %*% V[kk, kk]) * G)), tolerance = 1e-5)
  # one of the two set aside by name: the other stays conditional, so a
  # subject the fit saw is needed
  key <- statmodels7:::nested_random_terms(f@spec)$key
  expect_length(key, 2L)
  nd1 <- data.frame(Time = 2, Dose = 4.5, Subject = factor("1"))
  ua1 <- cf[["nl.lKa.random.1"]]
  uc1 <- cf[["nl.lCl.random.1"]]
  expect_equal(predict(f, what = "mu", newdata = nd1,
                       random = stats::setNames("zero", key[1L])),
               fol(2, p0, 0, uc1), tolerance = 1e-12)
  expect_equal(predict(f, what = "mu", newdata = nd1,
                       random = stats::setNames("zero", key[2L])),
               fol(2, p0, ua1, 0), tolerance = 1e-12)
  expect_error(predict(f, what = "mu", newdata = nd, random = "zero",
                       interval = "group"), "subformula")
})

test_that("random = 'zero' reaches a random break-point", {
  set.seed(3)
  id <- factor(rep(1:12, each = 10))
  t <- rep(seq(0, 10, length.out = 10), 12)
  psi_i <- 5 + rnorm(12, sd = 0.8)
  y <- 2 + 0.8 * t - 1.2 * pmax(t - psi_i[id], 0) + rnorm(120, sd = 0.3)
  f <- statmod(y ~ seg(t, psi ~ random(~ 1 | id)),
               distrib = distributions7::gaussian1_distrib(),
               data = data.frame(id = id, t = t, y = y))
  cf <- coef(f)$mu
  nd <- data.frame(t = c(2, 4.5, 6, 9), id = factor("50"))
  ref <- cf[["(Intercept)"]] + cf[["seg.beta"]] * nd$t +
    cf[["seg.gamma1"]] * pmax(nd$t - cf[["seg.psi1.(Intercept)"]], 0)
  expect_equal(predict(f, what = "mu", newdata = nd, random = "zero"), ref,
               tolerance = 1e-10)
})
