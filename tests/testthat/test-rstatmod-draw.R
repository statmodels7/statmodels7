# Everything a simulation carries is drawn, and what `par` names is held.
# The checks are against the prior the model declares -- which the draw is
# supposed to come from -- and against the fit that reads the data back.

test_that("a random effect is drawn from the prior random() declares", {
  set.seed(101)
  dd <- data.frame(id = factor(rep(1:150, each = 4)))
  s <- rstatmod(y ~ 1 + random(~1 | id), gaussian1_distrib(), dd)

  # the scale is drawn and reported, and the effects come from it
  expect_true(is.data.frame(s$hyper))
  expect_identical(nrow(s$hyper), 1L)
  expect_identical(s$hyper$name, "sigma")
  eff <- s$par$mu[-1L]
  expect_length(eff, 150L)
  expect_equal(stats::sd(eff), s$hyper$value[[1L]], tolerance = 0.15)

  # the negative control: it is NOT the plain draw of width sd, which is
  # what every coefficient got before the prior was read
  expect_false(isTRUE(all.equal(stats::sd(eff), 1, tolerance = 0.1)))
})

test_that("a Laplace prior gives effects of its own shape", {
  set.seed(102)
  dd <- data.frame(id = factor(rep(1:400, each = 2)))
  s <- rstatmod(
    y ~ 1 + random(~1 | id,
                   distrib = distributions7::fixed(
                     distributions7::laplace_distrib(), mu = 0)),
    gaussian1_distrib(), dd)
  eff <- s$par$mu[-1L]
  # a Laplace of the same variance is leptokurtic, so the draw is not the
  # Gaussian one and the excess kurtosis says so
  k <- mean((eff - mean(eff))^4) / stats::sd(eff)^4 - 3
  expect_gt(k, 1)
})

test_that("a structural term's own parameters are drawn, developments too", {
  set.seed(103)
  m <- 25
  dd <- data.frame(id = factor(rep(seq_len(m), each = 20)),
                   t = rep(seq_len(20), m))
  s <- rstatmod(y ~ 0 + gas(p = 1, q = 1, omega ~ 1 + random(~1 | id),
                            by = id, time = t),
                gaussian1_distrib(), dd)
  psi <- unlist(s$structural)
  dev <- psi[startsWith(names(psi), "omega.random.")]
  expect_length(dev, m)

  # THE DEFECT THIS CLOSES: every deviation used to be exactly zero, the
  # term's own starting value, so the panel had no heterogeneity at all
  expect_true(all(dev != 0))
  expect_gt(stats::sd(dev), 0.05)
  # and the level really does vary by group, which is what the formula says
  expect_gt(stats::sd(tapply(s$latent, dd$id, mean)), 0.05)
  # far more than the three distinct values the held-at-zero version gave
  expect_gt(length(unique(round(s$theta$mu, 8))), 100L)
})

test_that("par holds an equation, a coefficient, a group or a parameter", {
  set.seed(104)
  dd <- data.frame(x = stats::runif(60),
                   id = factor(rep(1:10, each = 6)))

  # a whole equation, as it always did
  a <- rstatmod(y ~ x, gaussian1_distrib(), dd,
                par = list(mu = c(1, 2), sigma = log(0.3)))
  expect_equal(unname(a$par$mu), c(1, 2))

  # one coefficient of it, the rest drawn
  b <- rstatmod(y ~ x, gaussian1_distrib(), dd,
                par = list("mu.(Intercept)" = 5))
  expect_equal(unname(b$par$mu[[1L]]), 5)
  expect_false(b$par$mu[[2L]] == 5)

  # a group, by the name its members extend at a dot
  cc <- rstatmod(y ~ 1 + random(~1 | id), gaussian1_distrib(), dd,
                 par = list(mu.random = function(k) rep(0.7, k)))
  expect_equal(unname(cc$par$mu[-1L]), rep(0.7, 10))

  # a structural parameter, on the scale a reader knows
  d <- rstatmod(y ~ 0 + gas(p = 1, q = 1, time = t), gaussian1_distrib(),
                data.frame(t = 1:60), par = list(alpha1 = 0.35, sigma = 0))
  expect_equal(d$structural$alpha1, 0.35)
  # and its neighbours are still drawn
  expect_false(d$structural$omega == 0)
})

test_that("a structural group is held by one name", {
  set.seed(105)
  m <- 8
  dd <- data.frame(id = factor(rep(seq_len(m), each = 10)),
                   t = rep(1:10, m))
  s <- rstatmod(y ~ 0 + gas(p = 1, q = 1, omega ~ 1 + random(~1 | id),
                            by = id, time = t),
                gaussian1_distrib(), dd,
                par = list(omega.random = function(k) seq_len(k) / 10))
  psi <- unlist(s$structural)
  expect_equal(unname(psi[startsWith(names(psi), "omega.random.")]),
               seq_len(m) / 10)
})

test_that("a hyperparameter the term holds is used rather than drawn", {
  set.seed(106)
  dd <- data.frame(x = stats::runif(80))
  s <- rstatmod(y ~ s(x, k = 6, lambda = 3), gaussian1_distrib(), dd)
  expect_equal(s$hyper$value[s$hyper$name == "lambda"], 3)
  expect_true(all(s$hyper$held))
})

test_that("the truth is drawn once and shared across replicates", {
  set.seed(107)
  dd <- data.frame(id = factor(rep(1:20, each = 5)))
  s <- rstatmod(y ~ 1 + random(~1 | id), gaussian1_distrib(), dd,
                n_sim = 3)
  expect_length(s$data, 3L)
  expect_identical(nrow(s$hyper), 1L)
  # one set of coefficients, three data sets
  expect_length(s$par$mu, 21L)
  expect_false(isTRUE(all.equal(s$data[[1L]]$y, s$data[[2L]]$y)))
})

test_that("a key that addresses nothing is refused with what there is", {
  dd <- data.frame(x = stats::runif(20))
  expect_error(rstatmod(y ~ x, gaussian1_distrib(), dd,
                        par = list(nonesuch = 1)),
               "addresses nothing")
  expect_error(rstatmod(y ~ x, gaussian1_distrib(), dd,
                        par = list(nonesuch = 1)),
               "mu, sigma")
  # a group that answers with the wrong count is refused rather than recycled
  expect_error(rstatmod(y ~ x, gaussian1_distrib(), dd,
                        par = list(mu = function(k) 1)),
               "addresses 2")
})

test_that("a simulated random effect is recovered by a fit of it", {
  set.seed(108)
  dd <- data.frame(id = factor(rep(1:60, each = 8)))
  s <- rstatmod(y ~ 1 + random(~1 | id), gaussian1_distrib(), dd,
                par = list(mu.random = function(k) stats::rnorm(k, 0, 0.8),
                           "mu.(Intercept)" = 1, sigma = log(0.5)))
  fit <- statmod(y ~ 1 + random(~1 | id), gaussian1_distrib(), s$data)
  b <- coef(fit, readable = FALSE)$mu
  expect_equal(unname(b[[1L]]), 1, tolerance = 0.25)
  expect_equal(stats::cor(unname(b[-1L]), unname(s$par$mu[-1L])), 1,
               tolerance = 0.1)
})

test_that("drawing everything leaves an ordinary model where it was", {
  # the control for the chart rule: a prior that does not ride a chart is
  # centred where it always was, so an ordinary shape is unchanged
  set.seed(109)
  dd <- data.frame(x = stats::runif(200))
  v <- replicate(30, {
    s <- rstatmod(y ~ s(x, k = 8), gaussian1_distrib(), dd)
    c(stats::sd(s$data$y), s$hyper$value[[1L]])
  })
  expect_gt(stats::median(v[2L, ]), 0.5)
  expect_lt(stats::median(v[2L, ]), 2.5)
  expect_true(all(is.finite(v)))
})

test_that("a covariance class is drawn group by group, not shuffled", {
  # The class's prior is over the PAIR one group carries across the two
  # equations, and its penalty reads those pairs interleaved. A scatter that
  # got the interleaving wrong would put one group's mean effect beside
  # another group's scale effect, and the sampled correlation would collapse
  # to zero whatever the prior says.
  m <- 1500
  dd <- data.frame(g = factor(rep(seq_len(m), each = 2)))
  f <- y ~ 1 + random(~1 | u | g) | sigma ~ 1 + random(~1 | u | g)
  sp <- statmod_spec(f, gaussian1_distrib(),
                     transform(dd, y = stats::rnorm(nrow(dd))))
  u <- statmod_penalized(sp, statmod_design(sp))[[1L]]
  # the branch under test: a class is addressed by stacked positions and has
  # no `cols` of its own
  expect_null(u$cols)
  expect_length(u$index, 2L * m)

  se <- 1 / sqrt(m)
  for (seed in 1:3) {
    set.seed(seed)
    s <- rstatmod(f, gaussian1_distrib(), dd)
    a <- unname(s$par$mu[-1L])
    b <- unname(s$par$sigma[-1L])
    th <- stats::setNames(as.list(s$hyper$value), s$hyper$name)
    # the covariance those same hyperparameters describe, built by the prior
    # rather than by the draw
    S <- as.matrix(distributions7::mv_sigma(u$penalty@parent, th))
    rho <- S[1L, 2L] / sqrt(S[1L, 1L] * S[2L, 2L])
    expect_equal(stats::sd(a), sqrt(S[1L, 1L]), tolerance = 0.06)
    expect_equal(stats::sd(b), sqrt(S[2L, 2L]), tolerance = 0.06)
    expect_lt(abs(stats::cor(a, b) - rho), 4 * se)
    # the negative control, where the prior says the two are related at all:
    # breaking the pairing must fail the same check
    if (abs(rho) > 0.15) {
      expect_gt(abs(stats::cor(a, sample(b)) - rho), 4 * se)
    }
  }
})
