ex_info <- function() distributions7::expected_info()

fvx <- function(sw) frequentist_vcov(information = ex_info(), sandwich = sw)

# (I - H)^p of a symmetric matrix, written here from the eigen decomposition
mpow <- function(H, p) {
  e <- eigen((H + t(H)) / 2, symmetric = TRUE)
  e$vectors %*% ((1 - e$values)^p * t(e$vectors))
}

lm_fit <- function(n = 60, seed = 11) {
  set.seed(seed)
  dd <- data.frame(x = runif(n), z = rnorm(n))
  dd$y <- 1 + 2 * dd$x - 0.5 * dd$z + rnorm(n, sd = 0.2 + dd$x)
  dd
}

test_that("the leverage of a Gaussian model equals hatvalues of lm", {
  dd <- lm_fit()
  f <- statmod(y ~ x + z, distributions7::gaussian1_distrib(), dd)
  lv <- statmod_leverage(f, information = ex_info())
  h <- unname(hatvalues(lm(y ~ x + z, dd)))
  expect_equal(unname(lv$leverage[, "mu"]), h, tolerance = 1e-10)
  # the block is symmetric and its trace sums to the number of coefficients
  expect_equal(lv$blocks[, 1, 2], lv$blocks[, 2, 1], ignore_attr = TRUE)
  expect_equal(sum(lv$trace), 4, tolerance = 1e-9)
  # the sigma equation: an intercept only, leverage 1 / n for each row
  expect_equal(unname(lv$leverage[, "sigma"]), rep(1 / nrow(dd), nrow(dd)),
               tolerance = 1e-9)
})

test_that("the leverage of one equation equals that of plot_influence()", {
  dd <- lm_fit(80, 5)
  f <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
  lv <- statmod_leverage(f, information = ex_info())
  inf <- influence_deltas(f, "mu", ex_info())
  expect_equal(unname(lv$leverage[, "mu"]), inf$leverage, tolerance = 1e-10)
  # the sigma coordinates of a Gaussian under REML are read there from the
  # marginal criterion, so the two differ by construction; a fit with a
  # Poisson mean has no such coordinate
  set.seed(9)
  dp <- data.frame(x = rnorm(90))
  dp$y <- rpois(90, exp(0.2 + 0.4 * dp$x))
  fp <- statmod(y ~ x, distributions7::poisson_distrib(), dp)
  lp <- statmod_leverage(fp, information = ex_info())
  expect_equal(unname(lp$leverage[, "mu"]),
               influence_deltas(fp, "mu", ex_info())$leverage, tolerance = 1e-10)
})

test_that("the block of a location-scale model is the written-out formula", {
  dd <- lm_fit(50, 8)
  f <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
  X <- cbind(1, dd$x)
  sig <- exp(drop(X %*% f@coefficients$sigma))
  n <- nrow(dd)
  # expected information of (mu, log sigma): diag(1 / sigma^2, 2)
  K <- matrix(0, 4, 4)
  for (i in seq_len(n)) {
    Zi <- rbind(c(X[i, ], 0, 0), c(0, 0, X[i, ]))
    K <- K + t(Zi) %*% diag(c(1 / sig[i]^2, 2)) %*% Zi
  }
  Ki <- solve(K)
  lv <- statmod_leverage(f, information = ex_info())
  for (i in c(1, 17, 50)) {
    Zi <- rbind(c(X[i, ], 0, 0), c(0, 0, X[i, ]))
    Wh <- diag(sqrt(c(1 / sig[i]^2, 2)))
    Hi <- Wh %*% Zi %*% Ki %*% t(Zi) %*% Wh
    expect_equal(unname(lv$blocks[i, , ]), unname(Hi), tolerance = 1e-8)
  }
  expect_equal(sum(lv$trace), 4, tolerance = 1e-8)
})

test_that("a penalized fit's leverage sums to its effective degrees of freedom", {
  set.seed(4)
  dd <- data.frame(x = runif(120))
  dd$y <- sin(2 * pi * dd$x) + rnorm(120, sd = 0.3)
  f <- statmod(y ~ s(x, basis7::bspline_smooth(k = 8)),
               distributions7::gaussian1_distrib(), dd)
  lv <- statmod_leverage(f, information = ex_info())
  expect_gt(sum(lv$trace), 3)
  expect_lt(sum(lv$trace), 12)
  expect_true(all(lv$leverage[, "mu"] > 0 & lv$leverage[, "mu"] < 1))
  expect_error(statmod_leverage("a"), "fitted")
})

test_that("HC2, HC3 and HC4 of a Gaussian linear model are the textbook formulas", {
  dd <- lm_fit(70, 2)
  f <- statmod(y ~ x + z, distributions7::gaussian1_distrib(), dd)
  m <- lm(y ~ x + z, dd)
  X <- model.matrix(m)
  e <- unname(resid(m))
  h <- unname(hatvalues(m))
  A <- solve(crossprod(X))
  tb <- function(w) A %*% crossprod(X * sqrt(w)) %*% A
  # the exponent of HC4 reads the trace over both equations: the scale
  # equation adds 1 / n to every row
  ht <- h + 1 / length(h)
  d4 <- pmin(4, ht / mean(ht))
  hand <- list(HC2 = tb(e^2 / (1 - h)), HC3 = tb(e^2 / (1 - h)^2),
               HC4 = tb(e^2 / (1 - h)^d4))
  sws <- list(HC2 = hc2_sandwich(), HC3 = hc3_sandwich(), HC4 = hc4_sandwich())
  for (k in names(sws)) {
    v <- vcov(f, vcov_type = fvx(sws[[k]]), param = "mu", readable = FALSE)
    expect_equal(unname(v), unname(hand[[k]]), tolerance = 1e-8)
  }
  expect_true(all(diag(hand$HC3) > diag(hand$HC2)))
})

test_that("HC2 and HC3 of a Poisson model equal the glm formulas", {
  set.seed(6)
  n <- 150
  x <- rnorm(n)
  dp <- data.frame(x = x, y = rpois(n, exp(0.4 + 0.5 * x) * rgamma(n, 3, 3)))
  f <- statmod(y ~ x, distributions7::poisson_distrib(), dp)
  g <- glm(y ~ x, family = poisson, data = dp)
  X <- model.matrix(g)
  mu <- unname(fitted(g))
  h <- unname(hatvalues(g))
  e <- dp$y - mu
  A <- solve(crossprod(X * sqrt(mu)))
  tb <- function(w) A %*% crossprod(X * sqrt(w)) %*% A
  v2 <- vcov(f, vcov_type = frequentist_vcov(sandwich = hc2_sandwich()),
             readable = FALSE)
  v3 <- vcov(f, vcov_type = frequentist_vcov(sandwich = hc3_sandwich()),
             readable = FALSE)
  v4 <- vcov(f, vcov_type = frequentist_vcov(sandwich = hc4_sandwich()),
             readable = FALSE)
  expect_equal(unname(v2), unname(tb(e^2 / (1 - h))), tolerance = 1e-6)
  expect_equal(unname(v3), unname(tb(e^2 / (1 - h)^2)), tolerance = 1e-6)
  expect_equal(unname(v4),
               unname(tb(e^2 / (1 - h)^pmin(4, h / mean(h)))), tolerance = 1e-6)
})

clustered <- function(sizes, seed = 3) {
  set.seed(seed)
  id <- rep(seq_along(sizes), sizes)
  n <- length(id)
  u <- rnorm(length(sizes), sd = 0.6)[id]
  dd <- data.frame(id = id, x = runif(n), z = rnorm(n))
  dd$y <- 1 + 1.5 * dd$x + 0.5 * dd$z + u + rnorm(n, sd = 0.5)
  dd
}

cr_hand <- function(m, id, p) {
  X <- model.matrix(m)
  e <- unname(resid(m))
  A <- solve(crossprod(X))
  Hm <- X %*% A %*% t(X)
  meat <- 0
  for (g in unique(id)) {
    ix <- which(id == g)
    Fg <- mpow(Hm[ix, ix, drop = FALSE], p)
    ug <- t(X[ix, , drop = FALSE]) %*% (Fg %*% e[ix])
    meat <- meat + ug %*% t(ug)
  }
  A %*% meat %*% A
}

test_that("CR2 and CR3 of a clustered linear model are the cluster-block formulas", {
  dd <- clustered(c(3, 5, 4, 6, 2, 5, 4, 7, 3, 5))
  f <- statmod(y ~ x + z, distributions7::gaussian1_distrib(), dd)
  m <- lm(y ~ x + z, dd)
  v2 <- vcov(f, vcov_type = fvx(hc2_sandwich(cluster = ~ id)), param = "mu",
             readable = FALSE)
  v3 <- vcov(f, vcov_type = fvx(hc3_sandwich(cluster = ~ id)), param = "mu",
             readable = FALSE)
  expect_equal(unname(v2), unname(cr_hand(m, dd$id, -1 / 2)), tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_equal(unname(v3), unname(cr_hand(m, dd$id, -1)), tolerance = 1e-8,
               ignore_attr = TRUE)
  # one observation for each cluster reduces CR2 / CR3 to HC2 / HC3
  dd$one <- seq_len(nrow(dd))
  a <- vcov(f, vcov_type = fvx(hc2_sandwich(cluster = ~ one, data = dd)),
            param = "mu", readable = FALSE)
  b <- vcov(f, vcov_type = fvx(hc2_sandwich()), param = "mu", readable = FALSE)
  expect_equal(unname(a), unname(b), tolerance = 1e-9, ignore_attr = TRUE)
  a3 <- vcov(f, vcov_type = fvx(hc3_sandwich(cluster = ~ one, data = dd)),
             param = "mu", readable = FALSE)
  b3 <- vcov(f, vcov_type = fvx(hc3_sandwich()), param = "mu", readable = FALSE)
  expect_equal(unname(a3), unname(b3), tolerance = 1e-9, ignore_attr = TRUE)
})

# the Satterthwaite degrees of freedom of Pustejovsky and Tipton (2018) for
# the contrast c, written with the full n x n hat matrix and a working model
# of independent errors: q_g = (I - H) E_g A_g X_g (X'X)^-1 c
pt_df <- function(m, id, cvec) {
  X <- model.matrix(m)
  n <- nrow(X)
  A <- solve(crossprod(X))
  Hm <- X %*% A %*% t(X)
  Q <- NULL
  for (g in unique(id)) {
    ix <- which(id == g)
    Fg <- mpow(Hm[ix, ix, drop = FALSE], -1 / 2)
    tg <- Fg %*% (X[ix, , drop = FALSE] %*% (A %*% cvec))
    Eg <- numeric(n)
    Eg[ix] <- tg
    Q <- cbind(Q, (diag(n) - Hm) %*% Eg)
  }
  Gm <- crossprod(Q)
  sum(diag(Gm))^2 / sum(Gm^2)
}

test_that("the CR2 degrees of freedom are those of Pustejovsky and Tipton", {
  for (sizes in list(rep(4, 12), c(3, 5, 4, 6, 2, 5, 4, 7, 3, 5, 6, 2))) {
    dd <- clustered(sizes, seed = 9)
    f <- statmod(y ~ x + z, distributions7::gaussian1_distrib(), dd)
    m <- lm(y ~ x + z, dd)
    V <- vcov(f, vcov_type = fvx(hc2_sandwich(cluster = ~ id)),
              readable = FALSE)
    fn <- attr(V, "df_fn")
    expect_true(is.function(fn))
    got <- fn(c("mu:(Intercept)", "mu:x", "mu:z"))
    want <- c(pt_df(m, dd$id, c(1, 0, 0)), pt_df(m, dd$id, c(0, 1, 0)),
              pt_df(m, dd$id, c(0, 0, 1)))
    expect_equal(got, want, tolerance = 1e-7)
    expect_true(all(got > 0 & got < length(sizes)))
  }
})

test_that("a balanced two-group comparison of clusters has G - 2 degrees of freedom", {
  set.seed(12)
  G <- 14
  id <- rep(seq_len(G), each = 5)
  grp <- as.numeric(id > G / 2)
  dd <- data.frame(id = id, g = grp)
  dd$y <- 0.3 * dd$g + rnorm(G)[id] + rnorm(length(id), sd = 0.5)
  f <- statmod(y ~ g, distributions7::gaussian1_distrib(), dd)
  V <- vcov(f, vcov_type = fvx(hc2_sandwich(cluster = ~ id)), readable = FALSE)
  expect_equal(attr(V, "df_fn")("mu:g"), G - 2, tolerance = 1e-8)
})

test_that("confint, summary and statmod_test read the t reference under CR2", {
  dd <- clustered(c(3, 5, 4, 6, 2, 5, 4, 7), seed = 21)
  f <- statmod(y ~ x + z, distributions7::gaussian1_distrib(), dd)
  vt <- fvx(hc2_sandwich(cluster = ~ id))
  vn <- fvx(hc2_sandwich(cluster = ~ id, reference = "normal"))
  V <- vcov(f, vcov_type = vt, readable = FALSE)
  dfs <- attr(V, "df_fn")(c("mu:(Intercept)", "mu:x", "mu:z"))
  ci <- confint(f, vcov_type = vt, readable = FALSE)
  est <- ci$estimate[1:3]
  se <- sqrt(diag(V))[1:3]
  expect_equal(unname(ci$upper[1:3]), unname(est + qt(0.975, dfs) * se),
               tolerance = 1e-9)
  # the normal reference gives narrower intervals
  cn <- confint(f, vcov_type = vn, readable = FALSE)
  expect_equal(unname(cn$upper[1:3]), unname(est + qnorm(0.975) * se),
               tolerance = 1e-9)
  expect_true(all(ci$upper[1:3] - ci$lower[1:3] > cn$upper[1:3] - cn$lower[1:3]))
  # other sandwich types stay on the normal
  c3 <- confint(f, vcov_type = fvx(hc3_sandwich(cluster = ~ id)),
                readable = FALSE)
  v3 <- vcov(f, vcov_type = fvx(hc3_sandwich(cluster = ~ id)), readable = FALSE)
  expect_equal(c3$upper[2], est[2] + qnorm(0.975) * sqrt(v3[2, 2]),
               tolerance = 1e-9)
  # summary: the p-values are t
  s <- summary(f, vcov_type = vt)
  expect_no_error(capture.output(print(s)))
  st <- statmod_test(f, param = "mu", coefname = "x", vcov_type = vt)
  zst <- (est[2]) / se[2]
  expect_equal(unname(st@p.value), unname(2 * pt(-abs(zst), dfs[2])),
               tolerance = 1e-8)
  stn <- statmod_test(f, param = "mu", coefname = "x", vcov_type = vn)
  expect_equal(unname(stn@p.value), unname(2 * pnorm(-abs(zst))),
               tolerance = 1e-8)
})

test_that("the constructors check their arguments", {
  expect_error(hc4_sandwich(cluster = ~ id), "unused argument")
  expect_error(hc2_sandwich(reference = "z"), "should be one of")
  expect_equal(hc3_sandwich(cluster = ~ id)@type, "HC3")
  expect_error(SandwichMethod(type = "HC4", cluster = ~ id), "clusters")
})

test_that("CR2 on a fit with a random intercept and n = 2000 is timed", {
  skip_on_cran()
  set.seed(31)
  G <- 100
  id <- rep(seq_len(G), each = 20)
  dd <- data.frame(id = id, x = rnorm(2000))
  dd$y <- 1 + 0.5 * dd$x + rnorm(G, sd = 0.7)[id] + rnorm(2000)
  f <- statmod(y ~ x + random(~ 1 | id), distributions7::gaussian1_distrib(),
               dd)
  t0 <- system.time(
    V <- vcov(f, vcov_type = fvx(hc2_sandwich(cluster = ~ id)),
              param = "mu", readable = FALSE))[["elapsed"]]
  expect_true(all(is.finite(diag(V))))
  expect_lt(t0, 60)
  message(sprintf("CR2 variance, n = 2000, 100 clusters: %.2f s", t0))
})

test_that("integer prior weights equal the data expanded to copies", {
  set.seed(14)
  dd <- clustered(c(4, 3, 5, 4, 6, 3), seed = 14)
  dd$w <- sample(1:3, nrow(dd), replace = TRUE)
  ex <- dd[rep(seq_len(nrow(dd)), dd$w), ]
  fw <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd, weights = dd$w)
  fe <- statmod(y ~ x, distributions7::gaussian1_distrib(), ex)
  for (sw in list(hc2_sandwich(), hc3_sandwich(), hc2_sandwich(cluster = ~ id),
                  hc3_sandwich(cluster = ~ id))) {
    vw <- vcov(fw, vcov_type = fvx(sw), param = "mu", readable = FALSE)
    ve <- vcov(fe, vcov_type = fvx(sw), param = "mu", readable = FALSE)
    expect_equal(unname(vw), unname(ve), tolerance = 1e-7, ignore_attr = TRUE)
  }
})
