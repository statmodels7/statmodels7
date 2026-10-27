g1 <- distributions7::gaussian1_distrib()

fsw_w <- function(sw) frequentist_vcov(sandwich = sw)

sw_data <- function(n = 120, seed = 5) {
  set.seed(seed)
  d <- data.frame(x = runif(n), id = rep(seq_len(n / 4), each = 4))
  d$y <- 1 + 2 * d$x + rnorm(n, sd = 0.2 + d$x)
  d$w <- exp(rnorm(n, sd = 0.6))
  d
}

test_that("the default reading of the weights is frequency", {
  expect_identical(hc0_sandwich()@weights, "frequency")
  expect_identical(hc1_sandwich(weights = "sampling")@weights, "sampling")
  expect_identical(hc2_sandwich(weights = "sampling")@weights, "sampling")
  expect_identical(hc3_sandwich(weights = "sampling")@weights, "sampling")
  expect_identical(hc4_sandwich(weights = "sampling")@weights, "sampling")
  expect_error(hc0_sandwich(weights = "survey"), "should be one of")
  expect_output(print(hc1_sandwich(weights = "sampling")),
                "sampling weights")
  d <- sw_data()
  f <- statmod(y ~ x, g1, d, weights = d$w)
  expect_identical(
    vcov(f, vcov_type = fsw_w(hc1_sandwich())),
    vcov(f, vcov_type = fsw_w(hc1_sandwich(weights = "frequency"))))
})

test_that("sampling weights equal frequency weights when every weight is one", {
  d <- sw_data()
  f <- statmod(y ~ x, g1, d)
  for (mk in list(hc0_sandwich, hc1_sandwich, hc2_sandwich, hc3_sandwich,
                  hc4_sandwich)) {
    a <- vcov(f, vcov_type = fsw_w(mk()), readable = FALSE)
    b <- vcov(f, vcov_type = fsw_w(mk(weights = "sampling")),
              readable = FALSE)
    expect_equal(a, b, tolerance = 1e-12)
  }
  for (mk in list(function(w) hc0_sandwich(cluster = ~ id, weights = w),
                  function(w) hc2_sandwich(cluster = ~ id, weights = w))) {
    a <- vcov(f, vcov_type = fsw_w(mk("frequency")), readable = FALSE)
    b <- vcov(f, vcov_type = fsw_w(mk("sampling")), readable = FALSE)
    expect_equal(a, b, tolerance = 1e-12)
  }
})

test_that("the sampling meat is the sum of w^2 g g' on a Gaussian mean", {
  d <- sw_data()
  n <- nrow(d)
  f <- statmod(y ~ x, g1, d, weights = d$w, outer_criterion = ml())
  X <- cbind(1, d$x)
  w <- d$w
  b <- drop(solve(crossprod(X, X * w), crossprod(X, w * d$y)))
  expect_equal(unname(f@coefficients$mu), b, tolerance = 1e-6)
  r <- drop(d$y - X %*% b)
  A <- solve(crossprod(X, X * w))
  hc0 <- A %*% crossprod(X * (w * r)) %*% A
  v0 <- vcov(f, vcov_type = fsw_w(hc0_sandwich(weights = "sampling")),
             param = "mu", readable = FALSE)
  v1 <- vcov(f, vcov_type = fsw_w(hc1_sandwich(weights = "sampling")),
             param = "mu", readable = FALSE)
  expect_equal(unname(v0), unname(hc0), tolerance = 1e-7)
  # n counts the rows and p is the effective degrees of freedom of the model
  p <- sum(f@edf$edf)
  expect_equal(unname(v1), unname(hc0 * n / (n - p)), tolerance = 1e-7)
  # the frequency reading is the other meat
  vf <- vcov(f, vcov_type = fsw_w(hc0_sandwich()), param = "mu",
             readable = FALSE)
  freq <- A %*% crossprod(X * sqrt(w) * r) %*% A
  expect_equal(unname(vf), unname(freq), tolerance = 1e-7)
  expect_gt(max(abs(vf - v0)), 1e-3)
})

test_that("the sampling variance does not move when the weights are rescaled", {
  d <- sw_data()
  f1 <- statmod(y ~ x, g1, d, weights = d$w, outer_criterion = ml())
  f2 <- statmod(y ~ x, g1, d, weights = 7 * d$w, outer_criterion = ml())
  for (mk in list(hc0_sandwich, hc1_sandwich, hc2_sandwich, hc3_sandwich,
                  hc4_sandwich)) {
    a <- vcov(f1, vcov_type = fsw_w(mk(weights = "sampling")),
              readable = FALSE)
    b <- vcov(f2, vcov_type = fsw_w(mk(weights = "sampling")),
              readable = FALSE)
    expect_equal(a, b, tolerance = 1e-6)
  }
  # the frequency reading moves with the scale of the weights
  a <- vcov(f1, vcov_type = fsw_w(hc0_sandwich()), readable = FALSE)
  b <- vcov(f2, vcov_type = fsw_w(hc0_sandwich()), readable = FALSE)
  expect_equal(a, 7 * b, tolerance = 1e-6)
})

test_that("with a cluster the two readings share the meat and differ in n", {
  d <- sw_data()
  f <- statmod(y ~ x, g1, d, weights = d$w, outer_criterion = ml())
  a <- vcov(f, vcov_type = fsw_w(hc0_sandwich(cluster = ~ id)),
            readable = FALSE)
  b <- vcov(f, vcov_type = fsw_w(hc0_sandwich(cluster = ~ id,
                                              weights = "sampling")),
            readable = FALSE)
  expect_equal(a, b, tolerance = 1e-12)
  a1 <- vcov(f, vcov_type = fsw_w(hc1_sandwich(cluster = ~ id)),
             readable = FALSE)
  b1 <- vcov(f, vcov_type = fsw_w(hc1_sandwich(cluster = ~ id,
                                               weights = "sampling")),
             readable = FALSE)
  n <- nrow(d)
  p <- sum(f@edf$edf)
  ratio <- ((n - 1) / (n - p)) / ((sum(d$w) - 1) / (sum(d$w) - p))
  expect_equal(b1, a1 * ratio, tolerance = 1e-10)
})

test_that("a row with weight zero is not counted under the sampling reading", {
  d <- sw_data()
  d$w[c(3, 17, 40)] <- 0
  f <- statmod(y ~ x, g1, d, weights = d$w, outer_criterion = ml())
  fs <- statmod(y ~ x, g1, d[d$w > 0, ], weights = w, outer_criterion = ml())
  for (mk in list(hc1_sandwich, hc3_sandwich)) {
    a <- vcov(f, vcov_type = fsw_w(mk(weights = "sampling")),
              readable = FALSE)
    b <- vcov(fs, vcov_type = fsw_w(mk(weights = "sampling")),
              readable = FALSE)
    expect_equal(a, b, tolerance = 1e-6)
  }
})
