# SCAD and MCP are scaled by the curvature of the likelihood in each of their
# coordinates (decision 4 of 2026-09-27), the penalty being
# sum_j c_j rho(b_j; lambda / c_j), so the knee sits at a lambda / c_j and the
# shape parameter means the same thing whatever the data's scale. Since
# 0.162.0 (decided 2026-09-29) the curvature is SELF-CONSISTENT: it is the one
# the coordinate descent steps with, the centered weighted sum of squares of
# each column, read again at every pass and at the fit.

set.seed(11)
xc <- matrix(stats::rnorm(120 * 6, sd = rep(c(0.5, 1, 4), each = 240)),
             120, 6)
dc <- data.frame(y = 2 + xc[, 1] - xc[, 4] + stats::rnorm(120, sd = 1.5))
dc$x <- xc

curv_of_fit <- function(fit, fam) {
  tn <- grep(fam, names(fit@spec@terms$mu), fixed = TRUE, value = TRUE)
  expect_length(tn, 1L)
  modelterms7::term_penalty(fit@spec@terms$mu[[tn]])
}

centered_ss <- function(x) colSums(sweep(x, 2L, colMeans(x))^2)

test_that("a gaussian mean's curvature is the centered sum of squares over the
          fitted variance", {
  # the reference is the information of a gaussian mean with its intercept
  # profiled out, written by hand and read at the fitted variance; it shares
  # no arithmetic with coord_colsq() or coord_working()
  for (f in list(y ~ scad(x, lambda = 3, a = 3.7),
                 y ~ mcp(x, lambda = 3, gamma = 3))) {
    fit <- statmod(f, distributions7::gaussian1_distrib(), dc)
    pen <- curv_of_fit(fit, if (grepl("scad", deparse(f))) "scad" else "mcp")
    s <- fit@fitted$sigma[1L]
    expect_equal(pen@curv, centered_ss(xc) / s^2, tolerance = 1e-8)
  }
})

test_that("under standardize the curvature is on the standardized scale", {
  fit <- statmod(y ~ scad(x, lambda = 3, a = 3.7, standardize = TRUE),
                 distributions7::gaussian1_distrib(), dc)
  pen <- curv_of_fit(fit, "scad")
  d <- as.numeric(Matrix::diag(pen@map))
  s <- fit@fitted$sigma[1L]
  expect_equal(pen@curv, centered_ss(xc) / s^2 / d^2, tolerance = 1e-8)
})

test_that("columns far from centered fit, the curvature being centered", {
  # With means of 30 standard deviations the uncentered sum of squares is
  # about 900 times the centered one the step works with, and a curvature
  # read uncentered failed the step condition on every fit (0.159.0 to
  # 0.161.0). The self-consistent curvature makes the scaled step one.
  set.seed(12)
  xr <- sweep(matrix(stats::rnorm(150 * 5), 150, 5), 2L, 30, "+")
  dr <- data.frame(y = 1 + xr[, 1] - xr[, 2] + stats::rnorm(150))
  dr$x <- xr
  for (f in list(y ~ scad(x), y ~ mcp(x))) {
    fit <- statmod(f, distributions7::gaussian1_distrib(), dr)
    expect_true(fit@converged)
    pen <- curv_of_fit(fit, if (grepl("scad", deparse(f))) "scad" else "mcp")
    s <- fit@fitted$sigma[1L]
    expect_equal(pen@curv, centered_ss(xr) / s^2, tolerance = 1e-8)
  }
})

test_that("a pre-scaled design and standardize are one model", {
  # the penalty is written on u = D beta either way, and the curvature on u
  # is the same number, so at one lambda the two fits agree
  z <- scale(xc)
  dz <- dc
  dz$z <- z
  a <- statmod(y ~ scad(x, lambda = 3, a = 3.7, standardize = TRUE),
               distributions7::gaussian1_distrib(), dc)
  b <- statmod(y ~ scad(z, lambda = 3, a = 3.7),
               distributions7::gaussian1_distrib(), dz)
  expect_equal(fitted(a), fitted(b), tolerance = 1e-6)
  expect_equal(unname(a@coefficients$mu[-1L]) * attr(z, "scaled:scale"),
               unname(b@coefficients$mu[-1L]), tolerance = 1e-6)
})

test_that("a lasso and a ridge are not touched", {
  for (f in list(y ~ lasso(x, lambda = 3), y ~ ridge(x, lambda = 3))) {
    fit <- statmod(f, distributions7::gaussian1_distrib(), dc)
    fam <- if (grepl("lasso", deparse(f))) "lasso" else "ridge"
    pen <- curv_of_fit(fit, fam)
    expect_false("curv" %in% S7::prop_names(pen))
  }
})

test_that("under standardize the fit does not depend on the units", {
  # a covariate rescaled by k has its coefficient divided by k, its spread
  # multiplied by k and its curvature by k^2, so on the standardized scale
  # both the slope and the knee are unchanged and so is the fit. Without
  # standardize the lasso part lambda |b| is not unit-free, as for a lasso.
  d2 <- dc
  d2$x <- xc * 10
  f <- y ~ scad(x, lambda = 3, a = 3.7, standardize = TRUE)
  a <- statmod(f, distributions7::gaussian1_distrib(), dc)
  b <- statmod(f, distributions7::gaussian1_distrib(), d2)
  expect_equal(fitted(a), fitted(b), tolerance = 1e-6)
  expect_equal(a@coefficients$mu[-1L], 10 * b@coefficients$mu[-1L],
               tolerance = 1e-6)
})

test_that("the table is built for the kept coordinates", {
  # the kernel reads row a for coordinate keep[a]; the restricted penalty's
  # rows are the full penalty's kept rows, for a map and for a curvature
  D <- Matrix::Diagonal(x = c(1, 2, 3, 4))
  st <- c(0.1, 0.2, 0.3, 0.4)
  keep <- c(2L, 4L)
  v <- 1 / st
  for (pen in list(penalties7::lasso_penalty(map = D),
                   penalties7::mcp_penalty(n_coef = 4L))) {
    th <- stats::setNames(as.list(rep(1, length(pen@params))), pen@params)
    if ("gamma" %in% pen@params) th$gamma <- 3
    full <- penalties7::penalty_prox_spec(pen, th, st)
    sub <- penalties7::penalty_prox_spec(coord_table_penalty(pen, keep, v), th,
                                         st[keep])
    expect_equal(sub$cut, full$cut[keep, , drop = FALSE])
    expect_equal(sub$icept, full$icept[keep, , drop = FALSE])
  }
  # the call as it was: the whole penalty with the kept steps reads
  # coordinates 1 and 2 where 2 and 4 are meant
  pen <- penalties7::lasso_penalty(map = D)
  old <- penalties7::penalty_prox_spec(pen, list(lambda = 1), st[keep])
  full <- penalties7::penalty_prox_spec(pen, list(lambda = 1), st)
  expect_false(isTRUE(all.equal(old$cut[1:2, 1], full$cut[keep, 1])))
})

test_that("a scaled table takes the curvature of the current step", {
  # c_j = v_j on the scale of u, so the scaled step is one and the step
  # condition holds for every admissible shape
  D <- Matrix::Diagonal(x = c(1, 2, 3, 4))
  pen <- penalties7::scad_penalty(map = D, curv = c(9, 9, 9, 9))
  v <- c(5, 7, 11, 13)
  keep <- c(1L, 3L, 4L)
  tp <- coord_table_penalty(pen, keep, v)
  expect_equal(tp@curv, v[keep] / c(1, 3, 4)^2)
  tab <- penalties7::penalty_prox_spec(tp, list(lambda = 1, a = 2.05),
                                       1 / v[keep])
  expect_false(is.null(tab))
})

test_that("every coordinate is screened against its own kink", {
  D <- Matrix::Diagonal(x = c(0.5, 1, 4))
  pen <- penalties7::lasso_penalty(map = D)
  k <- coord_kinks(pen, list(lambda = 2))
  expect_equal(k, 2 * c(0.5, 1, 4), tolerance = 1e-8)
  expect_equal(k[[1L]], kink_scale(pen, list(lambda = 2)))
})

test_that("the fitted coefficients are a KKT point of the penalty at the fit's
          curvature", {
  # What self-consistency means: the coefficients are optimal for the
  # penalty whose curvature is the one at those coefficients. The derivative
  # of c * scad(b; lambda / c) is written out by hand. Strong effects put two
  # coefficients between the knees, where the gradient reads the curvature:
  # undamped, the first alternated between 1.009 and 1.496 and the fit did
  # not converge, the KKT residual reading 2.28; with the curvature read
  # once at the start the step condition failed outright.
  set.seed(15)
  x <- matrix(stats::rnorm(200 * 8), 200, 8)
  d <- data.frame(y = stats::rbinom(200, 1, stats::plogis(
    2.5 * x[, 1] - 2 * x[, 2] + x[, 3])))
  d$x <- x
  lam <- 16
  a <- 3.7
  fit <- statmod(y ~ scad(x, lambda = lam, a = a),
                 distributions7::bernoulli_distrib(), d)
  expect_true(fit@converged)
  pen <- curv_of_fit(fit, "scad")
  cf <- fit@coefficients$mu
  b <- unname(cf[-1L])
  mu <- stats::plogis(cf[[1L]] + drop(x %*% b))
  g <- drop(crossprod(x, d$y - mu))
  dscad <- function(t, l) {
    t <- abs(t)
    ifelse(t <= l, l, ifelse(t <= a * l, (a * l - t) / (a - 1), 0))
  }
  act <- b != 0
  expect_true(any(act))
  mid <- abs(b) > lam / pen@curv & abs(b) < a * lam / pen@curv
  expect_true(any(mid))
  lj <- lam / pen@curv
  expect_equal(g[act], pen@curv[act] * dscad(b[act], lj[act]) * sign(b[act]),
               tolerance = 1e-6)
  expect_true(all(abs(g[!act]) <= lam * (1 + 1e-6)))
})

test_that("the table's curvature is damped towards the current step's", {
  pen <- penalties7::scad_penalty(n_coef = 3L, curv = c(4, 4, 4))
  tp <- coord_table_penalty(pen, 1:3, c(1, 4, 16), prev = c(4, 4, 4))
  expect_equal(tp@curv, c(2, 4, 8))
})
