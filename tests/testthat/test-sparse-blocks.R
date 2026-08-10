# A penalty with a kink is fitted by a scheme of its own, and what decides is
# penalty_kinks(). Everything here is a defect y ~ scad(x) exposed.

set.seed(1)
xs <- matrix(stats::rnorm(100 * 20), 100, 20)
# a matrix column, which is how a penalized term reads a design: data.frame()
# would split the matrix into twenty columns and leave `x` to be found in the
# calling environment instead
ds <- data.frame(y = stats::rnorm(100))
ds$x <- xs

# rho' for scad and mcp, transcribed from the published forms: the KKT check
# below has to share no arithmetic with penalties7
scad_deriv <- function(b, lam, a) {
  u <- abs(b)
  sign(b) * ifelse(u <= lam, lam,
                   ifelse(u <= a * lam, (a * lam - u) / (a - 1), 0))
}
mcp_deriv <- function(b, lam, g) {
  u <- abs(b)
  sign(b) * ifelse(u <= g * lam, lam - u / g, 0)
}

test_that("a kinked penalty is fitted by the proximal scheme", {
  # penalty_kinks() stopped on the shape of theta, the tryCatch read the
  # failure as "no kink", and scad and mcp went into the jointly fitted
  # system -- solved by the curvature of a function that has none. On this
  # design of pure noise the fit kept 19.00 effective degrees of freedom out
  # of 20, which is no selection at all, and reported convergence FALSE.
  for (f in list(y ~ scad(x), y ~ mcp(x), y ~ lasso(x))) {
    spec <- statmod_spec(f, distributions7::gaussian1_distrib(), ds)
    blocks <- statmod_blocks(spec, statmod_design(spec))
    expect_length(blocks$sparse, 1L)
  }
  # while a penalty that is twice differentiable stays in the joint system
  spec <- statmod_spec(y ~ ridge(x), distributions7::gaussian1_distrib(), ds)
  expect_length(statmod_blocks(spec, statmod_design(spec))$sparse, 0L)
})

test_that("a penalty that cannot answer is reported, not assumed smooth", {
  # reading an error as "no kink" assigns the term to the scheme for the
  # opposite property, silently
  Mute <- S7::new_class("Mute", parent = penalties7::penalty,
                        constructor = function() {
                          S7::new_object(S7::S7_object(),
                                         penalty_name = "mute", map = NULL,
                                         n_coef = 3, params = "lambda",
                                         params_bounds = list(lambda = c(0, Inf)),
                                         link_params = list(lambda = NULL),
                                         params_smooth = c(lambda = TRUE))
                        })
  expect_error(penalty_has_kink(Mute(), "The penalty of q(x)"),
               "cannot say whether it is differentiable")
  expect_error(penalty_has_kink(Mute(), "The penalty of q(x)"),
               "The penalty of q\\(x\\)")
})

test_that("scad and mcp reach a point that satisfies the KKT conditions", {
  # the reference is the subdifferential of the objective, written out here:
  # at a coefficient away from zero the penalized score vanishes, and at one
  # that is exactly zero the unpenalized score is inside the interval the
  # penalty's kink opens
  lam <- 4
  cases <- list(
    list(f = y ~ scad(x), nm = "scad(x)", d = scad_deriv,
         th = c(lambda = lam, a = 3.7)),
    list(f = y ~ mcp(x), nm = "mcp(x)", d = mcp_deriv,
         th = c(lambda = lam, gamma = 3)))
  for (cs in cases) {
    hy <- list(mu = stats::setNames(list(cs$th), cs$nm))
    fit <- statmod(cs$f, distributions7::gaussian1_distrib(), ds, hyper = hy)
    expect_true(fit@converged)
    b <- fit@coefficients$mu
    s <- exp(fit@coefficients$sigma[1L])
    g <- -as.numeric(crossprod(xs, ds$y - cbind(1, xs) %*% b)) / s^2
    sl <- b[-1L]
    nz <- abs(sl) > 1e-8
    # the penalty selects: on a design of pure noise at this lambda most of
    # the twenty columns are dropped, where the joint system dropped none
    expect_lt(sum(nz), 20L)
    expect_lt(max(abs(g[nz] + cs$d(sl[nz], lam, cs$th[[2L]]))), 1e-5)
    if (any(!nz)) expect_lt(max(abs(g[!nz])), lam)
  }
})

test_that("the shrinkage answers the smoothing parameter", {
  # measured on this design: 17, 9, 0 surviving columns
  keep <- vapply(c(1, 5, 20), function(lam) {
    fit <- statmod(y ~ scad(x), distributions7::gaussian1_distrib(), ds,
                   hyper = list(mu = list("scad(x)" = c(lambda = lam,
                                                        a = 3.7))))
    sum(abs(fit@coefficients$mu[-1L]) > 1e-8)
  }, integer(1))
  expect_true(all(diff(keep) < 0))
  expect_identical(keep[[3L]], 0L)
})

test_that("the point beats ncvreg's on the objective they share", {
  # both are stationary points of a non-convex problem, so the comparison is
  # of the objective. Measured, ours is lower on both -- 52.9966 against
  # 53.5323 for scad and 52.9948 against 53.6638 for mcp -- and starting our
  # iteration FROM their point moves 6e-2 back to ours. The same shape as
  # ncvfit's MCP point in the modelterms7 comparison.
  skip_if_not_installed("ncvreg")
  n <- nrow(ds)
  xc <- scale(xs, center = TRUE, scale = FALSE)
  yc <- ds$y - mean(ds$y)
  dc <- data.frame(y = yc)
  dc$x <- xc
  lam <- 4
  for (cs in list(list(nm = "scad", sh = 3.7, key = "a"),
                  list(nm = "mcp", sh = 3, key = "gamma"))) {
    th <- stats::setNames(c(lam, cs$sh), c("lambda", cs$key))
    tn <- sprintf("%s(x)", cs$nm)
    f <- stats::as.formula(sprintf("y ~ %s(x) - 1 | sigma ~ 1", cs$nm))
    fit <- statmod(f, distributions7::gaussian1_distrib(), dc,
                   hyper = list(mu = stats::setNames(list(th), tn)))
    b <- fit@coefficients$mu
    s <- exp(fit@coefficients$sigma[1L])
    pen <- modelterms7::term_penalty(fit@spec@terms$mu[[tn]])
    obj <- function(bb) {
      sum((yc - xc %*% bb)^2) / (2 * s^2) +
        penalties7::penalty_value(pen, bb, as.list(th))
    }
    cv <- ncvreg::ncvfit(xc, yc, penalty = toupper(cs$nm),
                         lambda = lam * s^2 / n, gamma = cs$sh,
                         eps = 1e-12, max.iter = 1e6)
    expect_lt(obj(b), obj(as.numeric(cv$beta)))
  }
})
