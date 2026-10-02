skip_on_cran()

# the lasso objective over an intercept and every level of g, solved exactly:
# for a fixed intercept each level is a soft threshold of its group mean, and
# the intercept is a one-dimensional minimization
lasso_levels_by_hand <- function(y, g, sig, lam) {
  ng <- as.numeric(table(g))
  yb <- tapply(y, g, mean)
  obj <- function(p) sum((y - p[1] - p[-1][g])^2) / (2 * sig^2) + lam * sum(abs(p[-1]))
  bb <- function(c0) {
    z <- yb - c0
    sign(z) * pmax(abs(z) - lam * sig^2 / ng, 0)
  }
  o <- optimize(function(c0) obj(c(c0, bb(c0))), range(y), tol = 1e-12)
  list(par = unname(c(o$minimum, bb(o$minimum))), obj = obj)
}

test_that("a lasso over every level beside an intercept reports every level", {
  set.seed(1)
  g <- factor(rep(LETTERS[1:8], each = 12))
  for (vm in list(c(200, 200, 260, 200, 200, 150, 200, 200),
                  c(200, 200, 200, 200, 200, 200, 200, 260))) {
    y <- vm[g] + rnorm(96, sd = 10)
    f <- statmod(y ~ 1 + lasso(~ g), distrib = distributions7::gaussian1_distrib(),
                 data = data.frame(g = g, y = y))
    # the level a pivot on the whole information named is identified by the
    # kink, and nothing is aliased
    expect_length(f@aliased, 0L)
    cf <- coef(f)$mu
    expect_false(anyNA(cf))
    ref <- lasso_levels_by_hand(y, g, exp(coef(f)$sigma[[1]]), hyper(f)$estimate)
    expect_equal(unname(cf), ref$par, tolerance = 1e-6)
    expect_equal(ref$obj(unname(cf)), ref$obj(ref$par), tolerance = 1e-9)
    # the summary prints, also where the last level is the only one selected
    expect_no_error(capture.output(print(summary(f))))
  }
})
