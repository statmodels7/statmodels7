test_that("a score loading has a second start on the scale of its score", {
  # The score is unscaled, so a loading of 0.1 is a weak response only where
  # the score is of order one. On the share of front-seat casualties in
  # Seatbelts a beta on the logit has a score of order phi, near 130, and the
  # filter exploded at the first step: the fit stopped, unconverged, at
  # -3326.66 where the static model reaches 413.58. The second start, 0.1/I,
  # is fitted from where the first fails.
  sb <- data.frame(month = seq_len(nrow(datasets::Seatbelts)),
                   share = datasets::Seatbelts[, "front"] /
                     (datasets::Seatbelts[, "front"] +
                        datasets::Seatbelts[, "rear"]),
                   law = datasets::Seatbelts[, "law"])
  spec <- statmod_spec(share ~ law + gas(p = 1, q = 1, time = month),
                       distributions7::beta1_distrib(), sb)
  info <- predictor_information(spec)
  st <- statmod_structural_state(statmod_design(spec))
  # the first start is the term's own, the second is 0.1/I
  expect_equal(exp(st$zeta[[1L]][["kappa1"]]), 0.1, tolerance = 1e-12)
  expect_equal(exp(st$info_start[[1L]][["kappa1"]]), 0.1 / info[["mu"]],
               tolerance = 1e-12)
  fit <- statmod(share ~ law + gas(p = 1, q = 1, time = month),
                 distributions7::beta1_distrib(), sb)
  fit0 <- statmod(share ~ law, distributions7::beta1_distrib(), sb)
  expect_identical(statmod_certificate(fit)$state, "converged")
  expect_gt(as.numeric(logLik(fit)), as.numeric(logLik(fit0)))

  # for a gaussian mean the information is 1/sigma^2, so the second start of
  # the loading is a tenth of the variance
  dd <- data.frame(t = 1:50, y = 900 + 100 * stats::rnorm(50))
  sp2 <- statmod_spec(y ~ gas(p = 1, q = 1, time = t),
                      distributions7::gaussian1_distrib(), dd)
  s0 <- exp(statmod_intercepts(sp2)$sigma)
  expect_equal(predictor_information(sp2)[["mu"]], 1 / s0^2,
               tolerance = 1e-10)
})

test_that("a free level fits the same model as the intercept spelling", {
  # flow ~ 0 + gas() searched omega and the persistence on the ridge
  # omega = L (1 - sum phi) and stopped after 200 Newton iterations at
  # -642.67; the joint fit now searches the level L itself
  nile <- data.frame(year = 1871:1970, flow = as.numeric(Nile))
  for (q in 1:2) {
    a <- statmod(flow ~ gas(p = 1, q = q, time = year),
                 distrib = distributions7::gaussian1_distrib(), data = nile)
    b <- statmod(flow ~ 0 + gas(p = 1, q = q, time = year),
                 distrib = distributions7::gaussian1_distrib(), data = nile)
    expect_equal(as.numeric(logLik(b)), as.numeric(logLik(a)), tolerance = 1e-8)
    expect_identical(statmod_certificate(b)$state, "converged")
    ph <- coef(a)$mu[grep("xi", names(coef(a)$mu))]
    expect_equal(coef(b)$mu[["gas.omega"]],
                 coef(a)$mu[["(Intercept)"]] * (1 - sum(ph)), tolerance = 1e-5)
  }
})

test_that("the level chart's derivatives are those of the composed objective", {
  nile <- data.frame(year = 1871:1970, flow = as.numeric(Nile))
  b <- statmod(flow ~ 0 + gas(p = 1, q = 2, time = year),
               distrib = distributions7::gaussian1_distrib(), data = nile)
  spec <- b@spec
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, b@hyper, design, FALSE, "bartlett")
  jp <- statmod_joint_pieces(spec, design, obj, b@hyper)
  lv <- joint_level_chart(spec, jp)
  expect_false(is.null(lv))
  u <- c(obj$stack(b@coefficients), jp$zeta()) + c(0.01, 3, 0.1, 0.2, -0.1)
  v <- lv$from(u)
  expect_equal(lv$to(v), u)
  f <- function(v) jp$fn(lv$to(v))
  g <- lv$grad(v, jp$gr(lv$to(v)))
  h <- lv$hess(v, jp$gr(lv$to(v)), jp$he(lv$to(v)))
  expect_equal(g, numDeriv::grad(f, v), tolerance = 1e-6, ignore_attr = TRUE)
  k <- 2:5
  expect_equal(h[k, k], numDeriv::hessian(f, v)[k, k], tolerance = 1e-5)
})
