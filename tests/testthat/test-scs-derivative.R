# plot_derivative() of modelterms7::scs(): the curve against a central
# difference of predict(), its sign under a monotone and a concave shape, and
# a shape named by an object in the call.

test_that("plot_derivative() of scs() is the derivative of the fitted effect", {
  set.seed(51)
  d <- data.frame(x = runif(200))
  d$y <- 1 - exp(-3 * d$x) + rnorm(200, sd = 0.08)
  sh <- modelterms7::shape_increasing_concave()
  fit <- statmod(y ~ modelterms7::scs(x, sh),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  nm <- names(fit@spec@terms$mu)[vapply(fit@spec@terms$mu, function(t)
    S7::S7_inherits(t, modelterms7::ScsTerm), logical(1))]
  mu <- function(x) predict(fit, what = "mu", newdata = data.frame(x = x))
  for (ord in 1:2) {
    p <- autoplot(fit, term = nm, param = "mu",
                  plot_type = plot_derivative(order = ord), data = d)
    cur <- NULL
    for (L in p$layers) {
      if (is.data.frame(L$data) && all(c("x", "fit", "se") %in% names(L$data))) cur <- L$data
    }
    if (is.null(cur)) cur <- p$data
    expect_true(all(is.finite(cur$fit)) && all(cur$se >= 0))
    h <- 1e-4
    xs <- cur$x[cur$x > min(d$x) + 0.01 & cur$x < max(d$x) - 0.01]
    num <- if (ord == 1) (mu(xs + h) - mu(xs - h)) / (2 * h)
           else (mu(xs + h) - 2 * mu(xs) + mu(xs - h)) / h^2
    est <- cur$fit[match(xs, cur$x)]
    expect_lt(max(abs(est - num)) / max(abs(num)), 1e-3)
    if (ord == 1) expect_true(all(cur$fit >= -1e-8))
    if (ord == 2) expect_true(all(cur$fit <= 1e-6))
  }
})
