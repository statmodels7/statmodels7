# A small fit with a linear part, a factor, a smooth and a random intercept
# in mu and a modeled sigma, shared by the tests of the graphs.
plot_test_fit <- function(n = 200, seed = 1) {
  set.seed(seed)
  dd <- data.frame(x = runif(n), z = rnorm(n),
                   g = factor(sample(c("a", "b", "c"), n, TRUE)),
                   id = factor(sample(1:12, n, TRUE)))
  dd$y <- sin(2 * pi * dd$x) + 0.5 * dd$z + c(0, 0.4, -0.3)[dd$g] +
    stats::rnorm(12, sd = 0.5)[dd$id] + stats::rnorm(n, sd = exp(-1 + dd$x))
  fit <- statmod(y ~ z + g + s(x) + random(~ 1 | id) | sigma ~ x,
                 distributions7::gaussian1_distrib(), dd)
  list(fit = fit, data = dd)
}
