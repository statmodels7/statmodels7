# The fast context a score-driven filter hands to modelterms7's kernel: the
# name distributions7's registry knows the distribution by and the constants
# its entries read after the parameters, both from distrib_scalar_route().

test_that("a binomial filter carries its size into the fast context", {
  set.seed(3)
  dd <- data.frame(id = factor(rep(1:4, each = 30)), t = rep(1:30, 4),
                   x = stats::rnorm(120))
  dd$y <- stats::rbinom(120, size = 10, prob = 0.4)
  spec <- statmod_spec(y ~ x + gas(p = 1, q = 1, by = id, time = t),
                       binomial_distrib(size = 10), dd)
  cb <- structural_callbacks(spec, list(mu = 0.4), "mu")
  expect_identical(cb$fast$family, "BinomialDistrib")
  expect_identical(names(cb$fast$theta), c("mu", "size"))
  expect_identical(cb$fast$theta$size, 10)
  expect_length(cb$fast$theta$mu, 120L)
})
