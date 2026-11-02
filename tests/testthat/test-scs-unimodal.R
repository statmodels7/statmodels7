# The search over the versions of a term: scs() with shape_unimodal() and no
# mode is fitted at the indices around the pilot's peak and the best criterion
# is kept; with a mode it is fitted once.

uni_data <- function() {
  set.seed(71)
  x <- runif(200)
  data.frame(x = x, y = exp(-(x - 0.4)^2 / (2 * 0.15^2)) + rnorm(200, sd = 0.2))
}

test_that("a unimodal scs() without a mode is fitted by a search", {
  d <- uni_data()
  fit <- statmod(y ~ modelterms7::scs(x, modelterms7::shape_unimodal()),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  va <- fit@history$variants
  expect_true(is.data.frame(va))
  expect_identical(sum(va$chosen), 1L)
  expect_true(nrow(va) >= 3L && nrow(va) <= 5L)
  # the chosen version has the best score among those fitted
  expect_identical(which(va$chosen), which.min(va$score))
  mu <- predict(fit, what = "mu", newdata = data.frame(x = seq(0.02, 0.98, length.out = 301)))
  df <- diff(mu); s <- sign(df[abs(df) > 1e-10])
  expect_lte(sum(diff(s) != 0), 1L)
  expect_match(paste(capture.output(print(fit)), collapse = "\n"), "chosen among")
  # the search does not leave its state behind
  expect_null(statmod_variants$active)
  expect_null(statmod_variants$override)
})

test_that("a unimodal scs() with a mode is fitted from its two starts", {
  d <- uni_data()
  fit <- statmod(y ~ modelterms7::scs(x, modelterms7::shape_unimodal(mode = 0.4)),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  expect_identical(fit@history$variants$version, c("start: sides", "start: spline", "start: sides, lambda 0.01", "start: spline, lambda 0.01"))
  mu <- predict(fit, what = "mu", newdata = data.frame(x = seq(0.02, 0.98, length.out = 301)))
  expect_lt(abs(seq(0.02, 0.98, length.out = 301)[which.max(mu)] - 0.4), 0.15)
})

test_that("variant_score() makes every criterion a minimum", {
  d <- uni_data()
  # REML is maximized, so its score is the criterion with its sign turned
  f1 <- statmod(y ~ x, distrib = distributions7::gaussian1_distrib(), data = d)
  expect_identical(variant_score(f1), -f1@criterion)
  # without an outer criterion the score is the penalized objective
  f2 <- statmod(y ~ x, distrib = distributions7::gaussian1_distrib(), data = d,
                outer_criterion = NULL)
  expect_identical(variant_score(f2), f2@objective)
})

test_that("a version's starting smoothing parameter reaches the outer search", {
  d <- uni_data()
  v <- modelterms7::scs(x, modelterms7::shape_unimodal(index = 4L))
  v@spec$lambda_start <- 0.01
  v@spec$start_rule <- "sides"
  spec <- statmod_spec(y ~ modelterms7::scs(x, modelterms7::shape_unimodal(index = 4L)),
                       distrib = distributions7::gaussian1_distrib(), data = d)
  spec@terms$mu[[2L]] <- modelterms7::term_build(v, d)
  h <- statmod_hyper_start(spec)
  expect_equal(unlist(h$mu)[[1L]], 0.01)
})
