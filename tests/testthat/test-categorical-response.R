# A binary outcome recorded as a factor.
#
# It is how these data are ordinarily written and glm() has always taken one.
# Before coerce_response() the response reached stats::dbinom() untouched and
# the run died on "Non-numeric argument to mathematical function", naming
# neither the variable nor the cause.

set.seed(1)
n <- 400
dd <- data.frame(
  x = rnorm(n),
  g = factor(sample(letters[1:3], n, TRUE)),
  stringsAsFactors = FALSE
)
dd$y01   <- rbinom(n, 1, plogis(-0.4 + 0.9 * dd$x))
dd$y_fac <- factor(ifelse(dd$y01 == 1, "si", "no"), levels = c("no", "si"))
dd$y_chr <- as.character(dd$y_fac)
dd$y_log <- dd$y01 == 1

test_that("factor, character, logical and 0/1 give the same fit", {
  rif <- coef(glm(y_fac ~ x + g, family = binomial, data = dd))
  for (v in c("y01", "y_fac", "y_chr", "y_log")) {
    f <- stats::reformulate(c("x", "g"), v)
    b <- coef(statmod(f, distributions7::bernoulli_distrib(), data = dd))$mu
    expect_equal(unname(b), unname(rif), tolerance = 1e-6, info = v)
  }
})

test_that("the first level is the failure, as in glm", {
  # the negative control: reversing the levels reverses the sign, and it has
  # to reverse by the same amount glm reverses it by
  dd$rev <- factor(dd$y_fac, levels = c("si", "no"))
  dritto  <- coef(statmod(y_fac ~ x, distributions7::bernoulli_distrib(), dd))$mu
  rovescio <- coef(statmod(rev ~ x, distributions7::bernoulli_distrib(), dd))$mu
  expect_equal(unname(dritto[2]), -unname(rovescio[2]), tolerance = 1e-5)
  expect_equal(unname(rovescio[2]),
               unname(coef(glm(rev ~ x, binomial, dd))[2]), tolerance = 1e-6)
})

test_that("the conversion reaches every family, where glm refuses", {
  # glm(gaussian) on a factor does not refuse cleanly: it signals
  # "NA/NaN/Inf in 'y'" after three warnings from Ops.factor. Here the
  # two-level factor has one numeric reading and the fit is the linear
  # probability model, which is what lm() on as.numeric() gives.
  b <- coef(statmod(y_fac ~ x, distributions7::gaussian1_distrib(), dd))$mu
  rif <- coef(lm(I(as.numeric(y_fac == "si")) ~ x, data = dd))
  expect_equal(unname(b), unname(rif), tolerance = 1e-6)
})

test_that("anything but two levels is refused, by name", {
  dd$y3 <- factor(sample(c("a", "b", "c"), n, TRUE))
  expect_error(
    statmod(y3 ~ x, distributions7::bernoulli_distrib(), dd),
    "factor with 3 levels")
  dd$y1 <- factor(rep("solo", n))
  expect_error(
    statmod(y1 ~ x, distributions7::bernoulli_distrib(), dd),
    "1 level")
})

test_that("a numeric response is returned untouched", {
  # the conversion must not put an attribute on a response that had none, or
  # the many places that compare the stored response against a fresh column
  # would stop matching
  spec <- statmod_spec(y01 ~ x, distributions7::bernoulli_distrib(), dd)
  expect_identical(spec@response, dd$y01)
})

test_that("other rows are coded with the mapping the fit used", {
  # A character vector does not carry its levels through a subset, so a fold
  # holding only one of the two values would be read the other way round and
  # the deviance would come back wrong without a word. The levels travel on
  # the fitted specification instead.
  fit <- statmod(y_chr ~ x, distributions7::bernoulli_distrib(), dd)
  solo_no <- dd[dd$y_chr == "no", ][1:20, ]
  spec <- statmod_respec(fit@spec, solo_no)
  expect_true(all(spec@response == 0))

  solo_si <- dd[dd$y_chr == "si", ][1:20, ]
  expect_true(all(statmod_respec(fit@spec, solo_si)@response == 1))

  # e la verosimiglianza a nuovi dati e' quella dei valori adattati li'
  expect_equal(as.numeric(logLik(fit)),
               as.numeric(logLik(fit, newdata = dd)), tolerance = 1e-8)
})

test_that("a value the fit never saw is refused rather than coded", {
  fit <- statmod(y_chr ~ x, distributions7::bernoulli_distrib(), dd)
  altro <- dd[1:10, ]
  altro$y_chr <- "forse"
  expect_error(statmod_respec(fit@spec, altro), "not among the levels")
})
