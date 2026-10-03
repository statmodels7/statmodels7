skip_on_cran()

test_that("a compartment names its random effects' sd without the entry prefix", {
  # two random developed parameters give the term two penalties, whose
  # hyperparameters are keyed term::entry; inside each compartment the
  # header already names the coefficient, so the row reads `effect sd`
  set.seed(3)
  id <- factor(rep(1:20, each = 15))
  t <- rep(seq(0, 10, length.out = 15), 20)
  psi_i <- 5 + rnorm(20, sd = 0.8)
  b0 <- 2 + rnorm(20, sd = 0.5)
  set.seed(6)
  g_i <- -1.2 + rnorm(20, sd = 0.3)
  y <- b0[id] + 0.8 * t + g_i[id] * pmax(t - psi_i[id], 0) + rnorm(300, sd = 0.4)
  d <- data.frame(id = id, t = t, y = y)
  f <- statmod(y ~ random(~ 1 | id) +
                 seg(t, psi ~ random(~ 1 | id), gamma1 ~ random(~ 1 | id)),
               distrib = distributions7::gaussian1_distrib(), data = d)
  out <- capture.output(print(summary(f)))
  # the names of an ordinary random() block, the one written in the equation
  # included (0.188.1)
  expect_length(grep("sigma [reml]", out, fixed = TRUE), 3L)
  expect_length(grep("effect sd", out, fixed = TRUE), 0L)
  expect_length(grep("::random", out, fixed = TRUE), 0L)
  # the keys hyper() reports are unchanged
  expect_true(all(grepl("::", hyper(f)$term[grepl("seg", hyper(f)$term)])))
})
