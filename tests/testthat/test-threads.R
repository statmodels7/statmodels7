# The thread policy: the result does not depend on the count, bit for bit.
# The assembly kernel decomposes over the elements of its output and each
# element is accumulated in full by one thread in the sequential order, so
# identical() is the right comparison and not a tolerance. Tests ask for 2
# threads, which is what CRAN's check machines have.

test_that("wcrossprod matches crossprod bit for bit, dense and above threshold", {
  set.seed(1)
  n <- 4000; p <- 30                       # n * p * p above the 2e6 gate
  A <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("a", 1:p)))
  B <- matrix(rnorm(n * p), n, p)
  w <- runif(n, 0.5, 2)
  expect_identical(wcrossprod(A, w, B, threads = 2L), crossprod(A * w, B))
  expect_identical(wcrossprod(A, w, B, threads = 2L),
                   wcrossprod(A, w, B, threads = 1L))
})

test_that("xtv and wxsq match their expressions bit for bit", {
  set.seed(9)
  n <- 8000; p <- 40                       # n * p above the 2e5 gate
  X <- matrix(rnorm(n * p), n, p)
  v <- rnorm(n)
  w <- runif(n, 0.5, 2)
  expect_identical(xtv(X, v, threads = 2L), as.numeric(crossprod(X, v)))
  expect_identical(wxsq(X, w, threads = 2L), as.numeric(crossprod(w, X^2)))
  expect_identical(xtv(X, v, threads = 2L), xtv(X, v, threads = 1L))
  expect_identical(wxsq(X, w, threads = 2L), wxsq(X, w, threads = 1L))
  expect_identical(xtx(X, threads = 2L), crossprod(X))
  expect_identical(xtx(X, threads = 2L), xtx(X, threads = 1L))
})

test_that("a sparse design keeps its Matrix route", {
  skip_if_not_installed("Matrix")
  set.seed(2)
  X <- Matrix::rsparsematrix(3000, 40, density = 0.01)
  w <- runif(3000, 0.5, 2)
  out <- wcrossprod(X, w, X, threads = 2L)
  ref <- crossprod(X * w, X)
  expect_identical(class(out), class(ref))
  expect_identical(out, ref)
})

test_that("a whole fit does not depend on threads, bit for bit", {
  set.seed(3)
  n <- 4000; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  b <- c(2, -1.5, 1, numeric(p - 3))
  d <- data.frame(y = as.numeric(X %*% b) + rnorm(n))
  d$X <- X
  f1 <- statmod(y ~ lasso(X), distributions7::gaussian1_distrib(), d)
  f2 <- statmod(y ~ lasso(X), distributions7::gaussian1_distrib(), d,
                threads = numericals7::n_threads(2))
  expect_identical(f1@coefficients, f2@coefficients)
  expect_identical(as.numeric(logLik(f1)), as.numeric(logLik(f2)))
  expect_identical(vcov(f1), vcov(f2))
})

test_that("a fit through the family's parallel kernels is bit-identical too", {
  set.seed(4)
  n <- 6000                                # above gamma1's kernel threshold
  dd <- data.frame(x = runif(n))
  mu <- exp(1 + sin(2 * pi * dd$x))
  dd$y <- rgamma(n, shape = 1 / 0.3, scale = mu * 0.3)
  f1 <- statmod(y ~ s(x, k = 10), distributions7::gamma1_distrib(), dd)
  f2 <- statmod(y ~ s(x, k = 10), distributions7::gamma1_distrib(), dd,
                threads = numericals7::n_threads(2))
  expect_identical(f1@coefficients, f2@coefficients)
  expect_identical(as.numeric(logLik(f1)), as.numeric(logLik(f2)))
})

test_that("a cross-validated fit does not depend on workers, bit for bit", {
  skip_on_cran()   # spawns worker processes; the sequential twin is covered
  set.seed(5)
  n <- 300; p <- 12
  X <- matrix(rnorm(n * p), n, p)
  b <- c(1.5, -1, 0.75, numeric(p - 3))
  d <- data.frame(y = as.numeric(X %*% b) + rnorm(n))
  d$X <- X
  set.seed(42)
  f1 <- statmod(y ~ lasso(X, n_lambda = 6), distributions7::gaussian1_distrib(),
                d, sparse_criterion = cv(nfolds = 3))
  set.seed(42)
  f2 <- statmod(y ~ lasso(X, n_lambda = 6), distributions7::gaussian1_distrib(),
                d, sparse_criterion = cv(nfolds = 3),
                threads = numericals7::n_threads(1, workers = 2))
  expect_identical(f1@coefficients, f2@coefficients)
  expect_identical(as.numeric(logLik(f1)), as.numeric(logLik(f2)))
  expect_identical(f1@hyper, f2@hyper)
})

test_that("a product-grid path does not depend on workers, bit for bit", {
  skip_on_cran()   # spawns worker processes; the sequential twin is covered
  # the combinations of the product are the independent unit: each restarts
  # its warm chain from the sweep's starting coefficients, so the same
  # bodies run whatever the count -- scad's grid is n_lambda x n_a runs
  set.seed(9)
  n <- 300; p <- 15
  X <- matrix(rnorm(n * p), n, p)
  b <- c(1.2, -0.9, 0.6, numeric(p - 3))
  d <- data.frame(y = as.numeric(X %*% b) + rnorm(n))
  d$X <- X
  f1 <- statmod(y ~ scad(X, n_lambda = 6, n_a = 3),
                distributions7::gaussian1_distrib(), d,
                sparse_criterion = bic())
  f2 <- statmod(y ~ scad(X, n_lambda = 6, n_a = 3),
                distributions7::gaussian1_distrib(), d,
                sparse_criterion = bic(),
                threads = numericals7::n_threads(1, workers = 2))
  expect_identical(f1@coefficients, f2@coefficients)
  expect_identical(f1@hyper, f2@hyper)
  expect_identical(f1@criterion, f2@criterion)
})

test_that("threads must be the n_threads() object, and the setting is restored", {
  d <- data.frame(x = rnorm(60))
  d$y <- 1 + d$x + rnorm(60)
  expect_error(statmod(y ~ x, distributions7::gaussian1_distrib(), d,
                       threads = 4),
               "n_threads")
  old <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)
  invisible(statmod(y ~ x, distributions7::gaussian1_distrib(), d,
                    threads = numericals7::n_threads(2)))
  expect_identical(Sys.getenv("RCPP_PARALLEL_NUM_THREADS",
                              unset = NA_character_), old)
})
