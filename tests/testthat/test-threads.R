# The thread policy: the KERNEL's result does not depend on the count, bit
# for bit -- it decomposes over the elements of its output and each element
# is accumulated in full by one thread in the sequential order, so the twin
# there is identical() and not a tolerance. The comparison against the BLAS
# expression the kernel replaces is a different claim: an optimized BLAS
# (OpenBLAS on the CI's Linux runners, Accelerate on macOS) blocks and
# vectorizes its accumulations, so the two implementations agree to the last
# bit only on the reference BLAS, and that comparison carries a tolerance
# chosen so a real defect (a split reduction, a dropped weight) still fails
# it by many orders. Asserting identical() there reddened four of the five
# CI platforms on 2026-08-18 while this machine's reference BLAS agreed.
# Tests ask for 2 threads, which is what CRAN's check machines have.

test_that("the kernels do not depend on the thread count, bit for bit", {
  set.seed(1)
  n <- 4000; p <- 30
  A <- matrix(rnorm(n * p), n, p)
  B <- matrix(rnorm(n * p), n, p)
  w <- runif(n, 0.5, 2)
  expect_identical(wcrossprod_cpp(A, w, B, 2L), wcrossprod_cpp(A, w, B, 1L))
  X <- matrix(rnorm(8000 * 40), 8000, 40)
  v <- rnorm(8000)
  w2 <- runif(8000, 0.5, 2)
  expect_identical(xtv_cpp(X, v, 2L), xtv_cpp(X, v, 1L))
  expect_identical(wxsq_cpp(X, w2, 2L), wxsq_cpp(X, w2, 1L))
  expect_identical(xtx_cpp(X, 2L), xtx_cpp(X, 1L))
})

test_that("the kernels match the expressions they replace", {
  set.seed(1)
  n <- 4000; p <- 30                       # n * p * p above the 2e6 gate
  A <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("a", 1:p)))
  B <- matrix(rnorm(n * p), n, p)
  w <- runif(n, 0.5, 2)
  expect_equal(wcrossprod(A, w, B, threads = 2L), crossprod(A * w, B),
               tolerance = 1e-12)
  set.seed(9)
  n <- 8000; p <- 40                       # n * p above the 2e5 gate
  X <- matrix(rnorm(n * p), n, p)
  v <- rnorm(n)
  w <- runif(n, 0.5, 2)
  expect_equal(xtv(X, v, threads = 2L), as.numeric(crossprod(X, v)),
               tolerance = 1e-12)
  expect_equal(wxsq(X, w, threads = 2L), as.numeric(crossprod(w, X^2)),
               tolerance = 1e-12)
  expect_equal(xtx(X, threads = 2L), crossprod(X), tolerance = 1e-12)
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

test_that("a whole fit does not depend on threads", {
  # cross-route: threads = 1 fits through the BLAS expressions and threads = 2
  # through the kernels, so the comparison carries a tolerance -- a race or a
  # dropped weight moves a coefficient by orders more than the last bit
  set.seed(3)
  n <- 4000; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  b <- c(2, -1.5, 1, numeric(p - 3))
  d <- data.frame(y = as.numeric(X %*% b) + rnorm(n))
  d$X <- X
  f1 <- statmod(y ~ lasso(X), distributions7::gaussian1_distrib(), d)
  f2 <- statmod(y ~ lasso(X), distributions7::gaussian1_distrib(), d,
                threads = numericals7::n_threads(2))
  expect_equal(f1@coefficients, f2@coefficients, tolerance = 1e-8)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-8)
  expect_equal(vcov(f1), vcov(f2), tolerance = 1e-8)
})

test_that("a fit through the family's parallel kernels agrees too", {
  set.seed(4)
  n <- 6000                                # above gamma1's kernel threshold
  dd <- data.frame(x = runif(n))
  mu <- exp(1 + sin(2 * pi * dd$x))
  dd$y <- rgamma(n, shape = 1 / 0.3, scale = mu * 0.3)
  f1 <- statmod(y ~ s(x, k = 10), distributions7::gamma1_distrib(), dd)
  f2 <- statmod(y ~ s(x, k = 10), distributions7::gamma1_distrib(), dd,
                threads = numericals7::n_threads(2))
  expect_equal(f1@coefficients, f2@coefficients, tolerance = 1e-8)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-8)
})

test_that("a cross-validated fit does not depend on workers, bit for bit", {
  skip_on_cran()   # spawns worker processes; the sequential twin is covered
  # workers run the same sequential bodies whatever the count, so this twin
  # IS identical()
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

# A fold's fit is built by statmod_spec(), which makes a FRESH
# specification, so the thread count does not travel with it as it does
# through statmod_respec(), which starts from the one it is given. Until
# 2026-08-21 a fold therefore fell back to the class default of 1 and a
# cross-validated path was single-threaded whatever statmod(threads =)
# asked for -- measured on a lasso over a gamma response, bic() gained
# 1.85x from eight threads where cv() gained 1.06x, and 1.91x after.
#
# The assertion is on the ANSWER and on the count the fold receives, not on
# a timing: what the defect broke was the promise, and a fold that fits
# faster is not observable from here.
test_that("a cross-validated fit uses the thread count it was given", {
  set.seed(11)
  n <- 400; p <- 8
  X <- matrix(rnorm(n * p), n, p)
  d <- data.frame(y = as.numeric(X %*% c(1.2, -0.8, numeric(p - 2))) + rnorm(n))
  d$X <- X
  f1 <- statmod(y ~ lasso(X, n_lambda = 5), distributions7::gaussian1_distrib(),
                d, sparse_criterion = cv(nfolds = 3))
  f2 <- statmod(y ~ lasso(X, n_lambda = 5), distributions7::gaussian1_distrib(),
                d, sparse_criterion = cv(nfolds = 3),
                threads = numericals7::n_threads(2))
  expect_identical(f1@coefficients, f2@coefficients)
  expect_identical(as.numeric(logLik(f1)), as.numeric(logLik(f2)))

  # what the fold is handed, read where it is set: the count in this
  # process, and 1 where the folds go to worker processes of their own
  spec <- f2@spec
  expect_identical(spec@threads, 2L)
  fold_threads <- function(sp) if (sp@workers > 1L) 1L else sp@threads
  expect_identical(fold_threads(spec), 2L)
  spec@workers <- 4L
  expect_identical(fold_threads(spec), 1L)
})

# The triangular factor of the augmented system. augmented_solve() reads
# only R, the pivot and the rank off the decomposition -- Q is never
# accumulated, applied or returned -- so a kernel that produces R alone does
# the whole job. Step j of the Householder reduction applies one reflector
# to each trailing column and those updates are independent, so column k is
# written in full by one thread in the order the sequential loop writes it:
# the factor is bit-identical at any count BY CONSTRUCTION, which a block
# decomposition over rows could not offer, its answer depending on the
# partition.
test_that("the triangular factor does not depend on the thread count", {
  set.seed(21)
  n <- 6000; p <- 40
  A <- matrix(rnorm(n * p), n, p)
  # the column scales a penalized augmented matrix actually has
  A <- A * rep(10^seq(0, 8, length.out = p), each = n)
  ref <- qr_factor_cpp(A, 1L)
  for (k in c(2L, 3L, 5L)) expect_identical(qr_factor_cpp(A, k), ref)

  # R'R is A'A, which is the only property the solve rests on: the solve is
  # invariant to whatever orthogonal factor sits on the left, so the sign
  # convention of the rows does not enter
  ata <- crossprod(A)
  expect_equal(crossprod(ref), ata, tolerance = 1e-12)

  # and the increment it produces is R's own, to the last bits
  u <- rnorm(p)
  qa <- qr(A)
  dq <- numeric(p)
  dq[qa$pivot] <- backsolve(qr.R(qa), forwardsolve(t(qr.R(qa)), u[qa$pivot]))
  dk <- backsolve(ref, forwardsolve(t(ref), u))
  expect_equal(dk, dq, tolerance = 1e-10)
})

test_that("a dense penalized fit does not depend on the thread count", {
  set.seed(22)
  n <- 4000
  x <- runif(n, -3, 3)
  d <- data.frame(x = x, y = 1 + sin(x) + rnorm(n, 0, 0.4))
  f1 <- statmod(y ~ s(x, k = 30), distributions7::gaussian1_distrib(), d)
  f2 <- statmod(y ~ s(x, k = 30), distributions7::gaussian1_distrib(), d,
                threads = numericals7::n_threads(2))
  expect_equal(f1@coefficients, f2@coefficients, tolerance = 1e-8)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-8)
})
