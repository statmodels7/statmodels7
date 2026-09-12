test_that("the Newton decrement is the rise a quadratic model still offers", {
  # It is defined as max_x (2 g'x - x'Ax) / 2, so the cheapest reference is
  # that maximisation done directly, which shares no arithmetic with the
  # solve.
  set.seed(7)
  for (p in 1:4) {
    B <- matrix(rnorm(p * p), p, p)
    A <- crossprod(B) + diag(p)
    g <- rnorm(p)
    dec <- joint_decrement(g, A)
    # the maximum of 2 g'x - x'Ax over x, located by optim rather than solved
    f <- function(x) -(2 * sum(g * x) - drop(crossprod(x, A %*% x)))
    o <- stats::optim(rep(0, p), f, method = "BFGS",
                      control = list(reltol = 1e-14, maxit = 2000))
    expect_equal(dec, -o$value / 2, tolerance = 1e-6)
  }
})


test_that("restricting a decrement to fewer coordinates cannot raise it", {
  # THIS IS WHAT KEEPS `edge` A LABEL AND NEVER A VERDICT. The certificate
  # takes the coordinates at a boundary out of the interior set before
  # reading the verdict, so the reading over what is left must be no larger
  # than the reading over everything -- otherwise removing a coordinate could
  # turn a certified point into an uncertified one and the threshold would be
  # deciding the state. It is an identity of the constrained maximum and not
  # a tolerance, so it is asserted over random matrices rather than on one.
  set.seed(11)
  for (i in 1:40) {
    p <- sample(2:6, 1)
    B <- matrix(rnorm(p * p), p, p)
    A <- crossprod(B) + diag(p)
    g <- rnorm(p)
    all <- joint_decrement(g, A)
    keep <- sample(seq_len(p), sample(1:(p - 1), 1))
    sub <- joint_decrement(g[keep], A[keep, keep, drop = FALSE])
    expect_lte(sub, all + 1e-9)
    # and the same statement one coordinate at a time, which is the reading
    # the boundary label uses
    expect_true(all(coord_decrement(g, A) <= all + 1e-9))
  }
})


test_that("a decrement is refused where the curvature is not definite", {
  # A point whose curvature is not negative definite in every direction under
  # test is not a maximum there, and no rise follows from a quadratic model
  # around it. NA rather than a number is the answer, and the certificate
  # then reports "unknown" with the reason.
  A <- matrix(c(2, 0, 0, -1), 2, 2)
  expect_true(is.na(joint_decrement(c(1, 1), A)))
  expect_true(is.na(joint_decrement(c(1, 1), matrix(c(1, 2, 2, 1), 2, 2))))
  expect_true(is.na(joint_decrement(c(1, NA), diag(2))))
  # an empty interior set is not a refusal: there is nothing left to rise in
  expect_identical(joint_decrement(numeric(0), matrix(0, 0, 0)), 0)

  # a coordinate whose own curvature is not positive is never called settled,
  # so the boundary label cannot excuse it from the verdict
  expect_identical(coord_decrement(c(1, 1), A)[[2L]], Inf)
  expect_equal(coord_decrement(c(3, 1), diag(c(2, 8)))[[1L]], 9 / 4)
})


test_that("resolution_summary keeps the largest usable reading", {
  # ⚠️ THE LARGEST AND NOT THE SMALLEST, which is what 0.127.0 changed.
  # criterion_resolution() already returns NA where the mode error exceeds
  # mode_error_limit(), so every reading reaching here comes from a mode the
  # layer would accept; what a resolution has to bound is how far the
  # criterion CAN move between two such modes, not how little it happened to
  # move at the luckiest evaluation. Measured on a pig1 smooth, the readings
  # taken during one search run 1.1e-10 to 1.2e-07 while the criterion's own
  # spread over six warm starts is 4.7e-07 -- the largest is within a factor
  # of four of it and the smallest 4300 times below.
  expect_identical(resolution_summary(c(3, 1, 2)), 3)
  expect_identical(resolution_summary(c(1e-10, 6.6e-08, 1.2e-07)), 1.2e-07)
  # non-finite and non-positive readings are not readings
  expect_identical(resolution_summary(c(1, NA, Inf, -2, 0, 4)), 4)
  expect_true(is.na(resolution_summary(numeric(0))))
  expect_true(is.na(resolution_summary(c(NA_real_, 0, -1))))
  # it can only grow as a search takes more readings, which is what makes a
  # run stop sooner rather than later as it learns
  set.seed(5)
  x <- abs(rnorm(20)) + 1e-6
  run <- vapply(seq_along(x), function(i) resolution_summary(x[1:i]), numeric(1))
  expect_false(is.unsorted(run))
})


test_that("the certificate's verdict is the decrement and its cost is stated", {
  skip_on_cran()
  set.seed(41)
  n <- 300
  d <- data.frame(x = runif(n))
  d$y <- sin(5 * d$x) + rnorm(n, 0, 0.3)
  fit <- statmod(y ~ s(x, bspline_smooth(k = 10)), gaussian1_distrib(), d,
                 outer_criterion = reml())
  ct <- statmod_certificate(fit)
  expect_identical(ct$state, "converged")
  expect_true(is.finite(ct$decrement))
  expect_lt(ct$decrement, 1e-2)
  # the gradient is still reported beside it, and it is a DIFFERENT number:
  # a verdict read off the decrement is not the same verdict read in other
  # units, which is the whole point of the change
  expect_true(is.finite(ct$gradient))
  expect_gt(ct$gradient, ct$decrement)
  expect_identical(ct$curvature, "analytic")

  # IT MUST BE ABLE TO REFUSE, or "converged" says nothing, and the refusal
  # now names the rise rather than the gradient
  strict <- statmod_certificate(fit, tol = ct$decrement / 10)
  expect_identical(strict$state, "not converged")
  expect_match(strict$reason[1], "would still rise")
})


test_that("the decrement certifies a healthy fit the gradient reading refuses", {
  # THE CASE THE CHANGE EXISTS FOR. The outer gradient of a criterion summed
  # over n observations and p penalized coefficients carries both, so one
  # absolute threshold cannot serve every shape: a random intercept over
  # n/10 groups beside a smooth reaches gradients of 0.02 to 0.4 at points
  # that are at their own optimum. Measured over 1350 fits against an
  # independently located optimum, `max|g| > 1e-2` flags 131 of 663 fits
  # within 1e-3 of theirs while the decrement at the same cut flags none.
  #
  # Here the premise is asserted rather than assumed: the fit is compared
  # against the same model searched hard, and only then is the pair of
  # readings compared. A fit that really were short would fail the premise
  # and skip rather than pass.
  skip_on_cran()
  set.seed(101)
  n <- 3000L
  m <- n %/% 10L
  x <- runif(n)
  id <- factor(rep(seq_len(m), length.out = n))
  b <- rnorm(m, 0, 0.5)
  sgn <- 2 * rbinom(n, 1, 0.5) - 1
  d <- data.frame(y = 1 + sin(2 * pi * x) + b[as.integer(id)] + sgn * rexp(n, 1),
                  x = x, id = id)
  f <- y ~ random(~ 1 | id) + s(x, bspline_smooth(k = 10))

  fit <- statmod(f, gaussian1_distrib(), d)
  ref <- statmod(f, gaussian1_distrib(), d,
                 outer_optimizer = optimizers7::newton(
                   criterion = optimizers7::crit_grad(1e-8),
                   line_search = optimizers7::armijo(max_step = 15),
                   maxit = 50L),
                 inner_optimizer = iwls(tol = 1e-10))
  gap <- ref@criterion - fit@criterion
  skip_if(!is.finite(gap) || gap > 1e-3,
          paste("the premise does not hold on this platform: the fit sits",
                format(gap, digits = 3), "below an independently located",
                "optimum, so it is not a healthy fit to certify"))

  ct <- statmod_certificate(fit)
  expect_identical(ct$state, "converged")
  expect_lt(ct$decrement, 1e-2)
  # and the reading it replaces would have refused this point
  expect_gt(ct$gradient, 1e-2)
})


test_that("a coordinate at an edge is named whatever its own curvature says", {
  # ⚠️ THE SECOND THING CI FOUND, and the deeper of the two. At a boundary the
  # curvature collapses along with the gradient, so the one-coordinate reading
  # g^2/(2A) is a ratio of two quantities going to zero and can come back
  # LARGE -- exactly where `edge` matters most. Reporting `boundary` on the
  # conjunction therefore left such a coordinate unnamed, and `boundary_key`
  # is what summary() reads to leave a standard error off it.
  #
  # What is REPORTED is now the value alone; what may be EXCLUDED from the
  # verdict keeps the conjunction, which is what stops `edge` deciding the
  # state. The degenerate curvature is built here rather than fitted for, the
  # model that produces one being weakly identified and landing differently on
  # every platform -- which is the whole reason this was not seen locally.
  skip_on_cran()
  set.seed(42)
  d <- data.frame(x = runif(300), y = rnorm(300))
  fit <- statmod(y ~ s(x, bspline_smooth(k = 10)), gaussian1_distrib(), d,
                 outer_criterion = reml())
  plain <- statmod_certificate(fit)
  skip_if(!length(plain$boundary), "this fit did not reach the chart's edge")
  expect_identical(plain$state, "boundary")

  # the same fit, with a curvature so nearly singular in that coordinate that
  # its own decrement is far past `tol`
  local_mocked_bindings(
    outer_curvature = function(...) list(A = matrix(1e-30, 1L, 1L),
                                         source = "analytic",
                                         why = character(0)))
  ct <- statmod_certificate(fit)
  # NAMED, exactly as before
  expect_identical(ct$boundary, plain$boundary)
  expect_identical(ct$boundary_key, plain$boundary_key)
  # and NOT excused: it stays under test, so the verdict is not "boundary"
  expect_gt(ct$decrement, 1e-2)
  expect_identical(ct$state, "not converged")
  expect_match(ct$reason[1], "would still rise")
})


test_that("with no curvature at all there is no verdict but still a boundary", {
  # ⚠️ THE REGRESSION CI FOUND AND THIS MACHINE HID. A verdict in criterion
  # units needs a curvature, and where neither the analytic route nor the
  # difference produces one there is none to give -- but `boundary_key` is
  # what summary() reads to leave a standard error off a coordinate pinned at
  # an edge, and returning early with none takes that away from a reader on
  # exactly the fits least able to spare it. The state stays "unknown" and
  # says why, so nothing is certified on a value alone.
  #
  # The branch is reached by mocking the curvature rather than by building
  # the model that reaches it -- a mixed covariance class inside a filter,
  # which costs a minute and lands differently on every platform.
  skip_on_cran()
  set.seed(42)
  d <- data.frame(x = runif(300), y = rnorm(300))
  fit <- statmod(y ~ s(x, bspline_smooth(k = 10)), gaussian1_distrib(), d,
                 outer_criterion = reml())
  plain <- statmod_certificate(fit)
  # the premise: this fit's one hyperparameter really is past the edge, or
  # the assertion below would hold for the wrong reason
  skip_if(!length(plain$boundary), "this fit did not reach the chart's edge")

  local_mocked_bindings(
    outer_curvature = function(...) list(A = NULL, source = NA_character_,
                                         why = "mocked: no curvature here"))
  ct <- statmod_certificate(fit)
  expect_identical(ct$state, "unknown")
  expect_identical(ct$boundary, plain$boundary)
  expect_identical(ct$boundary_key, plain$boundary_key)
  expect_true(is.finite(ct$gradient))
  expect_true(is.na(ct$decrement))
  expect_match(ct$reason, "no curvature here")
})


test_that("the curvature is differenced where the form has no analytic one", {
  # Two forms carry an exact outer GRADIENT and no analytic Hessian: a
  # criterion asked for on the expected information, and a separable penalty.
  # Both are certified, by differencing the exact gradient, and the route is
  # reported so a reader can tell which was taken.
  skip_on_cran()
  set.seed(41)
  n <- 300
  d <- data.frame(x = runif(n))
  d$y <- sin(5 * d$x) + rnorm(n, 0, 0.3)
  f <- y ~ s(x, bspline_smooth(k = 10))

  obs <- statmod(f, gaussian1_distrib(), d, outer_criterion = reml())
  exp_ <- statmod(f, gaussian1_distrib(), d, outer_criterion = reml("expected"))
  co <- statmod_certificate(obs)
  ce <- statmod_certificate(exp_)
  expect_identical(co$curvature, "analytic")
  expect_identical(ce$curvature, "differenced")
  expect_true(is.finite(ce$decrement))
  expect_identical(ce$state, "converged")
  # and the printed line says which route it rested on, so a reader of a
  # summary does not have to ask the object
  po <- capture.output(print(summary(obs)))
  pe <- capture.output(print(summary(exp_)))
  expect_false(any(grepl("curvature differenced", po, fixed = TRUE)))
  expect_true(any(grepl("curvature differenced", pe, fixed = TRUE)))

  # THE TWO ROUTES AGREE WHERE BOTH EXIST, which is what licenses the second
  # one. The analytic Hessian of the observed fit against one differenced
  # from its own exact gradient: the decrement read either way.
  spec <- obs@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  meth <- obs@methods$outer
  basis <- integrated_basis(spec, design, meth@kind)
  g <- statmod_marginal_grad(spec, design, obs@coefficients, obs@hyper, meth,
                             idx, basis)
  an <- outer_curvature(spec, design, obs@coefficients, obs@hyper, meth, idx,
                        basis, obs@methods$smooth)
  Hd <- statmod_hess_stencil(spec, design, obs@coefficients, obs@hyper, meth,
                             idx, basis, obs@methods$smooth)
  skip_if(is.null(Hd), "the stencil did not resolve a curvature here")
  expect_identical(an$source, "analytic")
  expect_equal(joint_decrement(g, an$A),
               joint_decrement(g, -(Hd + t(Hd)) / 2), tolerance = 1e-3)
})
