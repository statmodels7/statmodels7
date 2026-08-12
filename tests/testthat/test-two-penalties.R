# A term may carry more than one penalty, over different parameters of its
# own. Nothing shipped in modelterms7 that statmod() can fit does, so the
# term here is a subclass of a penalized block that splits itself in two and
# declares one penalty per half. What is under test is the layer: that the
# two penalties get keys, hyperparameters, an outer estimate, a degrees of
# freedom count and a summary row each.

TwoPen <- S7::new_class("TwoPen", parent = modelterms7::PenalizedTerm)

# S7::method()<- takes the generic by value and registers on it in place
term_penalties <- modelterms7::term_penalties
S7::method(term_penalties, TwoPen) <- function(term, ...) {
  k <- modelterms7::term_npar(term)
  h <- k %/% 2L
  list(list(name = "lo", index = seq_len(h),
            penalty = penalties7::ridge_penalty(n_coef = h)),
       list(name = "hi", index = h + seq_len(k - h),
            penalty = penalties7::ridge_penalty(n_coef = k - h)))
}

twopen <- function(x, label = "twopen") {
  base <- modelterms7::ridge(x, label = label)
  TwoPen(label = base@label, input = base@input,
         input_expr = base@input_expr,
         # the halves carry the penalties, so there is none over the whole
         factory = function(k) NULL,
         X = NULL, coef_names = character(0),
         blueprint = list(), penalty = NULL)
}

set.seed(4)
n <- 200
dd <- data.frame(z1 = stats::rnorm(n), z2 = stats::rnorm(n),
                 w1 = stats::rnorm(n), w2 = stats::rnorm(n))
# the first half carries signal and the second is noise, so the two halves
# have different amounts of shrinkage to ask for
dd$y <- 1 + 1.5 * dd$z1 - 1.2 * dd$z2 + stats::rnorm(n, sd = 0.5)

fml <- y ~ twopen(~ z1 + z2 + w1 + w2)
# a term is named by the call that produced it, not by its label
tn <- "twopen(~z1 + z2 + w1 + w2)"

test_that("a term's two penalties get a key each", {
  spec <- statmod_spec(fml, distributions7::gaussian1_distrib(), dd)
  keys <- statmod_penalty_keys(spec)
  expect_length(keys, 2L)
  expect_identical(vapply(keys, function(u) u$key, character(1)),
                   paste0(tn, c("::lo", "::hi")))
  expect_identical(vapply(keys, function(u) u$term, character(1)),
                   c(tn, tn))
  # the positions within the term, and then within the stacked coefficients
  units <- statmod_penalized(spec, statmod_design(spec))
  expect_identical(units[[1L]]$within, 1:2)
  expect_identical(units[[2L]]$within, 3:4)
  expect_false(any(units[[1L]]$index %in% units[[2L]]$index))

  # a term carrying one penalty keys under its own name, as it always did
  s1 <- statmod_spec(y ~ ridge(~ z1 + z2),
                     distributions7::gaussian1_distrib(), dd)
  expect_identical(vapply(statmod_penalty_keys(s1), function(u) u$key,
                          character(1)), "ridge(~z1 + z2)")
})

test_that("the two are estimated apart and counted together", {
  # The budget is raised above its default rather than the flag being
  # weakened. The alternation stops on a relative change, whose attainable
  # floor is platform arithmetic: at the default macOS reported
  # converged = FALSE where the other four platforms reported TRUE, with the
  # same hyperparameters and the same degrees of freedom to the digit, so what
  # differed was how close to the floor the run had to get and not the answer.
  fit <- statmod(fml, distributions7::gaussian1_distrib(), dd,
                 inner_method = iwls(maxit = 500L),
                 outer_method = reml())
  expect_true(fit@converged,
              info = sprintf("criterion %.6g, objective %.10g",
                             fit@criterion, fit@objective))
  h <- fit@hyper$mu
  expect_true(all(paste0(tn, c("::lo", "::hi")) %in% names(h)))
  # the half that carries no signal is shrunk harder, which is the whole
  # point of giving the two halves a hyperparameter each
  expect_lt(h[[paste0(tn, "::hi")]][["sigma"]],
            h[[paste0(tn, "::lo")]][["sigma"]])

  # one row per term, and the count is between the null space and the rank
  e <- fit@edf
  row <- e[e$parameter == "mu" & e$term == tn, ]
  expect_identical(nrow(row), 1L)
  expect_false(is.na(row$edf))
  expect_gt(row$edf, 0)
  expect_lt(row$edf, row$coefficients + 1e-8)

  # and it is the count modelterms7 gives from the same pieces, which is
  # what says the hyperparameters were routed to the penalties they belong
  # to rather than one being used for the whole block
  des <- statmod_design(fit@spec)
  tm <- fit@spec@terms$mu[[tn]]
  cols <- des$mu$blocks[[tn]]
  H <- statmod_information_at(fit@spec, fit@coefficients, des, TRUE, "bartlett")
  ref <- modelterms7::edf(
    tm, coef = fit@coefficients$mu[cols], hessian = H[cols, cols, drop = FALSE],
    theta = list(lo = as.list(h[[paste0(tn, "::lo")]]),
                 hi = as.list(h[[paste0(tn, "::hi")]])))
  # LOOSELY, because the two are not the same quantity: modelterms7 reads
  # the term's own block and statmod_edf reads that term's share of the
  # WHOLE model's smoother matrix, which differs by the coupling between
  # this block and the rest. A hyperparameter routed to the wrong penalty
  # moves the count by far more than that.
  expect_equal(row$edf, ref, tolerance = 1e-3)

  # and exactly, against the definition: the trace of this term's diagonal
  # block of (H + S)^-1 H over the coefficients of every equation
  S <- statmod_penalty_at(fit@spec, fit@coefficients, fit@hyper, des,
                          "hessian")
  S[!is.finite(S)] <- 0
  expect_equal(row$edf, sum(diag(solve(H + S, H))[cols]), tolerance = 1e-10)
})

test_that("the summary shows a hyperparameter per penalty, named for it", {
  fit <- statmod(fml, distributions7::gaussian1_distrib(), dd,
                 outer_method = reml())
  sm <- summary(fit)
  blk <- Filter(function(b) identical(b$term, tn), sm@tables$mu)
  expect_length(blk, 1L)
  expect_identical(blk[[1L]]$kind, "penalized")
  nms <- blk[[1L]]$table[[1L]]
  # two hyperparameters in one block are not the same number and cannot
  # appear under the same name
  expect_true(all(c("lo.sigma", "hi.sigma") %in% nms))
  expect_output(print(sm), "twopen")
})

test_that("a partially penalized term is not read as parametric", {
  set.seed(5)
  dx <- data.frame(x = sort(stats::runif(150, 0, 10)))
  dx$y <- 1 + 0.5 * dx$x + 2 * pmax(dx$x - 6, 0) + stats::rnorm(150, sd = 0.3)
  # seg penalizes its changes and nothing else, so term_penalty() is NULL
  # while term_penalties() is not: reading the first would file the term
  # under the unpenalized ones
  built <- modelterms7::term_build(modelterms7::seg(x, penalty = "lasso"), dx)
  expect_null(modelterms7::term_penalty(built))
  expect_identical(term_block_kind(built), "selection")
  expect_identical(
    term_block_kind(modelterms7::term_build(modelterms7::seg(x), dx)),
    "parametric")
  expect_identical(
    term_block_kind(modelterms7::term_build(
      modelterms7::seg(x, penalty = "ridge"), dx)), "penalized")
})
