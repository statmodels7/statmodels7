# A SMOOTH SHOWS ITS LINEAR PART WHATEVER CLASS ITS PENALTY IS.
#
# smooth_linear_cols() reads a penalty's matrix off the property "P", which a
# quadratic penalty has and an ADDITIVE one does not -- its components live
# in "mats", one per smoothing parameter. Read through the "P" question alone
# a penalty of that class fell into the branch written for a penalty carrying
# no matrix at all, a lasso or a heavy-tailed prior, which shrinks every
# coordinate it indexes; the block then reported its smoothing parameters and
# nothing else, losing the row a reader of a smooth most wants.
#
# It is the defect statmodels7 0.118.0 repaired for the SINGULAR question,
# one class further along.

test_that("a penalty carrying components reports its free leading columns", {
  set.seed(3)
  n <- 200
  d <- data.frame(x = sort(stats::runif(n)), z = stats::runif(n))
  d$y <- sin(2 * pi * d$x) + stats::rnorm(n, sd = 0.2)

  free_of <- function(tm) {
    bt <- modelterms7::term_build(tm, d)
    sum(smooth_linear_cols(bt, ncol(modelterms7::term_matrix(bt))))
  }

  # the family whose penalty is a sum of components
  ad <- modelterms7::s(x, basis7::adaptive_smooth(k = 20, m = 4))
  expect_identical(
    class(modelterms7::term_penalty(modelterms7::term_build(ad, d)))[[1L]],
    "penalties7::AdditivePenalty")
  expect_identical(free_of(ad), 1L)
  # at order 3 the null space is the lines and the parabolas, so TWO columns
  # are free: the count follows the construction rather than being fixed
  expect_identical(
    free_of(modelterms7::s(x, basis7::adaptive_smooth(k = 20, m = 4, diff = 3))),
    2L)
  expect_identical(
    free_of(modelterms7::s(x, basis7::adaptive_smooth(k = 20, m = 4,
                                                      null_space = "drop"))),
    0L)

  # THE READING IS WIDENED AND NOT REPLACED: the shapes that answered before
  # answer the same, which is what says the repair is confined
  expect_identical(free_of(modelterms7::s(x, basis7::bspline_smooth(k = 15))), 1L)
  expect_identical(free_of(modelterms7::s(x, basis7::pspline_smooth(k = 15))), 1L)
  # a tensor product is centered and has no free leading column under either
  # penalty class, which is why the defect was unreachable before the
  # adaptive family existed
  expect_identical(
    free_of(modelterms7::te(x, z, smooths = basis7::bspline_smooth(k = 4))), 0L)
  expect_identical(
    free_of(modelterms7::te(x, z, smooths = basis7::bspline_smooth(k = 4),
                            anisotropic = FALSE)), 0L)
})

test_that("the summary of an adaptive smooth carries its linear row", {
  set.seed(7)
  n <- 200
  d <- data.frame(x = sort(stats::runif(n)))
  d$y <- sin(2 * pi * d$x) + stats::rnorm(n, sd = 0.25)
  fit <- statmod(y ~ s(x, basis7::adaptive_smooth(k = 15, m = 3)),
                 distributions7::gaussian1_distrib(), d)
  sm <- summary(fit)
  blk <- Filter(function(b) grepl("adaptive_smooth", b$term, fixed = TRUE),
                sm@tables$mu)
  expect_length(blk, 1L)
  nms <- blk[[1L]]$table[[1L]]
  # the three smoothing parameters AND the unpenalized column, which carries
  # the term's LABEL as a prefix -- "s(x).lin" -- where the printed form
  # strips it. The label does not carry the construction; the key does.
  expect_true(all(c("lambda1", "lambda2", "lambda3") %in% nms))
  expect_true(any(endsWith(nms, ".lin")))
  expect_identical(sum(endsWith(nms, ".lin")), 1L)
})
