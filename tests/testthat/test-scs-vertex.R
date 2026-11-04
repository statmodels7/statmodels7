# A shape-constrained term at the vertex of its shape: the data go against
# the shape, every constrained coordinate runs to its bound and the block
# vanishes. The certificate reports it, the degrees of freedom count the free
# coordinates of the shape, and a smoothing parameter with no effect is not
# determined.

test_that("a term at the vertex is reported and counts its free coordinates", {
  set.seed(2)
  d <- data.frame(x = runif(300), g = gl(2, 150))
  d$yd <- -d$x + rnorm(300, sd = 0.1)
  d$yc <- -(d$x - 0.5)^2 + rnorm(300, sd = 0.05)
  d$yg <- ifelse(d$g == 1, d$x, -d$x) + rnorm(300, sd = 0.1)
  G <- distributions7::gaussian1_distrib()
  printed <- function(s) paste(utils::capture.output(print(s)), collapse = "
")
  for (cs in list(list(f = yd ~ modelterms7::scs(x, modelterms7::shape_increasing()), df = 2),
                  list(f = yc ~ modelterms7::scs(x, modelterms7::shape_convex()), df = 3))) {
    fit <- suppressWarnings(statmod(cs$f, distrib = G, data = d))
    ct <- statmod_certificate(fit)
    expect_identical(ct$state, "boundary")
    expect_length(ct$vertex, 1L)
    expect_match(paste(ct$reason, collapse = " "), "vertex of its shape")
    s <- suppressWarnings(summary(fit))
    expect_equal(s@df, cs$df, tolerance = 1e-6)
    expect_true(any(grepl("vertex of its shape", s@notes)))
    expect_match(printed(s), "lambda [not determined]", fixed = TRUE)
  }
  # a level at the vertex beside one that is not leaves the shared smoothing
  # parameter determined
  fit <- suppressWarnings(statmod(yg ~ g + modelterms7::scs(x, modelterms7::shape_increasing(), by = g),
                                  distrib = G, data = d))
  ct <- statmod_certificate(fit)
  expect_length(ct$vertex, 1L)
  expect_false(ct$vertex[[1L]]$complete)
  expect_false(grepl("not determined", printed(suppressWarnings(summary(fit))), fixed = TRUE))
  # and a fit inside the cone has no vertex
  d$y <- d$x + rnorm(300, sd = 0.1)
  fit <- statmod(y ~ modelterms7::scs(x, modelterms7::shape_increasing()), distrib = G, data = d)
  expect_length(statmod_certificate(fit)$vertex, 0L)
})
