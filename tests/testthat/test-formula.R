# Splitting the multi-parameter formula.
#
# The tests carry the PARSE TREE and not only the result. R's `~` binds
# looser than `|` and associates to the left, so the shape the recovery
# depends on is a fact about the parser; a change in it would otherwise be
# discovered as a model fitted to the wrong equations rather than as an error.

test_that("the parse tree is the one the recovery assumes", {
  f <- y ~ a | p2 ~ b | p3 ~ c
  # the whole formula's right-hand side is the LAST piece alone, which is why
  # splitting it on `|` returns one piece and drops the rest
  expect_identical(f[[3L]], quote(c))
  expect_identical(f[[1L]], quote(`~`))
  # and the left operand is itself a `~` call, twice over
  expect_identical(f[[2L]][[1L]], quote(`~`))
  expect_identical(f[[2L]][[2L]][[1L]], quote(`~`))
  # the naive reading, recorded so that its failure is visible
  naive <- function(e) {
    if (is.call(e) && identical(e[[1L]], quote(`|`))) {
      c(naive(e[[2L]]), list(e[[3L]]))
    } else {
      list(e)
    }
  }
  expect_length(naive(f[[3L]]), 1L)
})

test_that("every equation is recovered, in the family's order", {
  r <- statmod_equations(y ~ x1 + x2 | sigma ~ z | nu ~ 1,
                         c("mu", "sigma", "nu"))
  expect_identical(r$response, quote(y))
  expect_named(r$equations, c("mu", "sigma", "nu"))
  expect_identical(r$given, c("mu", "sigma", "nu"))
  expect_identical(r$equations$mu[[2L]], quote(x1 + x2))
  expect_identical(r$equations$sigma[[2L]], quote(z))
  expect_identical(r$equations$nu[[2L]], quote(1))
})

test_that("the equations may be given out of order", {
  r <- statmod_equations(y ~ x | nu ~ w | sigma ~ z, c("mu", "sigma", "nu"))
  expect_named(r$equations, c("mu", "sigma", "nu"))
  expect_identical(r$equations$sigma[[2L]], quote(z))
  expect_identical(r$equations$nu[[2L]], quote(w))
  # `given` records what the formula supplied, in the order it supplied it
  expect_identical(r$given, c("mu", "nu", "sigma"))
})

test_that("a parameter with no equation gets an intercept", {
  r <- statmod_equations(y ~ x, c("mu", "sigma", "nu"))
  expect_named(r$equations, c("mu", "sigma", "nu"))
  expect_identical(r$equations$sigma[[2L]], quote(1))
  expect_identical(r$equations$nu[[2L]], quote(1))
  expect_identical(r$given, "mu")
})

test_that("a bar inside a call is untouched", {
  # this is what makes random(1 | id) and gas(by = ~ random(1 | id)) work: the
  # walk descends only through `~` and the top-level `|`
  r <- statmod_equations(y ~ x1 + random(1 | id) | sigma ~ random(1 | site),
                         c("mu", "sigma"))
  expect_identical(r$equations$mu[[2L]], quote(x1 + random(1 | id)))
  expect_identical(r$equations$sigma[[2L]], quote(random(1 | site)))

  g <- statmod_equations(
    y ~ gas(p = 1, q = 1, by = ~ random(1 | id)) | theta ~ random(1 | id),
    c("mu", "theta"))
  expect_identical(g$equations$mu[[2L]],
                   quote(gas(p = 1, q = 1, by = ~ random(1 | id))))
  expect_identical(g$equations$theta[[2L]], quote(random(1 | id)))
})

test_that("the response survives whatever it is", {
  r <- statmod_equations(cens(y, lwr = 0) ~ x | sigma ~ 1, c("mu", "sigma"))
  expect_identical(r$response, quote(cens(y, lwr = 0)))
  r2 <- statmod_equations(log(y) ~ x, c("mu", "sigma"))
  expect_identical(r2$response, quote(log(y)))
})

test_that("the equations carry the formula's own environment", {
  # a term's symbols must resolve where the user wrote them
  local({
    marker <- 42
    f <- y ~ x | sigma ~ z
    r <- statmod_equations(f, c("mu", "sigma"))
    for (eq in r$equations) {
      expect_identical(get("marker", envir = environment(eq)), 42)
    }
  })
})

test_that("what the formula gets wrong is named", {
  expect_error(statmod_equations(y ~ x | wrong ~ z, c("mu", "sigma")),
               "not a parameter")
  expect_error(statmod_equations(y ~ x | wrong ~ z, c("mu", "sigma")),
               "mu, sigma")
  expect_error(statmod_equations(y ~ x | mu ~ z, c("mu", "sigma")),
               "more than one equation")
  expect_error(statmod_equations(~ x, c("mu")), "left-hand side")
  expect_error(statmod_equations("y ~ x", c("mu")), "must be a formula")
})

test_that("modelterms7 shadows the search path", {
  env <- statmodels7:::terms_first(baseenv())
  # a term call resolves to ours whatever is attached
  expect_identical(get("s", envir = env), modelterms7::s)
  expect_identical(get("ridge", envir = env), modelterms7::ridge)
  expect_identical(get("seg", envir = env), modelterms7::seg)
  # and everything behind it stays visible
  expect_identical(get("sum", envir = env), base::sum)
})
