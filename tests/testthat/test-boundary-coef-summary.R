test_that("a coefficient estimated on the criterion at a boundary has no standard error", {
  set.seed(1)
  xb <- data.frame(x = runif(n = 400))
  sim <- rstatmod(y ~ x | sigma ~ 1, distrib = gaussian1_distrib(), data = xb,
                  par = list(mu = c(1, 2), sigma = log(0.5)))
  fit <- statmod(y ~ x | sigma ~ 1, distrib = student_t1_distrib(),
                 data = sim$data)
  cert <- statmod_certificate(fit)
  expect_identical(cert$state, "boundary")
  expect_identical(cert$boundary_coef, "nu:(Intercept)")
  s <- summary(fit)
  tab <- s@tables$nu[[1L]]$table
  expect_true(is.finite(tab$estimate))
  expect_true(all(is.na(tab[, c("se", "statistic", "p_value", "lower",
                                "upper")])))
  expect_true(any(grepl("nu:(Intercept) was estimated on the outer criterion",
                        s@notes, fixed = TRUE)))
  # the coefficients that are not at the boundary keep their standard errors
  expect_true(all(is.finite(s@tables$mu[[1L]]$table$se)))
  expect_true(all(is.finite(s@tables$sigma[[1L]]$table$se)))
})
