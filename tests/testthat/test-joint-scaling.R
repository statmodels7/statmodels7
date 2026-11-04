# The compiled derivatives of a joint filter's scaled driving vector
# (src/joint_scaling.cpp) against the R implementation they replaced
# (helper-joint-scaling-r.R), at every order and both scalings, on a
# two-parameter and a three-parameter family; with a control that moves
# theta by one part in a million.

jsc_case <- function(form, dist, dat, theta, order = 5L) {
  spec <- statmod_spec(form, dist, dat)
  tm <- spec@terms$mu$gas_joint
  y <- spec@response
  fam <- list(
    gl = distributions7::distrib_gradient(dist, y, theta, scale = "link"),
    H = distributions7::distrib_hessian(dist, y, theta, scale = "link"),
    D3 = distributions7::distrib_deriv3(dist, y, theta, scale = "link"),
    D4 = distributions7::distrib_deriv4(dist, y, theta, scale = "link"),
    D5 = distributions7::distrib_deriv5(dist, y, theta, scale = "link"),
    D6 = distributions7::distrib_deriv6(dist, y, theta, scale = "link"))
  list(spec = spec, tm = tm, theta = theta, fam = fam)
}

jsc_run <- function(cs, fun, theta = cs$theta, order = 5L) {
  f <- cs$fam
  fun(cs$spec, cs$tm, theta, f$gl, f$H, f$D3, f$D4, f$D5, order, D6 = f$D6)
}

jsc_rel <- function(a, b) {
  nm <- c("u", paste0("U", 1:5))
  vapply(nm, function(k) {
    if (is.null(b[[k]])) return(NA_real_)
    max(abs(a[[k]] - b[[k]])) / max(abs(b[[k]]))
  }, numeric(1))
}

jsc_data <- function(n = 60L, seed = 4L) {
  set.seed(seed)
  data.frame(t = seq_len(n), x = rnorm(n), y = rnorm(n, 1, 0.8),
             yp = rgamma(n, shape = 3, rate = 2),
             yt = 1 + 0.7 * rt(n, df = 6))
}

jsc_check <- function(form, dist, dat, theta) {
  for (d in c(1, 0.5)) {
    txt <- gsub("scaling = [0-9.]+", paste0("scaling = ", d), deparse1(form))
    cs <- jsc_case(stats::as.formula(txt), dist, dat, theta)
    expect_identical(cs$tm@terms[[1L]]@scaling, d)
    new <- jsc_run(cs, joint_drive_scaled)
    ref <- jsc_run(cs, joint_drive_scaled_r)
    rel <- jsc_rel(new, ref)
    expect_true(all(rel < 1e-12), info = paste(d, format(rel, digits = 3),
                                               collapse = " "))
    # every lower order alone, which the callers also request
    for (o in 1:4) {
      lo <- jsc_run(cs, joint_drive_scaled, order = o)
      expect_identical(lo[paste0("U", seq_len(o))], new[paste0("U", seq_len(o))])
      expect_true(all(vapply(lo[paste0("U", (o + 1L):5)], is.null, logical(1))))
    }
    # control: the reference at theta moved by one part in a million (every
    # parameter, the gaussian information not depending on mu)
    th2 <- lapply(theta, function(v) v * (1 + 1e-6))
    ctl <- jsc_rel(new, jsc_run(cs, joint_drive_scaled_r, theta = th2))
    expect_true(max(ctl, na.rm = TRUE) > 1e-9, info = format(ctl, digits = 3))
  }
}

test_that("the compiled scaling derivatives match the R route, gaussian1", {
  dat <- jsc_data()
  n <- nrow(dat)
  theta <- list(mu = 1 + 0.3 * sin(seq_len(n)),
                sigma = exp(-0.2 + 0.25 * cos(seq_len(n) / 3)))
  jsc_check(y ~ x + gas(p = 1, q = 1, time = t, scaling = 1) |
              sigma ~ gas(p = 1, q = 1, time = t, scaling = 1),
            distributions7::gaussian1_distrib(), dat, theta)
})

test_that("the compiled scaling derivatives match the R route, gamma1", {
  dat <- jsc_data()
  n <- nrow(dat)
  theta <- list(mu = exp(0.4 + 0.3 * sin(seq_len(n))),
                phi = exp(-0.5 + 0.2 * cos(seq_len(n) / 4)))
  jsc_check(yp ~ x + gas(p = 1, q = 1, time = t, scaling = 1) |
              phi ~ gas(p = 1, q = 1, time = t, scaling = 1),
            distributions7::gamma1_distrib(), dat, theta)
})

test_that("the compiled scaling derivatives match the R route, student_t1", {
  dat <- jsc_data()
  n <- nrow(dat)
  theta <- list(mu = 1 + 0.2 * sin(seq_len(n)),
                sigma = exp(-0.4 + 0.2 * cos(seq_len(n) / 3)),
                nu = 4 + sin(seq_len(n) / 5)^2)
  dist <- distributions7::student_t1_distrib()
  # two filtered parameters among three
  jsc_check(yt ~ x + gas(p = 1, q = 1, time = t, scaling = 1) |
              sigma ~ gas(p = 1, q = 1, time = t, scaling = 1) | nu ~ 1,
            dist, dat, theta)
  # all three filtered
  jsc_check(yt ~ x + gas(p = 1, q = 1, time = t, scaling = 1) |
              sigma ~ gas(p = 1, q = 1, time = t, scaling = 1) |
              nu ~ gas(p = 1, q = 1, time = t, scaling = 1),
            dist, dat, theta)
})

test_that("a non-positive information is rejected by the compiled route", {
  dat <- jsc_data()
  n <- nrow(dat)
  cs <- jsc_case(y ~ x + gas(p = 1, q = 1, time = t, scaling = 1) |
                   sigma ~ gas(p = 1, q = 1, time = t, scaling = 1),
                 distributions7::gaussian1_distrib(), dat,
                 list(mu = rep(1, n), sigma = rep(1, n)))
  I0 <- array(0, c(n, 2L, 2L))
  I0[, 1L, 1L] <- 1
  I0[, 2L, 2L] <- 1
  I0[5L, 2L, 2L] <- -1
  expect_false(jsc_s0_cpp(I0, n, 2L, 1)$ok)
  I0[5L, 2L, 2L] <- 1
  expect_true(jsc_s0_cpp(I0, n, 2L, 1)$ok)
})
