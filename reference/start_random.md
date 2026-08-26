# Start From a Random Draw

Draws every coefficient from a generator of the caller's choosing and,
by default, adds the draw to the intercept-only fit. Built for running a
model from several starts and comparing where they land, which is how a
multimodal likelihood is explored.

## Usage

``` r
start_random(fn = stats::rnorm, ..., center = TRUE)
```

## Arguments

- fn:

  A generator called as `fn(n, ...)` and returning `n` values. Defaults
  to [`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html). Any
  function of that shape serves.

- ...:

  Further arguments to `fn`, captured at construction and stored on the
  object: `sd = 0.1`, or `min` and `max` for
  [`stats::runif()`](https://rdrr.io/r/stats/Uniform.html).

- center:

  A single logical. `TRUE`, the default, adds the draw to the
  intercept-only start; `FALSE` uses the draw alone.

## Value

A `StartRandom` object, inheriting from
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md),
with properties `fn`, `args` and `center`.

## Why the draw is centered by default

A coefficient drawn from a standard normal is a sensible perturbation
and a hopeless absolute value: the intercept of a location equation is
on the scale of the response, which may be a thousand. Adding the draw
to
[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md)'s
answer keeps that scale and randomizes the direction, leaving the
several starts something to differ in.

`center = FALSE` uses the draw alone, which is right when the generator
is already scaled to the problem.

## The random stream

The draw comes from the session's generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) governs the result
and a fit begun this way repeats only alongside its seed. Nothing is
seeded internally.

## See also

[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md)
for what the draw is added to,
[`optimizers7::multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.html)
for running several starts inside an optimizer.

## Examples

``` r
start_random()
#> <start> random around the intercept-only fit
start_random(stats::runif, min = -2, max = 2)
#> <start> random around the intercept-only fit
start_random(sd = 0.1)@args
#> $sd
#> [1] 0.1
#> 

set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 100 + 2 * dd$x + rnorm(40, sd = 0.3)
spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
design <- statmod_design(spec)

# Centered, the intercept stays on the response's scale.
set.seed(2)
start_at(start_random(sd = 0.5), spec, design, NULL)$mu
#> [1] 100.60044079   0.09242459

# Uncentered, it is the raw draw, which is nowhere near 100.
set.seed(2)
start_at(start_random(sd = 0.5, center = FALSE), spec, design, NULL)$mu
#> [1] -0.44845727  0.09242459
```
