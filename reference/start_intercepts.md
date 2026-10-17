# Start at the Intercept-Only Fit

This is the default strategy: each equation's intercept starts at the
maximum likelihood estimate of the model that has no covariates. Where
the first parameter is a location or a mean on the identity link, the
slopes of its parametric columns start at a least-squares fit, and a
term with a starting value of its own starts at the value that the term
declares. The remaining coefficients start at zero.

## Usage

``` r
start_intercepts()
```

## Value

A `StartIntercepts` object, inheriting from
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md).

## Details

The model with no covariates is the model in which every slope is zero,
so its maximum likelihood estimate is a natural starting point for the
intercepts of the model with covariates. It requires one small fit, on
an intercept-only design.

Three further rules apply:

- A location or mean on the identity link whose equation has an
  intercept and at least one other parametric column starts its
  parametric coefficients at a weighted least-squares fit net of the
  offset, and the intercepts of the other parameters are fitted to the
  residuals. This applies when every term of that equation is a
  parametric, smooth, random or penalized term.

- A term whose
  [`modelterms7::term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.html)
  method returns values that are not all zeros, such as a break-point or
  a nonlinear term, starts at those values.

- All other coefficients start at zero, which is where the penalty of a
  penalized block is smallest.

Two adjustments apply to the intercepts:

- The mean of an offset is subtracted from the intercept before the
  value is written, because the intercept-only fit is computed without
  the offset and the predictor would otherwise start displaced by it.

- The value is written to a parametric intercept and to nothing else. A
  term such as
  [`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
  also names the intercept of each of its own parameters `(Intercept)`,
  but those coefficients lie on the chart of that parameter and not on
  the scale of the predictor.

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
uses this strategy when `start` is `NULL`. The constructor exists so
that the default can be named, compared with other strategies and passed
on.

## See also

[`start_origin()`](https://statmodels7.github.io/statmodels7/reference/start_origin.md),
[`start_random()`](https://statmodels7.github.io/statmodels7/reference/start_random.md),
[`start_search()`](https://statmodels7.github.io/statmodels7/reference/start_search.md)

## Examples

``` r
start_intercepts()
#> <start> intercept-only fit

set.seed(1)
dd <- data.frame(x = runif(60, 0, 10))
dd$y <- 500 + 20 * dd$x + rnorm(60, sd = 5)
spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)

# The intercept starts at the response's own scale, and the slope at its
# least-squares value.
start_at(start_intercepts(), spec, statmod_design(spec), NULL)
#> $mu
#> [1] 500.88837  19.94518
#> 
#> $sigma
#> [1] 1.4661
#> 
```
