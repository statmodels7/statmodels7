# What the Marginal Criterion Can Resolve at This Fit

The size of the difference in the criterion that carries no information,
measured at the point. Nothing about it is assumed.

## Usage

``` r
criterion_resolution(st, spec, design, method, criterion_at)
```

## Arguments

- st:

  An evaluation's stored state, carrying `par`, `split`, `score`, `cf`,
  `hy`, `ctx` and `value`.

- spec:

  The specification.

- design:

  The design.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- criterion_at:

  A function of `(cf, hy, par, ctx)` returning what **this** search is
  running, a prediction-error criterion or a marginal one. It is passed
  in and never chosen here: reading the marginal criterion of a fit
  whose search is
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  answers for a quantity the search never sees, and the number that came
  back stopped two such fits short of their own optimum.

## Value

A single positive number, or `NA_real_` where the pieces are not
available.

## Details

The criterion is read at the penalized mode, and the inner fit stops
short of that mode by whatever its rule allows. Writing \\g\\ for the
score it stopped at and \\K = H + S\\ for the penalized information, the
mode is out by \\\delta\beta = K^{-1} g\\, and the criterion read at
\\\hat\beta - \delta\beta\\ instead of at \\\hat\beta\\ differs by an
amount that is exactly what an evaluation from a different warm start
would differ by. One assembly of the criterion at given coefficients
answers it, with no refit.

**The quadratic form alone is not enough**, and the reason is structural
before it is a matter of accuracy: the criterion carries
\\-\tfrac{1}{2}\log\lvert K(\beta)\rvert\\, which is not stationary in
\\\beta\\, so a mode error enters it at **first** order. Measured,
\\\tfrac{1}{2} g' K^{-1} g\\ is right to one per cent at an inner
tolerance of `1e-4`, where the second-order term dominates, and
undershoots by 50 to 1000 times at `1e-6` and below, where the first
order does.

Measured against the spread of the criterion at one hyperparameter
reached from six different warm starts, over four shapes and five inner
tolerances, the displaced reading tracks it across six orders of
magnitude, from `8.6e-2` down to `2.2e-7`, at a ratio between 0.05 and
0.99, and separates shapes that a formula cannot: at an inner tolerance
of `1e-6` it reads `1.7e-4` on a random intercept over 500 levels
against `7.9e-8` on a gaussian smooth.

It reads low by two to three times, the spread being a range over six
paths where this is one deviation, and that is the side to err on. A
resolution smaller than the truth leaves the search where it was; one
larger stops a healthy search short.

## See also

[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md),
[`optimizers7::armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.html)
