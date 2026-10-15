# What the Outer Criterion Can Resolve at This Fit

The size of a difference in the criterion that carries no information,
computed from the fit at the point.

## Usage

``` r
criterion_resolution(st, spec, design, method, criterion_at)
```

## Arguments

- st:

  An evaluation's stored state, carrying `ok`, `par`, `split`, `score`,
  `cf`, `hy`, `ctx` and `value`.

- spec:

  The specification.

- design:

  The design.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- criterion_at:

  A function of `(cf, hy, par, ctx)` returning what the search is
  running, a prediction-error criterion or a marginal one. It is passed
  in and not chosen here, because the marginal criterion of a fit whose
  search runs
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  is a quantity that the search does not evaluate.

## Value

A single positive number, or `NA_real_` where the pieces are not
available. Where the mode error exceeds
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
the `NA` carries the attribute `mode_error`, which holds that error.

## Details

The criterion is read at the penalized mode, and the inner fit stops
short of that mode by whatever its rule allows. Writing \\g\\ for the
score at which the inner fit stopped and \\K = H + S\\ for the penalized
information, the mode is out by \\\delta\beta = K^{-1} g\\, and the
criterion read at \\\hat\beta - \delta\beta\\ instead of at
\\\hat\beta\\ differs by an amount that equals the difference between
two evaluations from different warm starts. One assembly of the
criterion at given coefficients gives it, with no refit.

The criterion carries \\-\tfrac{1}{2}\log\lvert K(\beta)\rvert\\, which
is not stationary in \\\beta\\, so a mode error enters it at first order
and the quadratic form \\\tfrac{1}{2} g' K^{-1} g\\ alone does not bound
the change. The criterion is therefore evaluated at the displaced
coefficients.

The reading estimates the spread of the criterion across warm starts. It
is usually below that spread, because the spread is a range over several
paths and the reading is one deviation. A resolution smaller than the
true one leaves the search where it was, and a larger one stops a search
before it has converged.

The reading is `NA` where \\g' K^{-1} g / 2\\ exceeds
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
because a displacement computed from an inner fit that is far from its
mode is not a resolution.

## See also

[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md),
[`optimizers7::armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.html)
