# The Newton Decrement of a Scoring Step

`iwls_decrement()` is \\\tfrac{1}{2} g^\top (H+S)^{-1} g\\ over the free
coordinates, read from the step that
[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
takes: the decrease in the penalized log-likelihood that a full scoring
step predicts.

## Usage

``` r
iwls_decrement(pc, g, method, damp, frozen)
```

## Arguments

- pc:

  Scoring pieces at the current coefficients.

- g:

  The gradient of the objective there.

- method:

  An
  [`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

- damp:

  The current Levenberg damping.

- frozen:

  Positions held out of the step.

## Value

A number, `NA` where the solve fails.

## Details

It is in log-likelihood units, so it does not depend on the scale of the
response, where the score per observation does.
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)'s
built-in rule reads the mode as located when it is at most the
`decrement_limit` of
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
whose default \\10^{-6}\\ is three orders under
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
the reading used by the outer search to determine whether a criterion is
usable.
