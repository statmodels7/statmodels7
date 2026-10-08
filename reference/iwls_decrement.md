# The Newton Decrement of a Scoring Step

`iwls_decrement()` is \\\tfrac{1}{2} g^\top (H+S)^{-1} g\\ over the free
coordinates, read from the step
[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
takes: the decrease in the penalized log-likelihood that a full scoring
step predicts. `iwls_decrement_limit()` is the value under which
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)'s
built-in rule reads the mode as located.

## Usage

``` r
iwls_decrement(pc, g, method, damp, frozen)

iwls_decrement_limit()
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

`iwls_decrement()`: a number, `NA` where the solve fails.
`iwls_decrement_limit()`: a number.

## Details

It is in log-likelihood units, so it does not depend on the scale of the
response, where the score per observation does. The limit is
\\10^{-6}\\, three orders under
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
which is the reading the outer search uses to decide whether a criterion
is usable.
