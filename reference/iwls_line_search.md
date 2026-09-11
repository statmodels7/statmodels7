# One Scoring Step and Its Line Search

Solves for the scoring increment from the given pieces and halves the
step until Armijo's condition holds or the budget of halvings is spent.

## Usage

``` r
iwls_line_search(obj, beta, value, g, pc, method, damp, frozen, guard = FALSE)
```

## Arguments

- obj:

  The objective.

- beta:

  The current coefficients.

- value:

  The objective there.

- g:

  The gradient there.

- pc:

  The pieces
  [`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
  reads.

- method:

  The
  [`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

- damp:

  The Levenberg damping.

- frozen:

  The held positions.

- guard:

  Whether an error at a trial point is read as a rejection.

## Value

A list of `ok`, `cand`, `vnew`, `step_used`, `sol` and `errors`, the
number of trial points whose objective raised.

## Details

It is the body
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
ran inline, moved into a function so that the step can be tried twice at
one iterate, first on the observed information and then on the expected
one, where
[`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md)
settled the method with the fallback. With `guard = FALSE` the
operations are the inline ones, in the same order. With `guard = TRUE` a
trial point whose objective raises an error is read as rejected rather
than stopping the fit, and counted.

## See also

[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md),
its caller.
