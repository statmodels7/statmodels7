# Has the Step's Stopping Rule Been Met?

Applies the stopping rule with which the step was configured: the
built-in one, the score per observation against `tol`, or the caller's
optimizers7 criterion read on the same state.

## Usage

``` r
iwls_met(method, state)
```

## Arguments

- method:

  An
  [`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

- state:

  The iteration state, as
  [`optimizers7::crit_met`](https://statmodels7.github.io/optimizers7/reference/crit_met.html)
  documents it.

## Value

A single logical.

## Details

Both routes are implemented here, although the loop evaluates the rule
at two places. This keeps the state that a criterion receives identical
at both places.

In the state that a rule receives, `state$gradient` is already the score
per observation, so `crit_grad(t)` and `tol = t` apply the same
threshold to the score; the built-in rule also requires a small Newton
decrement, which a criterion does not receive. `state$f_new` is the
negative penalized log-likelihood, not divided by the sample size, which
is the scale on which the penalty is added.
