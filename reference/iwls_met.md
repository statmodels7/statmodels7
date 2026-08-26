# Has the Step's Stopping Rule Been Met?

Applies whichever stopping rule the step was configured with: the
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

Both routes live here, in one place, although the loop asks the question
at two. That is what keeps the state a criterion is shown identical in
both.

Two properties of that state matter to a caller writing a rule.
`state$gradient` is already the score per observation, so `crit_grad(t)`
and `tol = t` are the same rule. `state$objective` is the penalized
log-likelihood unaveraged, that being the scale the penalty is added on.
