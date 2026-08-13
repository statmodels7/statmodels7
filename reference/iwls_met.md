# Has the Step's Stopping Rule Been Met?

The built-in rule, the score per observation against `tol`, or the
caller's optimizers7 criterion read on the same state.

## Usage

``` r
iwls_met(method, state)
```

## Arguments

- method:

  An
  [`Iwls`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

- state:

  The iteration state, as
  [`optimizers7::crit_met`](https://statmodels7.github.io/optimizers7/reference/crit_met.html)
  documents it.

## Value

A single logical.

## Details

The two routes are here rather than at the two places the loop asks, so
that what a rule is shown is written once. `state$gradient` is already
the score PER OBSERVATION, which is what makes `crit_grad(t)` and
`tol = t` the same rule; the objective in the state is the penalized
log-likelihood unaveraged, the scale the penalty is added on.
