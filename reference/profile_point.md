# One Refit of a Profile

Fits the model with the labelled coefficients held at `val`. Without an
outer search it is the restricted fit of
[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md)
with several coefficients held. With one, the search runs over the
coordinates that the criterion estimates and, when `hyper = "refit"`,
over the hyperparameters, with the held coefficients excluded from it.

## Usage

``` r
profile_point(ctx, val, start = NULL, maxit = Inf)
```

## Arguments

- ctx:

  The context of the call.

- val:

  The values of the held coefficients, in the order of `ctx$which`.

- start:

  `NULL` to start at the coefficients of the fit, or the list returned
  by a previous call, whose `par` and `hyper` are the start.

## Value

A list with `par`, `hyper`, `objective`, `loglik`, `criterion` and
`converged`.

## Details

A break-point of a continuous
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
term is held as in
[`pin_held_breakpoints()`](https://statmodels7.github.io/statmodels7/reference/pin_held_breakpoints.md),
with one refit for each break-point that the held position can stand
for, the one with the smallest objective being kept.
