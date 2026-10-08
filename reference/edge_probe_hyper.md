# Probe the Hyperparameters at the Edge of Their Chart From Inside

Moves each hyperparameter whose free value has run past `edge` back
towards its starting value, with the other coordinates held, and returns
the free vector of the best probe where it improves the criterion by
more than `tol`.

## Usage

``` r
edge_probe_hyper(par, start, nh, value, fn, edge = 8, tol = mode_error_limit())
```

## Arguments

- par:

  The point the search reported: the free values of the hyperparameters
  first, then any coefficients the criterion estimates.

- start:

  The starting point of the search, of the same length.

- nh:

  The number of hyperparameters at the head of `par`.

- value:

  The value of the search's objective at `par`, which it minimizes.

- fn:

  The search's objective.

- edge:

  The free value past which a coordinate is at the edge, the default of
  [`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md).

- tol:

  The gain, in the criterion's units, that counts:
  [`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md).

## Value

The free vector of the best probe, or `NULL` where no coordinate is at
the edge or no probe gains more than `tol`.

## Details

A hyperparameter lives on a chart such as \\\sigma = e^\eta\\ or \\\rho
= \tanh(\eta/2)\\, and the slope of the criterion in \\\eta\\ is its
slope in the bounded quantity times \\d\sigma/d\eta\\ or
\\d\rho/d\eta\\, which tend to zero at the edge. A search can therefore
stop at the edge with a vanishing gradient where the criterion still
rises towards the inside. A finite step reads that rise where a
derivative cannot: on the free scale the step to \\\eta = \pm 8\\ is
still at the edge (a standard deviation of 3.4e-4), so each coordinate
is tried at three quarters, one half, one quarter and the whole of the
way back to its starting value. A coordinate that really sits at an
edge, such as the smoothing parameter of a smooth of noise, loses
criterion at every probe and is left where it is; what that costs is
four evaluations of the criterion, one inner fit each.
