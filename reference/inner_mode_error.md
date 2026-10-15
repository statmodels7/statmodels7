# How Far Above Its Mode the Inner Fit Stopped

\\\tfrac12 g'K^{-1}g\\, the decrease the penalized likelihood's own
Newton correction predicts at the point the inner fit returned: how much
log-likelihood is still on the table there.

## Usage

``` r
inner_mode_error(ctx, spec, design, coef, hyper, score, expected = FALSE)
```

## Arguments

- ctx:

  The evaluation context, so the penalized factorization is the one the
  criterion will read rather than a second copy.

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients the inner fit returned.

- hyper:

  The hyperparameters.

- score:

  The inner objective's gradient at those coefficients.

- expected:

  Whether the penalized information is the expected one.

## Value

A single number, or `NA` where the penalized system could not be read
there, which is itself a reason to call the point unavailable.

## Details

It decides whether a point is available, which is a different question
from the one that the inner optimizer's flag answers. The flag reports
whether a stopping rule fired. Availability is whether the criterion, a
Laplace expansion at the mode, is valid at this point.

Availability is a matter of distance and has a natural scale,
log-likelihood units. The flag is a boolean about a threshold on a score
whose size depends on the model. The two can differ: an inner fit may
stop on the objective-stall guard with the score slightly above the
absolute tolerance and the penalized likelihood already at its mode to
many digits. Such a point is available by this reading and unavailable
by the flag.

It is used to add points and not to remove any: a run whose flag reports
convergence stays usable whatever this reads.

## See also

[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
