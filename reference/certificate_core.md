# The Certificate's Readings Before the Edge Check

Everything
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
computes except the check at the edge of a chart, which that function
applies to the result.

## Usage

``` r
certificate_core(fit, tol, flat, edge)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- tol:

  The largest rise, in the criterion's own units, that a certified point
  may still have available to it. ⚠️ Until 0.127.0 this was a threshold
  on the outer gradient. The default has not moved and its meaning has:
  a caller who set one should read
  [`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md).

- flat:

  The largest curvature \\\lvert A\_{jj}\rvert\\ of the outer criterion
  in a hyperparameter's free value at which that hyperparameter is
  reported as sitting at a boundary; for a coefficient the criterion
  estimates, the curvature times \\\max(1, \gamma_j^2)\\. It is also the
  largest eigenvalue, in absolute value, of a flat direction of the
  equilibrated curvature. It decides the label alone and never the
  verdict; see the details.

- edge:

  The free value beyond which a hyperparameter is reported as sitting at
  a boundary where no curvature can be read, and only there.

## Value

The list
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
returns.
