# The Step and the Tolerance of the Outer Hessian's Stencil

The step at which
[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md)
differences, the largest relative disagreement between its two readings,
and the tolerance of the probes' inner fit.

## Usage

``` r
hess_stencil_step()

hess_stencil_tol()

hess_stencil_inner_tol()
```

## Value

A single number.

## Details

`hess_stencil_step()` is `1e-3`, on the free scale of the
hyperparameters. The default of
[`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html),
about 6e-6, is intended for a function evaluated to machine precision.
Here the gradient is computed by refitting a mode, and the
reproducibility of that refit bounds the step from below.

`hess_stencil_tol()` is `1e-3`. It is compared with the largest absolute
difference between the matrices at `h` and `3 * h`, divided by the
largest absolute entry of the first. A relative error of `1e-3` in the
curvature corresponds to a relative error of about `5e-4` in a standard
error.

`hess_stencil_inner_tol()` is `1e-10`, the largest tolerance that the
probes use when the inner optimizer is
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md).
At the tolerance of a fit a probe can stop at its starting point, and
the difference then misses the movement of the mode. The resulting bias
does not grow as \\1/h\\, so the comparison of the two steps does not
detect it.

## See also

[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md)
