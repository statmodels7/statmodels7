# What a Hyperparameter's Own Uncertainty Adds to a Variance

\\J V\_\theta J'\\, the movement of the penalized mode under the
estimated hyperparameters, as a matrix to be added to \\V_b\\.

## Usage

``` r
hyper_correction(spec, design, coef, hyper, method, Vb, keep, nz = 0L)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

- hyper:

  The hyperparameters.

- method:

  The outer method that estimated them, or `NULL`.

- Vb:

  The bayesian variance the correction is to be added to, over the
  coordinates `keep` names followed by any structural tail.

- keep:

  A logical vector over the stacked coefficients, saying which ones `Vb`
  spans.

- nz:

  How many structural parameters `Vb` carries beyond them.

## Value

A list with `C`, the correction or `NULL`, `n_hyper`, how many
hyperparameters a differentiable criterion estimated, and `complete`,
whether every one of them contributed.

## Details

The two matrices
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
reports by default are conditional on the hyperparameters: both are read
at the value the outer search stopped at, as though it had been known.
It was estimated from the same data, and what that costs is \$\$V' =
V_b + J V\_\theta J', \qquad J = -(H + S)^{-1} \frac{\partial^2
\rho}{\partial\beta \partial\theta},\$\$ the delta method applied to the
map from the hyperparameter to the mode (Wood, Pya and Safken, 2016). It
is the same quantity
[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md)
contracts against the information to obtain a count of parameters; here
the matrix itself is what is wanted.

\\V_b\\ is PASSED IN rather than recomputed, and that is what keeps the
two halves of the sum describing one model. The caller has already
settled which information \\H\\ is, which coordinates are held and which
are aliased; a correction built on a second inverse, regularized
differently, would not be the movement of the mode whose variance it is
added to.

## Where there is nothing to add, and where it cannot be read

The correction is exactly zero where no hyperparameter was estimated by
a differentiable criterion. A kinked penalty's is the argument of a
minimum over a grid, which
[`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
skips, and the map from it to the mode turns a corner whenever a
coefficient joins or leaves the active set, so there is no derivative to
propagate. That is a property of the model and not a failure, and
`n_hyper` is zero.

It is UNAVAILABLE, with `n_hyper` positive and `C` `NULL`, where the
criterion's own Hessian cannot be read: over a shared hyperparameter,
whose curvature would be that of the wrong function, which is the gap
[`statmod_hyper_vcov()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_vcov.md)
refuses for the same reason, and where the search left a coordinate at
the edge of its range.

It is PARTIAL, with `complete` false, where some of it could be read and
some could not: a hyperparameter
[`hyper_variance()`](https://statmodels7.github.io/statmodels7/reference/hyper_variance.md)
held contributes nothing, and so does a penalty over a structural term's
own parameters. The matrix returned is then a lower bound on the
correction rather than the whole of it.

## References

Wood, S. N., Pya, N. and Safken, B. (2016). Smoothing parameter and
model selection for general smooth models. *Journal of the American
Statistical Association*, 111(516), 1548–1563.

## See also

[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md),
[`statmod_hyper_vcov()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_vcov.md)
