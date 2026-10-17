# What a Hyperparameter's Own Uncertainty Adds to a Variance

Returns \\J V\_\theta J'\\, the movement of the penalized mode under the
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

  The bayesian variance to which the correction is added, over the
  coordinates that `keep` marks.

- keep:

  A logical vector over the joint vector (the stacked coefficients
  followed by any structural tail) that marks the coordinates that `Vb`
  spans. It is the whole vector of the caller and not its head, because
  a tail coordinate dropped as a flat direction leaves `Vb` narrower
  than the count alone would indicate.

- nz:

  How many of those coordinates are the structural tail.

## Value

A list with `C`, the correction or `NULL`, `n_hyper`, the number of
hyperparameters that
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
estimated, and `complete`, whether every one of them contributed and the
function skipped no penalty.

## Details

The matrices that
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
returns for `type = "bayesian"` and `type = "frequentist"` are
conditional on the hyperparameters: each is read at the value at which
the outer search stopped, as though that value were known. The value was
estimated from the same data, and the cost of that is \$\$V' = V_b + J
V\_\theta J', \qquad J = -(H + S)^{-1} \frac{\partial^2
\rho}{\partial\beta \partial\theta},\$\$ the delta method applied to the
map from the hyperparameter to the mode (Wood, Pya and Safken, 2016).
The same quantity is contracted against the information by
[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md)
to obtain a count of parameters; this function returns the matrix.

\\V_b\\ is passed in and is not recomputed. The caller has already
settled which information \\H\\ is, which coordinates are held and which
are aliased, and a correction built on a second inverse, regularized
differently, would not describe the movement of the mode to whose
variance it is added.

## Availability of the correction

The function returns `C = NULL` with `n_hyper = 0` where no
hyperparameter was estimated by
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md).
That includes a fit whose outer criterion is
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
and a kinked penalty, whose hyperparameter is the argument of a minimum
over a grid that
[`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
skips. The map from that hyperparameter to the mode turns a corner
whenever a coefficient joins or leaves the active set, so the map has no
derivative to propagate.

The correction is unavailable, and then `n_hyper` is positive and `C` is
`NULL`, where the criterion's own Hessian cannot be read, or where
[`hyper_variance()`](https://statmodels7.github.io/statmodels7/reference/hyper_variance.md)
cannot invert it because the coupling of a held coordinate to the kept
ones is too large to ignore.

The correction is partial, and `complete` is false, where some of it
could be read and some could not: a hyperparameter that
[`hyper_variance()`](https://statmodels7.github.io/statmodels7/reference/hyper_variance.md)
held contributes nothing, and a penalty for which
[`hyper_mode_cross()`](https://statmodels7.github.io/statmodels7/reference/hyper_mode_cross.md)
could not return a column (counted in its `skipped`) contributes
nothing. The matrix returned then omits those contributions.

A penalty over the own parameters of a structural term does not make the
correction partial. The mode of such a model moves in the joint vector
of the coefficients and the free parameters of the term, which `Vb`
spans, and
[`hyper_mode_cross()`](https://statmodels7.github.io/statmodels7/reference/hyper_mode_cross.md)
is called with `joint = TRUE` to return rows for the whole vector.

## References

Wood, S. N., Pya, N. and Safken, B. (2016). Smoothing parameter and
model selection for general smooth models. *Journal of the American
Statistical Association*, 111(516), 1548–1563.

## See also

[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md),
[`statmod_hyper_vcov()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_vcov.md)
