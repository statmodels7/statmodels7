# Whether the Default Search Steers by the Exact Hessian

Whether a search this package chooses should take
[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)'s
order-2 answer as its direction, which is a narrower question than
whether the Hessian exists.

## Usage

``` r
outer_newton_ok(spec, design, exact2)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, from
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md).

- exact2:

  Whether the criterion has an exact Hessian at all.

## Value

`TRUE` or `FALSE`.

## Details

The two are separated because a criterion may have an exact second
derivative that is a poor thing to steer by. Over a COVARIANCE CLASS the
criterion's Hessian is strongly indefinite wherever the chart's angle
approaches its boundary – measured on ten groups of ten whose truth
carries a correlation of exactly one, the eigenvalues of −H at the point
the search reports are 1.13e+07, -1.48e+10 and -4.20e+17 – so `newton()`
takes its eigen-floor branch and the repaired step lands elsewhere. On
that panel it ends at an angle of -11.35 against -8.14, at a criterion
of -127.6696 against -126.3168, reporting `not converged` where the
shorter run reports `boundary`, and there the whole outer curvature is
unreadable so NO coordinate keeps a standard error rather than only the
correlation. Away from the boundary the same shape is fine and merely
dearer: on twenty groups `newton()` reaches the same criterion in 32
evaluations against 85 and 13.4 seconds against 7.0.

What the Hessian's availability still governs is untouched: it is
supplied to an optimizer the caller NAMES, it is what
[`statmod_hyper_vcov()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_vcov.md)
and `vcov(type = "unconditional")` read, and it is recorded as
`exact_hessian`. Only the default direction is held back, and only for
the shape the measurement names.

The block width is asked of the penalty with
[`S7::prop_names()`](https://rconsortium.github.io/S7/reference/prop_names.html)
rather than assumed, a branch that does not carry one being univariate.

## See also

[`outer_default_optimizer()`](https://statmodels7.github.io/statmodels7/reference/outer_default_optimizer.md),
[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
