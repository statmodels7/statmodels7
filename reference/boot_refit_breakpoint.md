# Refit One Replicate of a Model With a Break-Point Term

Refits the model on the replicate from the break-points and coefficients
of the fit, and then runs the restart cycle of
[`statmod_boot_restart()`](https://statmodels7.github.io/statmodels7/reference/statmod_boot_restart.md)
and the settling of the sharp break-points, as
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
does at the end of a fit.

## Usage

``` r
boot_refit_breakpoint(object, spec_b, refit, method, cfg, beta0)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec_b:

  The specification of the replicate.

- refit, method, cfg, beta0:

  As in
  [`bootstrap_refit()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_refit.md).

## Value

A list with `coefficients`, `hyper` and `design`.

## Details

The refit from the fit's positions is
[`bootstrap_refit()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_refit.md).
A start from the positions of the fit alone would find the local optimum
nearest to them, which is another estimator than the one that produced
the fit; the restart cycle searches for a better optimum of the
replicate with the budget `n_boot` that the terms set, and the settling
puts a sharp break-point at the minimum of its profile. The
hyperparameters are held at the fit's values with
`refit = "coefficients"` and selected again with `"full"`, before the
restart.
