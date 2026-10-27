# The Status of One Refitted Cluster Replicate With Regrouped Effects

Returns 0 for a replicate that is kept, 2 when its mode error exceeds
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md)
or cannot be read, and 3 when it leaves unidentified a coordinate that
the fit identifies.

## Usage

``` r
boot_cluster_validate(object, spec_b, r, cfg, ref, cb)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec_b:

  The specification of the replicate.

- r:

  The refit, as
  [`bootstrap_refit()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_refit.md)
  returns it.

- cfg:

  The inner settings of the fit.

- ref:

  The reading of the fit itself, from
  [`boot_read_mode()`](https://statmodels7.github.io/statmodels7/reference/boot_read_mode.md).

- cb:

  The result of the `build` function of
  [`boot_cluster_random_plan()`](https://statmodels7.github.io/statmodels7/reference/boot_cluster_random_plan.md).

## Value

A single integer.

## Details

It is
[`boot_validate()`](https://statmodels7.github.io/statmodels7/reference/boot_validate.md)
with the stacked positions translated: the pinned coordinates and the
unidentified ones of the fit are mapped to the positions of the
replicate, and only the coordinates that the replicate shares with the
fit are compared. The effects of the regrouped terms have no counterpart
in the fit and are not compared.
