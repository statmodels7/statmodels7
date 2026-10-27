# The Status of One Refitted Replicate

Returns 0 for a replicate that is kept, 2 when its mode error exceeds
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md)
or cannot be read, and 3 when it leaves a coordinate unidentified that
the fit identifies, or has an information that cannot be factorized
where the fit's can.

## Usage

``` r
boot_validate(object, spec_b, r, cfg, ref)
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

## Value

A single integer.
