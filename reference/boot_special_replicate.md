# One Replicate of a Special Model

The body of the refit of a replicate for a model with a structural term
or a break-point term, as
[`boot_body()`](https://statmodels7.github.io/statmodels7/reference/boot_body.md)
runs it.

## Usage

``` r
boot_special_replicate(
  object,
  spec_b,
  pl,
  mode,
  refit,
  method,
  cfg,
  beta0,
  keep_fit
)
```

## Arguments

- object, spec_b, pl:

  The fit, the specification and the draw of the replicate.

- mode:

  `"structural"` or `"refresh"`.

- refit, method, cfg, beta0:

  As in
  [`bootstrap_refit()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_refit.md).

- keep_fit:

  Whether to return the refit.

## Value

A list as
[`boot_body()`](https://statmodels7.github.io/statmodels7/reference/boot_body.md)
returns.
