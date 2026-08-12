# The Parameters Under One Regime

The distribution's parameters at every observation with the regime
term's equation shifted by one state's level.

## Usage

``` r
statmod_theta_shifted(spec, eta_static, param, shift)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- eta_static:

  The static predictors.

- param:

  The equation the term sits in.

- shift:

  The level to add.

## Value

A named list of per-observation parameters.
