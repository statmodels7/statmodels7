# The Penalized Objective at the Unrestricted Joint Mode

The minimum of the penalized objective over every coefficient at the
fitted hyperparameters, which is the unrestricted side of the likelihood
ratio
[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md)
computes.

## Usage

``` r
joint_objective(fit)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

## Value

A single number, on the scale of `fit@objective`.

## Details

The restricted fit of
[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md)
holds one coefficient and takes every other one to the joint mode of the
penalized likelihood. The unrestricted side has to be read at the same
kind of point, or the statistic compares two fits whose other
coefficients were chosen by different rules. Where the outer criterion
estimated some coefficients itself – by default,
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
estimates a dispersion's unpenalized coefficients on the criterion
([`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md))
– the point the fit returned is not the joint mode, and `fit@objective`
is not its minimum. Measured on a gaussian linear regression at \\n =
40\\, differencing against it gave a likelihood ratio of 1.0458 against
the 1.1643 that [`stats::lm()`](https://rdrr.io/r/stats/lm.html) gives,
and -0.118 at the estimate itself.

The full model is therefore refitted with nothing held, from the fitted
coefficients, which are already close to the joint mode, at the fitted
hyperparameters. Where the criterion estimated nothing the fit IS at the
joint mode and `fit@objective` is returned unchanged, so every such
statistic is what it was.

## See also

[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md),
[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md).
