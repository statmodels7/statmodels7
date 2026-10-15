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
estimated some coefficients itself (by default
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
estimate on the criterion the unpenalized coefficients of every
parameter but the first, see
[`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md)),
the point that the fit returned is not the joint mode, and
`fit@objective` is not its minimum. A likelihood ratio differenced
against `fit@objective` then contains the difference between the
objective at the fitted point and at the joint mode.

The full model is therefore refitted with nothing held, from the fitted
coefficients, which are already close to the joint mode, at the fitted
hyperparameters. Where the criterion estimated nothing, the fit is at
the joint mode and `fit@objective` is returned unchanged.

## See also

[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md),
[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md).
