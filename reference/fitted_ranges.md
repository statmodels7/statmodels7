# What Each Distribution Parameter Reached

One line per parameter of the distribution, giving the range of its
fitted values, for a fit that did not converge.

## Usage

``` r
fitted_ranges(x)
```

## Arguments

- x:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

## Value

A single string, empty when the parameters cannot be read. It is a note
of
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
rather than a line of
[`print.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/print.StatmodFit.md):
it qualifies the fit rather than describing it, and it is read once when
something looks wrong.

## Details

The measured case this exists for: a lasso at a fixed hyperparameter,
with a free scale and a design the model can interpolate. Fitting the
coefficients shrinks the residuals, which shrinks the scale, which
raises the working weights, which makes the penalty count for relatively
less, which lets more coefficients in. At 200 observations and 400
columns the scale reached 3.8e-15 and 380 of the 400 coefficients
survived, where the same block fitted at a held scale kept the five that
were real.

Nothing here diagnoses that. It reports where the parameters ended up,
which is a fact, and a scale at 1e-15 says the rest on its own. Naming a
cause would mean picking a threshold for what counts as running away,
and the same fit at 100 columns converges to a scale of 0.77 that is
nothing of the kind.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
