# The Expected Information of Each Predictor at the Intercept-Only Fit

For each distribution parameter, the mean over the observations of the
expected information of its linear predictor, \\-\mathrm{E}\[\partial^2
\ell / \partial \eta^2\]\\, at the parameters of the intercept-only fit.

## Usage

``` r
predictor_information(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A named numeric vector, one entry per distribution parameter, or `NULL`
where the intercept-only fit is not available.

## Details

It sets the scale of a score-driven term's second start
([`structural_start_fixups()`](https://statmodels7.github.io/statmodels7/reference/structural_start_fixups.md)):
the score of a predictor has variance equal to this information, so a
loading of \\0.1/\mathcal{I}\\ moves the predictor by a tenth of the
score's natural unit whatever the family.
