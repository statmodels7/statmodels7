# How Far From Its Mode a Fit Is, and What It Does Not Identify

Returns \\\tfrac12 g'K^{-1}g\\ at the coefficients, as the certificate
reads it, and the coordinates that the penalized information does not
identify there.

## Usage

``` r
boot_read_mode(
  object,
  spec,
  design,
  cf,
  hy,
  cfg,
  pinned = object@methods$pinned
)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md),
  for its criterion and pinned coordinates.

- spec, design:

  The specification and design where the point is read.

- cf, hy:

  The coefficients and hyperparameters.

- cfg:

  The inner settings of the fit.

- pinned:

  The pinned coordinates at the point read, in its own stacked
  positions; those of the fit by default.

## Value

A list with `me` (a number, `NA` if the point could not be read), `def`
(stacked positions) and `singular` (whether the penalized information
could not be factorized).

## Details

The coefficients that the marginal criterion estimates are held, as the
certificate holds them, and the coordinates of a kinked penalty that are
zero are left out of the score. The information is that of the inner
fit. The unidentified coordinates are those that
[`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md)
pins, apart from the ones pinned for another reason: a coefficient whose
observations all have a zero weight is among them, although the inner
fit does not report it as aliased.
