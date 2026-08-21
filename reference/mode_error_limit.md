# How Far Above Its Minimum an Inner Fit May Sit and Still Report a Resolution

The limit
[`criterion_resolution`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
puts on \\g'K^{-1}g/2\\, the decrease the mode's own Newton correction
predicts, before it refuses to report a resolution at all.

## Usage

``` r
mode_error_limit()
```

## Value

A single number.

## Details

The resolution is read by displacing the coefficients to where the inner
score says the mode is and asking how far the criterion moved. That is a
resolution while the displacement is a correction and something else
entirely once it is not: an inner fit that stopped far from a mode
produces a large displacement, a large movement, and a number that
[`crit_abs_obj`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.html)
would read as licence to stop.

The quantity tested is in log-likelihood units rather than in the
coefficients', so one limit serves every shape. Measured over whole fits
it reaches 2.8e-08 on a smooth and 1.4e-10 on a random intercept,
against 20.9 to 22.9 on a hierarchical break-point model whose inner fit
reports convergence at a score of 247.8. Nine orders separate them, so
the value below has five orders of room on either side and is not what
decides anything.

## See also

[`criterion_resolution`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
