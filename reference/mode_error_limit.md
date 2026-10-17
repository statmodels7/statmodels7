# How Far Above Its Minimum an Inner Fit May Sit and Still Report a Resolution

Returns the limit that
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
places on \\g'K^{-1}g/2\\, the decrease that the mode's own Newton
correction predicts. Above this limit,
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
reports no resolution.

## Usage

``` r
mode_error_limit()
```

## Value

A single number.

## Details

The resolution is read by displacing the coefficients to the mode
implied by the inner score and measuring how far the criterion moves.
That is a resolution while the displacement is a correction and
something else entirely once it is not: an inner fit that stopped far
from a mode produces a large displacement, a large movement, and a
number that
[`optimizers7::crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.html)
would read as license to stop.

The quantity tested is in log-likelihood units and not in the units of
the coefficients, so one limit serves every model.

## See also

[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
