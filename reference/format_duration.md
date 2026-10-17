# Format a Duration in a Readable Unit

Converts a number of seconds to a string in microseconds, milliseconds,
seconds, minutes, hours or days, whichever keeps the number readable.

## Usage

``` r
format_duration(seconds, digits = 3L)
```

## Arguments

- seconds:

  A number of seconds. Vectorized: a vector in gives a vector out, and
  each element is in its own unit.

- digits:

  How many significant digits to keep, passed to
  [`signif()`](https://rdrr.io/r/base/Round.html). Defaults to 3.

## Value

A character vector of the length of `seconds`. Each element is a number
and a unit separated by a space, as `"340 us"` or `"2.6 h"`. A
non-finite value gives `"NA"` and a negative value gives `"negative"`.

## Details

Fit times range from a few hundred microseconds to several hours, so the
unit is chosen from the magnitude of the value.

The thresholds are the usual ones: a value reads in microseconds below a
millisecond, in milliseconds below a second, in seconds below a minute,
and then in minutes, hours and days. The unit is chosen after the value
is rounded to `digits` significant digits, so that a value just below a
threshold reads in the larger unit (59.99 seconds reads `"1 min"`).

A duration of zero is reported as `"0 s"`. On a coarse clock a fast fit
can measure exactly zero seconds, and suppressing that value would make
the printed output depend on the timer resolution of the platform.

## See also

[`print.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/print.StatmodFit.md),
which reports the elapsed time of a fit with this function.

## Examples

``` r
format_duration(c(3.4e-4, 0.25, 90, 7200))
#> [1] "340 us"  "250 ms"  "1.5 min" "2 h"    

# A duration of zero is reported in seconds, a missing value as "NA".
format_duration(c(0, NA))
#> [1] "0 s" "NA" 
```
