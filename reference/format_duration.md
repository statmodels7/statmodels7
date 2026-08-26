# Format a Duration in the Unit It Deserves

Turns a number of seconds into a string in microseconds, milliseconds,
seconds, minutes, hours or days, whichever keeps it readable.

## Usage

``` r
format_duration(seconds, digits = 3L)
```

## Arguments

- seconds:

  A number of seconds. Vectorized: a vector in gives a vector out, each
  element in its own unit.

- digits:

  How many significant digits to keep, passed to
  [`signif()`](https://rdrr.io/r/base/Round.html). Defaults to 3.

## Value

A character vector the length of `seconds`, each element a number and a
unit separated by a space, as `"340 us"` or `"2.6 h"`.

## Details

A fit that took 340 microseconds and one that took 2.6 hours must both
read at a glance, and no single unit does that.

The thresholds are the obvious ones: below a millisecond it reads in
microseconds, below a second in milliseconds, below a minute in seconds,
then minutes, hours and days.

Zero is a reading. On a coarse clock a fast fit measures exactly zero
seconds, and it is reported as `"0 s"` rather than suppressed. A guard
of the form `if (elapsed > 0)` would claim that zero cannot be measured,
which for a duration is false, and would make what the object prints
depend on the platform's timer resolution.

## See also

[`print.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/print.StatmodFit.md),
which reports a fit's elapsed time with this.

## Examples

``` r
format_duration(c(3.4e-4, 0.25, 90, 7200))
#> [1] "340 us"  "250 ms"  "1.5 min" "2 h"    

# Zero is a reading, not a missing value.
format_duration(0)
#> [1] "0 s"
```
