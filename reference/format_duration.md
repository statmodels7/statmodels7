# Format a Duration in the Unit It Deserves

Turns a number of seconds into a string in microseconds, milliseconds,
seconds, minutes, hours or days, whichever keeps it readable.

## Usage

``` r
format_duration(seconds, digits = 3L)
```

## Arguments

- seconds:

  A number of seconds.

- digits:

  Significant digits.

## Value

A single string.

## Details

A fit that took 340 microseconds and one that took 2.6 hours must both
read at a glance, which a single unit cannot do.

Zero is a reading and not the absence of one: on a coarse clock a fast
fit measures exactly zero seconds, so it is reported as such rather than
suppressed. A guard of the form `if (elapsed > 0)` would be a claim that
zero cannot be measured, and for a duration that claim is false.

## Examples

``` r
format_duration(c(3.4e-4, 0.25, 90, 7200))
#> [1] "340 us"  "250 ms"  "1.5 min" "2 h"    
```
