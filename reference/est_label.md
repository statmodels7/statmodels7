# The Name of an Estimate and Its Interval in a Legend

Writes the name of an estimate followed by its interval on a second
line, as in "estimated effect\n(95% CI)", for a legend entry that draws
the estimate and its band together.

## Usage

``` r
est_label(what, level, interval = "CI")
```

## Arguments

- what:

  The name of the estimate.

- level:

  A number in (0, 1).

- interval:

  The name of the interval.

## Value

A single string.
