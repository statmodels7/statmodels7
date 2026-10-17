# A Vector of a Simulated Truth, Short Enough to Read

Returns the values where there are six or fewer, and otherwise the first
two values with the count and the standard deviation.

## Usage

``` r
sim_values(v)
```

## Arguments

- v:

  A numeric vector.

## Value

A single string.

## Details

A random effect over a hundred groups is a hundred numbers, and printing
them hides the few numbers that a reader needs. For a block of that size
the printed form gives the count and the spread.

## See also

[`print.StatmodSim()`](https://statmodels7.github.io/statmodels7/reference/print.StatmodSim.md)
