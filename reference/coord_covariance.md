# Which Way of Holding the Gradient Is Cheaper

`TRUE` for the covariance form, `FALSE` for the residual.

## Usage

``` r
coord_covariance(n, m)
```

## Arguments

- n:

  The number of observations.

- m:

  How many coordinates are visited.

## Value

A single logical.

## Details

The covariance form replaces an \\O(n)\\ read with an \\O(m)\\ one, and
pays for it by building a column of \\X'WX\\ at \\O(nm)\\ the first time
a coordinate moves off zero. It is therefore worth it only when \\m\\ is
small next to \\n\\, and the measurement is unambiguous: at 5000
observations with 200 columns screened to 200, taking the covariance
form cost 70 milliseconds against 55 for the residual, the Gram columns
being dearer than the residual passes they replaced. The threshold is
set where the two cross rather than at a rule of thumb.
