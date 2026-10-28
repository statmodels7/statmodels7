# The Size of a Negative Eigenvalue That Is Not Rounding

Returns the factor \\c\\ such that a curvature whose smallest eigenvalue
is below \\-c\\\varepsilon\\\lambda\_{\max}\\ is read as having a
direction of ascent, with \\\varepsilon\\ the machine epsilon and
\\\lambda\_{\max}\\ the largest eigenvalue.

## Usage

``` r
curvature_ascent_limit()
```

## Value

A single number.

## Details

[`joint_reading()`](https://statmodels7.github.io/statmodels7/reference/joint_reading.md)
applies it to the joint curvature of the coefficients and a structural
term's parameters. An information that is positive semidefinite and
singular, such as one with a duplicated column, has a smallest
eigenvalue that rounding makes slightly negative, and the factor is
large enough for that eigenvalue not to be read as a direction of
ascent.

## See also

[`joint_reading()`](https://statmodels7.github.io/statmodels7/reference/joint_reading.md)
