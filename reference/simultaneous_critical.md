# The Critical Value of a Simultaneous Band

Returns the `level` quantile of \\\max_x \|D(x) e\| / \mathrm{se}(x)\\
over the draws \\e\\, with the maximum over the rows of `D` whose
standard error is positive.

## Usage

``` r
simultaneous_critical(D, se, draws, level)
```

## Arguments

- D:

  The matrix of the rows of the derivative of the block.

- se:

  The standard errors of the rows.

- draws:

  The draws of the coefficients, one in each row.

- level:

  The level.

## Value

A single number.
