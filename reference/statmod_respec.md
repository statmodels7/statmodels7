# The Same Model Read on Other Rows

A specification carrying the fitted terms and a new data frame, so that
every block is reapplied to those rows rather than rebuilt from them.

## Usage

``` r
statmod_respec(spec, data, need_response = TRUE)
```

## Arguments

- spec:

  The fitted
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- data:

  The rows to read it on.

- need_response:

  Whether the response has to be there.

## Value

A
[`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
whose `newdata` is set.

## Details

A term records how its block was made – a factor's levels and contrasts,
a spline's knots, a basis reparametrization – and
[`term_predict`](https://statmodels7.github.io/modelterms7/reference/term_predict.html)
reapplies that record. Rebuilding instead gives a block of the same
shape, multiplying the same coefficients, that means something else:
measured on `y ~ s(x, k = 10)` at 200 observations, predicting on 40 of
the rows the model was fitted to differed from the fitted values there
by 0.237, and on the 51 rows with \\\|x\| \< 0.5\\, where the rebuilt
knots move furthest, by 1.19. The whole data handed back agrees exactly,
which is why nothing noticed.

## See also

[`statmod_design`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
