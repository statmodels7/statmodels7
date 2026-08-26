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
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- data:

  The rows to read the model on.

- need_response:

  `TRUE` where the response must be present, as for a log-likelihood;
  `FALSE` for a prediction.

## Value

A
[`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
carrying the fitted terms, with `newdata` set to `data`, `n_obs` its row
count, and the offsets and response evaluated there. Every other
property is `spec`'s.

## Details

A term records how its block was made: a factor's levels and contrasts,
a spline's knots, a basis reparametrization.
[`modelterms7::term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.html)
reapplies that record to new rows.

Rebuilding instead gives a block of the same shape, multiplying the same
coefficients, that means something else. Measured on `y ~ s(x, k = 10)`
at 200 observations: predicting on 40 of the rows the model was fitted
to differed from the fitted values there by 0.237, and on the 51 rows
with \\\|x\| \< 0.5\\, where the rebuilt knots move furthest, by 1.19.
Handing back the whole data agrees exactly, which is why nothing
noticed.

The offsets are re-evaluated against `data` rather than carried across,
since a vector of the fitting data's length says nothing about other
rows.

## See also

[`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md),
which reapplies the terms,
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md),
the caller.
