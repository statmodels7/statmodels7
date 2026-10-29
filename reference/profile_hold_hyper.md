# Hold Hyperparameters on a Specification

Writes the hyperparameters of the rows `hyp` of the table of
[`profile_hyper_table()`](https://statmodels7.github.io/statmodels7/reference/profile_hyper_table.md),
at the values `eta` on the scale of their links, into `held_hyper` of
the specification, under the key of every member of a shared
hyperparameter.

## Usage

``` r
profile_hold_hyper(spec, hyp, eta)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- hyp:

  The rows of the table, a data frame with the attributes of
  [`profile_hyper_table()`](https://statmodels7.github.io/statmodels7/reference/profile_hyper_table.md)
  and a column `row`, the row of the full table.

- eta:

  The values on the scale of the links, one for each row of `hyp`.

## Value

`spec`, with `held_hyper` set.
