# The Columns of the Linear Part That Involve a Variable

Returns the positions, within the coefficients of a linear term, of the
columns whose model-matrix term involves `var`, read from the `assign`
attribute of the model matrix of `data`.

## Usage

``` r
linpar_columns(tm, var, data)
```

## Arguments

- tm:

  A
  [modelterms7::LinparTerm](https://statmodels7.github.io/modelterms7/reference/LinparTerm.html).

- var:

  A single string.

- data:

  The data frame of the fit.

## Value

An integer vector.
