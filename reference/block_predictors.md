# The Predictors a Coefficient Direction Induces

\\(X_k v_k)\_i\\ for each distribution parameter, which is how a
movement of the coefficients is felt by the log-density.

## Usage

``` r
block_predictors(design, params, npar, offs, v)
```

## Arguments

- design:

  The design.

- params:

  The parameter names.

- npar, offs:

  The block sizes and offsets.

- v:

  A stacked coefficient vector.

## Value

A list of numeric vectors, one per parameter.
