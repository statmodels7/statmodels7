# The Leverage Block of Each Observation

Computes \\H\_{ii} = W_i^{1/2} Z_i V_b Z_i' W_i^{1/2}\\ for every
observation, from the pieces of
[`leverage_parts()`](https://statmodels7.github.io/statmodels7/reference/leverage_parts.md)
and the inverse penalized information \\V_b\\.

## Usage

``` r
leverage_blocks(parts, Vb, chunk = 4000L)
```

## Arguments

- parts:

  The list that
  [`leverage_parts()`](https://statmodels7.github.io/statmodels7/reference/leverage_parts.md)
  returns.

- Vb:

  The inverse penalized information over the coordinates of `parts$Z`, a
  dense matrix.

- chunk:

  How many observations to process at a time.

## Value

An array of dimension \\n \times d \times d\\.
