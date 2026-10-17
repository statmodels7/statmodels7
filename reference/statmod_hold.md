# A Specification With Coefficients Held

Writes values for the coefficients that
[`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md)
returned into the specification's `held_coef`, which the inner fit
enforces.

## Usage

``` r
statmod_hold(spec, coords, values)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coords:

  The result of
  [`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md).

- values:

  The values, one per position in `coords$where`.

## Value

`spec`, in which `held_coef` is set. The holds already present are kept
for coefficients that `coords` does not name.
