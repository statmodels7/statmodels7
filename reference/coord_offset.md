# The Offset of One Equation

Returns the offset of one distribution parameter's equation, evaluated
in the fitting data and recycled to the sample size. The function
returns zeros for an equation with no offset, so the result can be added
unconditionally.

## Usage

``` r
coord_offset(spec, p, n)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md),
  read for its `offsets` list.

- p:

  Which distribution parameter, a string naming one of the family's.

- n:

  The number of observations to which the offset is recycled.

## Value

A numeric vector of length `n`: the offset, or zeros where the equation
has none.

## Details

The offsets are stored on the specification as evaluated vectors, one
per parameter, with `NULL` for an equation that has none. This function
turns that `NULL` into zeros at the point of use, so callers do not have
to test for it.

An offset shorter than `n` is recycled with
[`rep_len()`](https://rdrr.io/r/base/rep.html), so a single number is a
constant offset. That is the shape produced by
`offsets = list(mu = log(2))`.

## See also

[`eval_offsets()`](https://statmodels7.github.io/statmodels7/reference/eval_offsets.md),
which evaluates the expressions that this reads,
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
for the `offsets` argument.
