# Add Two Sets of Offsets

Combines the offsets a formula names with those the `offsets` argument
supplies, per parameter.

## Usage

``` r
add_offsets(a, b)
```

## Arguments

- a, b:

  Two named lists of offsets keyed by distribution parameter, either
  list's entries possibly `NULL`. Either list may itself be `NULL`.

## Value

A named list with one entry per parameter named in either input, `NULL`
where neither supplied one.

## Details

Where both name an offset for the same parameter the two are **summed**,
as [`stats::glm()`](https://rdrr.io/r/stats/glm.html) does with a
formula offset and an `offset` argument together. Where only one does,
that one is taken.

## See also

[`eval_offsets()`](https://statmodels7.github.io/statmodels7/reference/eval_offsets.md)
for the formula's own,
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
for the `offsets` argument.
