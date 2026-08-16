# Add Two Sets of Offsets

Combines the offsets a formula names with those the `offsets` argument
supplies, per parameter.

## Usage

``` r
add_offsets(a, b)
```

## Arguments

- a, b:

  Two named lists of offsets, either entry possibly `NULL`.

## Value

A named list.

## Details

They are SUMMED where both are given, which is what
[`glm`](https://rdrr.io/r/stats/glm.html) does with a formula offset and
an `offset` argument together.
