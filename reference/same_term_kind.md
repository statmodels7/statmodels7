# Whether Two Terms Are of the Same Kind

`TRUE` when two built terms have the same class, which is what keeps the
stem of
[`block_stem()`](https://statmodels7.github.io/statmodels7/reference/block_stem.md)
from pairing a smooth with something else that happens to name its
coefficients the same way.

## Usage

``` r
same_term_kind(a, b)
```

## Arguments

- a, b:

  Two built terms, either of which may be `NULL`.

## Value

A single logical.

## See also

[`block_stem()`](https://statmodels7.github.io/statmodels7/reference/block_stem.md),
the other half of the pairing.
