# The Values One Entry of `par` Stands For

A numeric vector of the length the key addresses, from a vector, a
single number or a function of the count.

## Usage

``` r
par_values(e, k, key)
```

## Arguments

- e:

  The entry.

- k:

  How many values the key addresses.

- key:

  The name, for the message.

## Value

A numeric vector of length `k`.

## Details

A function answering with the wrong count is reported rather than
recycled: R would recycle it without a word and the simulation would be
of another model.

## See also

[`rstatmod_named()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_named.md)
