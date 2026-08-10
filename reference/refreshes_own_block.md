# Does a Term Recompute Its Own Block?

`TRUE` when the term registers a
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.html)
method of its own rather than inheriting the identity registered on
`model_term`.

## Usage

``` r
refreshes_own_block(term)
```

## Arguments

- term:

  One built term.

## Value

A single logical.

## Details

The owning class of a method is `attr(m, "signature")[[1]]`, and it is
compared by name and package rather than by
[`identical()`](https://rdrr.io/r/base/identical.html): an S7 class
re-created from the same definition is not identical to the original,
which is what happens whenever a package's code is re-evaluated rather
than loaded.

## See also

[`unfittable_reason`](https://statmodels7.github.io/statmodels7/reference/unfittable_reason.md)
