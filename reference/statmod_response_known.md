# Whether New Rows Carry the Response

Returns `TRUE` where every row has the response, `FALSE` where every row
lacks it, and signals an error where some rows have it and others do
not.

## Usage

``` r
statmod_response_known(y)
```

## Arguments

- y:

  The response as the specification carries it.

## Value

`TRUE` or `FALSE`.

## Details

The distinction separates a re-reading of a model on another series from
a continuation of the series to which the model was fitted, and it has
to be all-or-nothing. The recursion of a filter at one row reads the
rows before it, so a frame that carries the response on some rows only
describes neither operation, and any answer would require choosing one
of the two readings for the caller.

## See also

[`statmod_eta_continued()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta_continued.md)
