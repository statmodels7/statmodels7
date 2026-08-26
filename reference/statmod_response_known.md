# Whether New Rows Carry the Response

`TRUE` where every row has one, `FALSE` where none does, and an error
where some do.

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

It is what separates a re-reading of a model on another series from a
continuation of the one it was fitted to, and the separation has to be
all-or-nothing: a filter's recursion at one row reads the rows before
it, so a frame carrying the response on some rows only describes neither
operation, and answering it would mean choosing a reading the caller did
not ask for.

## See also

[`statmod_eta_continued()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta_continued.md)
