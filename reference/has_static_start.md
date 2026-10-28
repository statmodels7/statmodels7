# Whether a Term Carries a Statically Started Filter

Whether a Term Carries a Statically Started Filter

## Usage

``` r
has_static_start(tm)
```

## Arguments

- tm:

  A model term.

## Value

A single logical: `TRUE` for a
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
term with `init = "static"`, or a joint filter with such a member.
