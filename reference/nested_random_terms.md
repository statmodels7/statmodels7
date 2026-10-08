# The Random-Effect Terms Written Inside a Subformula

Lists the
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
sub-terms that develop a parameter of another term, such as
`Asym ~ 1 + random(~ 1 | Tree)` inside
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
or `psi ~ random(~ 1 | id)` inside
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html).

## Usage

``` r
nested_random_terms(spec)
```

## Arguments

- spec:

  The fit's specification.

## Value

A data frame with columns `param`, `term`, `comp`, `s` (the sub-term's
place in the parameter's subformula) and `key`, one row per nested
random-effect sub-term, or zero rows.

## Details

The terms are read through
[`modelterms7::term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.html),
which every term carrying developed parameters answers. The key of a row
is the outer term's key, the parameter and the sub-term joined by
`"::"`, which is the key
[`hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
reports where a term carries more than one penalty.
