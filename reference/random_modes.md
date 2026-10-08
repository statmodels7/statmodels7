# The Random-Effect Terms a Prediction Can Set Aside

Lists the
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
terms written in the model's equations, with the mode
[`predict()`](https://rdrr.io/r/stats/predict.html) reads each at.

## Usage

``` r
random_modes(spec, random)
```

## Arguments

- spec:

  The fit's specification.

- random:

  The mode or modes, as
  [`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
  takes them.

## Value

A data frame with columns `param`, `key` and `mode`, one row per
random-effect term, or zero rows where the model has none, with the
nested terms in the attribute `"nested"`.

## Details

`random` is a single string, applied to every such term, or a character
vector named by the terms' keys, one mode per term named, the rest
staying conditional. A key is matched with its white space removed, so
`"random(~1|g)"` finds the term written `random(~ 1 | g)`. A name that
matches no term signals an error listing those there are.

Only a term written in an equation is listed. One written inside a
subformula develops another term's parameter; those are listed by
[`nested_random_terms()`](https://statmodels7.github.io/statmodels7/reference/nested_random_terms.md)
and returned in the attribute `"nested"`, with their modes. An unnamed
mode applies to them as well, and a name may be the key of either kind.
