# Which Hyperparameters the Terms Hold

Returns one key per hyperparameter that a term fixed in its constructor,
composed of `parameter`, the penalty's key and the hyperparameter's name
joined by a carriage return.

## Usage

``` r
statmod_held(spec, design = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

## Value

A character vector, possibly empty.

## Details

Which hyperparameters are estimated is a property of the terms and not
of the criterion, because the term is where the penalty is named.
[`path_rows()`](https://statmodels7.github.io/statmodels7/reference/path_rows.md)
and the summary's account of what was estimated and what was held read
this enumeration. The outer index reads the same `fixed` entries of
[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md).

## See also

[`modelterms7::term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.html),
[`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
