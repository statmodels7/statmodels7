# Which Hyperparameters the Terms Hold

One key per hyperparameter a term fixed in its constructor, as
`parameter`, the penalty's key and its name joined by a carriage return.

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

which hyperparameters are estimated is a property of the terms, not of
the criterion: the term is where the penalty is named, and a criterion
argument saying otherwise was read by nothing when the two disagreed.
Everything here consults this one enumeration: the outer index, the
path, and the summary's account of what was estimated and what was
given.

## See also

[`modelterms7::term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.html),
[`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
