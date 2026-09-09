# Which Structural Term a Model Carries, If Any

The single term of the structural branch, or `NULL` where the model
carries none.

## Usage

``` r
structural_term_of(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

## Value

A built term, or `NULL`.

## Details

At most one structural term is admitted per formula, so a caller wanting
to ask that term a question – does it answer
[`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html),
does it answer
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html)
– has one to ask. It is looked up through the design's structural
attribute rather than by walking the terms, which is where the fit
records what it built.
