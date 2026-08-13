# Does a Term Supply Its Third Derivative?

Whether the term implements
[`term_third`](https://statmodels7.github.io/modelterms7/reference/term_third.html),
read from the class the method is registered on rather than from a list
of class names, so a term written later is covered without an edit here.

## Usage

``` r
answers_term_third(term)
```

## Arguments

- term:

  A built term.

## Value

A single logical.

## Details

A structural term that has not written it inherits a method that signals
an error, and an additive term inherits one that returns zero, which is
the right answer for a predictor that is a block of columns. The
question is therefore whether the owning class is the refusing base.

## See also

[`outer_gradient_ok`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
