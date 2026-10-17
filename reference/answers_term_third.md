# Does a Term Supply Its Third Derivative?

Determines whether the term implements
[`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html)
from the class on which the method is registered, never from a list of
class names, so a term written later is covered without an edit here.

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

A structural term that has not implemented the method inherits one that
signals an error, and an additive term inherits one that returns zero,
which is correct for a predictor that is a block of columns. The
function therefore returns `FALSE` where the class that owns the method
is the base class `structural_term`.

## See also

[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
