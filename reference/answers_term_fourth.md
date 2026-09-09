# Does a Term Supply Its Fourth Derivative?

[`answers_term_third()`](https://statmodels7.github.io/statmodels7/reference/answers_term_third.md)'s
question one order up, and read the same way: from the class the method
is registered on, never from a list of class names.

## Usage

``` r
answers_term_fourth(term)
```

## Arguments

- term:

  A built term.

## Value

A single logical.

## Details

The criterion's own second derivative reads a fourth order through the
recursion, so a term that has written the third and not the fourth
supplies an exact gradient and no exact Hessian. That is the state
`regime()` is in: it implements neither, and a term implementing only
the third would be answered here with `FALSE` and left to
[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md),
which is the same fallback its search already has.

## See also

[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md),
[`statmod_structural_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_hess.md)
