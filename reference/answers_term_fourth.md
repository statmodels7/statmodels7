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
recursion, so a term that implements the third derivative and not the
fourth supplies an exact gradient and no exact Hessian. For such a term
the function returns `FALSE` and the Hessian is computed by
[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md).
A term of the likelihood shape (`regime()` and the marginal break-point
term) implements neither derivative, and
[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
returns `FALSE` for it at both orders.

## See also

[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md),
[`statmod_structural_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_hess.md)
