# Does a Term Supply Its Fourth Derivative?

Determines whether the term implements
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html),
which is the question of
[`answers_term_third()`](https://statmodels7.github.io/statmodels7/reference/answers_term_third.md)
one order up, and reads it the same way: from the class on which the
method is registered, never from a list of class names.

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
