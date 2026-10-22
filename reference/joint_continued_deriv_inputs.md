# The Inputs to the Derivative of a Continued Joint Filter

Returns the derivatives of the levels and of the driving scores at the
observed rows, and of the filter's parameters, in the coordinates of the
variance matrix, as
[`modelterms7::term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.html)
takes them for a joint filter.

## Usage

``` r
joint_continued_deriv_inputs(ospec, odesign, coef, ost)
```

## Arguments

- ospec, odesign:

  The fit's specification and design.

- coef:

  The coefficients.

- ost:

  The structural state.

## Value

A list with `deriv`, `cols` (for each filtered parameter, the columns of
its equation's coefficients) and `key`.

## Details

The level of a filtered equation is its predictor less the static part,
so its derivative row is the equation's row less its design. The driving
score of parameter \\l\\ depends on the predictors of every equation, so
its derivative row is \\\sum_b U\_{l;b} V_b\\.

## See also

[`statmod_eta_continued()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta_continued.md)
