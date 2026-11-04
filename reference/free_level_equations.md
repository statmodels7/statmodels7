# The Free Levels That reml() Integrates Beside the Coefficients

Returns the equations whose structural level is free and is integrated
by
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
as one more column of that equation's design, or `character(0)` where
nothing is added.

## Usage

``` r
free_level_equations(spec, design, method = NULL, basis = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md),
  or `NULL` to skip the check of its kind.

- basis:

  The integrated subspace, `NULL` under
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md).

## Value

A character vector of distribution parameters, one per free level.

## Details

A score-driven level, the first level of
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
and each level of a joint filter add a constant to their equation's
predictor. Where an intercept is present the level is held
([`statmod_held_levels()`](https://statmodels7.github.io/statmodels7/reference/statmod_held_levels.md))
and the intercept is integrated as an unpenalized coefficient. Where the
level is free, it plays the intercept's part, and the criterion
integrates it in the same way: as the constant that it adds to the
predictor, a column of ones in its equation. For
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
that constant is the stationary level \\\omega / (1 - \sum_j \xi_j)\\,
so `y ~ x + random(~ 1 | g) + gas()` and
`y ~ 0 + x + random(~ 1 | g) + gas()` give the same criterion at every
hyperparameter.

The rule applies to
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
on the route over the coefficients, where no penalty covers the term's
own parameters; where one does,
[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
already integrates every free parameter of the term.
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
integrates the penalized directions only and is unchanged. A developed
level contributes its constant coordinate, the one that
[`modelterms7::term_level_param()`](https://statmodels7.github.io/modelterms7/reference/term_level_param.html)
names. A model carrying a Jeffreys penalty is left on its own route,
whose curvature is factored in two parts.

## See also

[`level_augmented()`](https://statmodels7.github.io/statmodels7/reference/level_augmented.md),
[`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md)
