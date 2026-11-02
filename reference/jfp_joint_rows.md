# The Derivative of Every Predictor Beside a Filter

Returns, for a model carrying a filter, the derivative of each
equation's predictor in the joint vector of the coefficients and the
filter's free parameters, with what the second derivative of a filtered
predictor needs.

## Usage

``` r
jfp_joint_rows(spec, design, coef, ev, extra = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design as it arrived, carrying the structural state.

- coef:

  The coefficients.

- ev:

  What
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  returns at `coef`.

- extra:

  `NULL`, or one matrix per distribution parameter with the static rows
  of further unknowns, placed after the coefficients and before the
  filter's parameters; every coordinate is then kept.

## Value

A list with `V` (one dense matrix per distribution parameter over the
kept coordinates) and `Vfull` (the same over every coordinate, held
levels included), `keep`, `driven` (the equations whose predictor the
filter writes), `curv`, a function of the per-observation weights (a
named list keyed by parameter) returning \\\sum_t w_t\sum_c
g\_{c,t}\\\partial^2\eta\_{c,t}/\partial u\\ \partial u\\ over the kept
coordinates, `curv_full` the same over every coordinate, `third`, a
function of the weights and of a direction \\d\\ over every coordinate
returning a list with `curv`, the third derivative contracted against
\\w_tg\_{c,t}\\ and \\d\\, and `dphi`, the second derivative along
\\d\\, one matrix per parameter, both over every coordinate
([`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html)),
`higher`, the same for one, two or three directions (a list) returning
`curv`, the derivative of order two plus their number contracted against
the weights and the directions, and `rows`, the derivative of order one
plus their number along them
([`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html),
[`modelterms7::term_fifth()`](https://statmodels7.github.io/modelterms7/reference/term_fifth.html)),
and `nb`, `mf`, `zn` and `term`.

## Details

The rows are those of
[`joint_design_rows()`](https://statmodels7.github.io/statmodels7/reference/joint_design_rows.md),
with the equation that the filter drives replaced by the forward
Jacobian of the recursion
([`filter_curvature()`](https://statmodels7.github.io/statmodels7/reference/filter_curvature.md),
or
[`joint_filter_curvature()`](https://statmodels7.github.io/statmodels7/reference/joint_filter_curvature.md)
for a joint filter), and a level that an intercept holds dropped as
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
drops it.
