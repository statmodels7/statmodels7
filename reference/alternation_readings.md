# Whether the Alternation Has Reached Its Point

Reads, at the point the alternation reached, how far the smooth part is
from its mode and how far each coefficient a kinked penalty set to zero
is from its optimality condition, both in log-likelihood units, and
compares them with
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md).

## Usage

``` r
alternation_readings(
  spec,
  design,
  obj,
  beta,
  hyper,
  expected,
  approx,
  aliased = integer(0)
)

alternation_settled(
  spec,
  design,
  obj,
  beta,
  hyper,
  expected,
  approx,
  aliased = integer(0)
)
```

## Arguments

- spec, design, obj, beta, hyper:

  The fit's specification, design, objective, stacked coefficients and
  hyperparameters.

- expected, approx:

  Which information the fit uses.

- aliased:

  The positions the scoring step's pivot left out.

## Value

`alternation_readings()` a named vector of `mode` and `kkt`, the second
`NA` where no coefficient sits at a kink; `alternation_settled()` a
single logical.

## Details

The smooth part is every coordinate that is not a zero of a kinked
penalty, not in a frozen working block, not aliased, not held and not
named by
[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)
on the information over the others, together with a structural term's
own free parameters where the model carries one. Its reading is the
Newton decrement \\\tfrac12 g^\top K^{-1} g\\, with \\g\\ the gradient
of the penalized objective and \\K = H + S\\ the penalized information
over those coordinates (the joint one, through
[`statmod_joint_pieces()`](https://statmodels7.github.io/statmodels7/reference/statmod_joint_pieces.md),
beside a structural term). A coefficient a kinked penalty holds at zero
is optimal where the pull of the likelihood, \\s_j\\, does not exceed
the size \\\kappa_j\\ of the kink, and its reading is the KKT decrement
\$\$\frac{\max(0,\\ \lvert s_j\rvert - \kappa_j)^2}{2\\c_j},\$\$ the
gain a Newton step on that coordinate alone would predict, with \\c_j\\
the diagonal of the unpenalized information. The largest such decrement
is reported.

A reading that cannot be computed is `NA`, and `alternation_settled()`
reads `NA` on the smooth part as not settled.

## See also

[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md),
[`zero_readings()`](https://statmodels7.github.io/statmodels7/reference/zero_readings.md),
[`inner_mode_error()`](https://statmodels7.github.io/statmodels7/reference/inner_mode_error.md)
