# The Compiled Inner Fit

`engine_inner_fit()` runs the inner penalized fit of
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
in compiled code (`src/engine_iwls.cpp`) when every component of the
model has a compiled route, and returns `NULL` otherwise, in which case
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
calls
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md).
The compiled fit performs the operations of
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
in the same order, with the families and links read through the
registries of distributions7 and linkfunctions7. `engine_route()`
decides whether a model is covered and assembles the inputs.

## Usage

``` r
engine_inner_fit(
  spec,
  design,
  obj,
  beta,
  hyper,
  method,
  groups,
  frozen,
  verbose,
  damp_on_reject,
  guard,
  idx = NULL
)

engine_route(spec, design, obj, beta, hyper, method)

engine_model(spec, design, coef, hyper, how)

engine_penalty_units(spec, design, units, hyper)

engine_kink_unit(pen, th, index, mp)

engine_matrix(x)
```

## Arguments

- spec, design, obj, hyper:

  As in
  [`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md).

- beta:

  The full starting vector, with the held values written in.

- method:

  A settled
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md).

- groups, frozen:

  As
  [`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
  passes them to
  [`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md).

- verbose:

  Whether to print the iterations.

- damp_on_reject, guard:

  As in
  [`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md).

- idx:

  `NULL` for the whole vector, or the positions of the block that
  [`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
  fits with the others held, on the assembled route that
  [`subset_pieces()`](https://statmodels7.github.io/statmodels7/reference/subset_pieces.md)
  takes.

- coef:

  The coefficients, split by distribution parameter.

- how:

  The solver of the inner fit, `"qr"`, `"svd"` or `"chol"`, which
  decides whether the square root of the penalty is formed.

- units:

  The penalized units,
  [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md).

- pen, th, index, mp:

  A penalized unit's penalty, its hyperparameters, its stacked positions
  and its map.

- x:

  A base matrix or a `Matrix` object.

## Value

`engine_inner_fit()`: the list that
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
returns, or `NULL`. `engine_route()`: `NULL`, or a list with the
elements `model`, `pen` and `method` that the compiled fit reads.

## Details

A model is covered when `spec@engine` is `"auto"`; the method is an
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
without a `criterion` whose expected information is the closed form
(`method = NULL`),
[`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html),
[`distributions7::expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.html)
or
[`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html)
(the same object for the information and its backup); the block is the
whole coefficient vector, with no break-point kink, no moving block and
no nonlinear curvature unit; the family is a univariate family of the
registry without a wrapper and every link has a compiled route; no term
carries a structural state, a Jeffreys penalty or a mixed penalty unit;
and every penalized unit is one that `engine_penalty_units()` describes.

Where the R route would signal an error the compiled fit stops and
`engine_inner_fit()` returns `NULL`, so that
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
reproduces it. With
[`distributions7::expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.html)
the random stream is then restored to its state before the compiled fit,
so that
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
makes the same draws as when it runs alone.

`engine_model()` returns the data of the model as the compiled routes
read them, `list(model, pen, covered)`, or `NULL` where the model is not
covered: the family, its links and the design blocks, and the penalty at
`hyper` (its units, its Hessian and, for the QR and SVD solvers, its
square root). `covered` is the highest derivative order of the family in
the registry.

`engine_penalty_units()` describes each penalized unit to the compiled
fit by the quantities that its `penalty_value()` and
`penalty_gradient()` methods read at the given hyperparameters, so that
the compiled value and gradient repeat those methods operation by
operation: for a `QuadraticPenalty` the map, `P`, `penalty_matrix` and
the three constants; for an `AdditivePenalty` the summed matrix of
`additive_sum()`, the matrices and their smoothing parameters; for a
`DistribPenalty` with a gaussian prior of fixed mean, the mean and the
standard deviation, or, for a multivariate one, the mean, the factor of
[`parameters7::param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.html),
the log-determinant and the inverse covariance. It returns `NULL` for
any other unit.

`engine_kink_unit()` describes a penalty with a kink to the compiled
fit: kind 5 for the lasso (a
[`penalties7::DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.html)
whose `fixed()` parent is a `laplace2` at zero), 6 for the elastic net
(an `enet` parent at zero), 7 for
[`penalties7::scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.html)
and 8 for
[`penalties7::mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.html),
with the diagonal of the map, the rate and the constants of the density,
or the curvature. It returns `NULL` for a penalty without a kink, and
`FALSE` for a kinked one that the compiled fit does not describe (a map
that is not diagonal, a block of more than one coordinate, a location
other than zero).

`engine_matrix()` returns a penalty matrix as the compiled fit reads it:
a base matrix as it is, any `Matrix` object as a `dgCMatrix`, since the
R route multiplies by it as a sparse matrix.

## See also

[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md),
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
(argument `engine`)
