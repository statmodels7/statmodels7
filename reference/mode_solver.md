# The Movement of the Mode Under Held Coordinates and Kinks

Builds the solve that gives the movement of the penalized mode under a
perturbation \\c\\: \\v = -N (N^\top (K + D_m) N)^{-1} N^\top c\\, with
\\N\\ a basis of the directions that the held coordinates and the kinks
leave free.

## Usage

``` r
mode_solver(spec, design, coef, K, Kinv, Dm, npar, offs)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design at the mode.

- coef:

  The coefficients at the mode, a named list.

- K:

  The mode's penalized information, as
  [`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md)
  returns it with `laplace = FALSE`.

- Kinv:

  Its inverse.

- Dm:

  The block's own curvature from
  [`mode_curvature()`](https://statmodels7.github.io/statmodels7/reference/mode_curvature.md),
  or `NULL`.

- npar, offs:

  The coefficient counts and offsets per parameter.

## Value

A list with `solve`, a function of \\c\\ returning \\v\\ without the
sign; `J`, the curvature solved with; `Dm`, the curvature with the held
rows and columns set to zero; and `hold`, the held positions.

## Details

The held coordinates are those pinned in \\K\\ other than at a boundary
(a coefficient that the criterion estimates, a direction that the model
does not identify, a coefficient that a kinked penalty holds at zero).
The kinks are the break-points that sit on an observation at the mode
([`kink_constraints()`](https://statmodels7.github.io/statmodels7/reference/kink_constraints.md)):
there the one-sided slopes have opposite signs, so a small perturbation
leaves the break-point on the observation, as a lasso leaves a
coefficient at zero. Where the break-point's position is a single
coefficient that coefficient does not move; where it is developed over a
design the mode moves along the hyperplane on which the break-point
stays on the observation. The solve is the bordered system of the
equality constrained problem, written in the basis of its null space.
The criterion's own determinant is not affected: it keeps the side
convention of the break-point term.

Where \\N^\top (K + D_m) N\\ is not positive definite, \\K\\ alone
replaces \\K + D_m\\.

## See also

[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md),
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md),
[`null_basis()`](https://statmodels7.github.io/statmodels7/reference/null_basis.md)
