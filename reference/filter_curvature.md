# A Filter's Second-Order Recursion at the Fitted Point

Runs
[`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html)
at the parameters at which the filter was fitted. The family's
derivatives are looked up instead of evaluated again, and the function
returns the forward Jacobian of the filtered predictor beside the
contracted second derivative.

## Usage

``` r
filter_curvature(spec, design, f, ap, Vs, gl, H, D3, theta, g = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- f:

  The filter's entry of
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)'s
  `filters`.

- ap:

  Which distribution parameter the filter sits in.

- Vs:

  The static rows, one matrix per distribution parameter over the
  coefficients followed by the term's parameters.

- gl, H, D3:

  The family's first three derivatives on the link scale at the fitted
  predictors.

- theta:

  The per-observation parameters at which they were read, from which the
  expected information is read where the filter scales its score.

- g:

  The per-observation weights against which the second derivative of the
  filtered predictor is contracted, or `NULL` (the default) for the
  observation weights times the score of the filtered parameter, the
  term of the observed information. The Jeffreys penalty
  ([`jfp_hessian()`](https://statmodels7.github.io/statmodels7/reference/jfp_hessian.md))
  passes its own.

## Value

The list that
[`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html)
returns.

## Details

Three callers read it at the same point:
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md),
for the joint information;
[`filter_joint_movement()`](https://statmodels7.github.io/statmodels7/reference/filter_joint_movement.md),
for the rows along which the exact outer gradient reads how the mode's
movement reaches the determinant; and
[`filter_joint_jacobian()`](https://statmodels7.github.io/statmodels7/reference/filter_joint_jacobian.md),
for the standard error of a filtered predictor. All go through the one
memo slot, so the recursion runs once at a point.

The recursion is run again here at parameters at which it has already
been run, so the predictor that it reaches is the one at which the
derivatives were read, and the callbacks read the stored values instead
of evaluating the family again. In the filter itself the score is
evaluated at a predictor that the recursion has just produced, which
cannot be known in advance. The seed is part of the memo's key, so a
caller whose layout differed would miss the cache and not take the wrong
shape.

The assembly of the callback is written once, in
[`.structural_blocks()`](https://statmodels7.github.io/statmodels7/reference/dot-structural_blocks.md):
the exact gradient of a penalized filter needs the same pieces at one
order higher. The same pieces are also passed as data, with which an
eligible term runs its second-order recursion compiled, and the callback
is kept for the cases that the kernel does not handle.
