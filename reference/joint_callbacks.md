# The Driving Vector of a Joint Filter, One Observation at a Time

Returns the callbacks that
[`modelterms7::term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.html)
reads for a
[`modelterms7::GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.html):
the driving vector \\u =
\mathcal{I}\_{GG}^{-d}\\\partial\ell/\partial\eta_G\\ at the filtered
predictors of one row, the other predictors held at their static values,
and its derivative in the filtered predictors.

## Usage

``` r
joint_callbacks(spec, tm, eta)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- tm:

  The joint filter.

- eta:

  The static predictors, a named list over the distribution parameters.

## Value

A list with `score` and `curvature`, functions `(e, i)`, `score_y` and
`curvature_y`, the same at a response `y` given as their first argument,
and `fast`, the compiled route's context or `NULL`.

## Details

With \\g\\ the link-scale gradient, \\H\\ the Hessian and \\S =
\mathcal{I}\_{GG}^{-d}\\, \$\$\frac{\partial u_l}{\partial \eta_k} =
\sum_c S\_{lc} H\_{ck} + \sum_c \frac{\partial
S\_{lc}}{\partial\eta_k}\\ g_c .\$\$ The two callbacks are called in
turn at the same point, so the family is evaluated once per point.

Where the response is a vector and the C registries of distributions7
and linkfunctions7 know the family and every link, the result also
carries `fast`, the context with which
[`modelterms7::term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.html)
runs the recursion compiled and reads the family through the joint
registry; a family that the joint registry does not cover leaves that
context inert, and the callbacks run.
