# Start the Matrix Prior of a Random Effect at the Scale of the Data

Moves the starting free values of every estimated matrix parameter (a
covariance, precision or scale matrix carried by a parameters7 chart) so
that the matrix describes effects whose variance is half the variance of
the equation's working response at the intercept-only fit. The pattern
of the matrix at the neutral start (the correlations, which are zero
there) is kept, and only its scale changes.

## Usage

``` r
statmod_hyper_scale_start(spec, design, hyper)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

- hyper:

  The hyperparameters as
  [`statmod_hyper_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_start.md)
  returns them.

## Value

`hyper`, with the free values of each matrix parameter moved.

## Details

The neutral start of a chart is the identity matrix, a variance of one
on the link scale. On a response measured in thousands that is a prior
under which the effects are negligible, the marginal criterion is flat
in the scale there, and an outer search started there can leave along a
ridge instead of reaching the interior optimum.

The working response of equation \\p\\ is \\z_i = \eta_p + U\_{ip} /
I\_{ip}\\, with \\U\_{ip}\\ the score and \\I\_{ip}\\ the expected
information of observation \\i\\ on the link scale, both at the
intercept-only fit of
[`statmod_intercepts()`](https://statmodels7.github.io/statmodels7/reference/statmod_intercepts.md).
For a gaussian location on the identity link it is the response itself.
Half its weighted variance \\v\\ is given to the effects: the target
matrix is \\c M_0\\ for a covariance or a scale matrix and \\M_0 / c\\
for a precision, with \\M_0\\ the matrix at the neutral start and \\c\\
chosen so that the mean variance of the effects is \\v\\. The free
values are
[`parameters7::param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.html)
of that target.

A matrix parameter is left at its neutral start when any of its free
values is held, when the unit belongs to a structural term or to a
covariance class, when the response is multivariate, and when any step
above fails or returns a value that is not finite.

## See also

[`statmod_hyper_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_start.md),
which supplies the neutral start.
