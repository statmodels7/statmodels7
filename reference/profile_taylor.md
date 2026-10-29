# Profile by Local Expansion

Selects a method of
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md)
that does not refit. The coefficients that are not held are developed as
a polynomial in the held ones around the estimate, and the
log-likelihood and the penalty are computed exactly at the coefficients
that the polynomial gives.

## Usage

``` r
ProfileTaylor(name = character(0), order = integer(0))

profile_taylor(order = 2)
```

## Arguments

- name:

  The short name of the method, a single string. The class constructor
  `ProfileTaylor()` takes it, and `profile_taylor()` sets it to
  `"taylor"`.

- order:

  The highest order of the derivatives of the objective used, an integer
  from 2 to 5.

## Value

A `ProfileTaylor` object, inheriting from
[profile_method](https://statmodels7.github.io/statmodels7/reference/profile_method.md),
with the properties `name` and `order`.

## Details

With \\\theta\\ the held coefficients and \\\eta\\ the others, the
refitted \\\hat\eta(\theta)\\ solves \\\nabla\_\eta F(\eta, \theta) =
0\\ for the objective \\F\\ of the fit. Along the segment from the
estimate to a point of the grid, the derivatives of \\\hat\eta\\ are
obtained one order after another: the first is \\-H\_{\eta\eta}^{-1}
H\_{\eta\theta}\\\delta\\, and each following one is minus
\\H\_{\eta\eta}^{-1}\\ times the derivative of that order, along the
segment, of the gradient evaluated at the polynomial of the lower
orders. That derivative is taken by one central stencil of numericals7
on the analytic gradient. The polynomial has degree `order - 1`, so that
`order` is the highest order of the derivatives of the objective that
the expansion uses. The Hessian is the observed information of the fit,
with, for a model with nonlinear or smoothed break-point terms, the term
of the block that moves with its coefficients. With the log-likelihood
and the penalty computed exactly at the predicted coefficients, the
penalized log-likelihood of every row is a lower bound of the profile.
The column `converged` is `NA`.

Quantities that the outer criterion estimates (coordinates, and
hyperparameters with `hyper = "refit"`) are developed as in
[`profile_continuation()`](https://statmodels7.github.io/statmodels7/reference/profile_continuation.md),
from the derivatives of the criterion in those quantities, and the
coefficients that are neither held nor estimated by the criterion are
then the penalized mode at the predicted values (one inner fit for each
point). The lower bound does not hold in that case, because the exact
method maximizes the criterion and not the penalized log-likelihood over
those quantities. A hyperparameter in `which` is held at the values of
its grid on the scale of its link, and the column `outer_criterion` is
the profile.

The requirements on the model are those of
[`profile_continuation()`](https://statmodels7.github.io/statmodels7/reference/profile_continuation.md).

## See also

[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md),
[`profile_grid()`](https://statmodels7.github.io/statmodels7/reference/profile_grid.md),
[`profile_continuation()`](https://statmodels7.github.io/statmodels7/reference/profile_continuation.md).

## Examples

``` r
profile_taylor()
#> <statmodels7::ProfileTaylor>
#>  @ name : chr "taylor"
#>  @ order: int 2
profile_taylor(order = 4)
#> <statmodels7::ProfileTaylor>
#>  @ name : chr "taylor"
#>  @ order: int 4
```
