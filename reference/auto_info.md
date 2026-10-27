# Information Chosen Against the Family

Selects an information that
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
settles once the family is known. It is the default of the argument
`inner_information` of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
and
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
accept it as well.

In `inner_information`, the expected information is chosen where the
family writes it exactly and at an acceptable cost, and the observed
information otherwise. With
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
a run on the expected information continues on the observed information
as soon as Fisher scoring contracts the score slowly, and a run on the
observed information takes the step on the expected information where
the repaired observed step finds no acceptable point. With an
optimizers7 optimizer the settled matrix is the Hessian that the
optimizer receives. In
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
the expected information is chosen for a Student t, a skew t or a Cauchy
response (wrapped or not), whose observed information is not positive
definite at an outlier, and the observed information is chosen for every
other family. The fitted model records the settled choice.

[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and the variance objects
[`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md),
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
and
[`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md)
do not accept it, since the quantity that each of them computes is
defined on one named matrix.
[information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)
compares the three information objects.

## Usage

``` r
auto_info(method = NULL, repair_floor = 1e-08)
```

## Arguments

- method:

  `NULL` (the default) or an object that inherits from
  [distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html);
  it is used as the `method` of the
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  that the settled choice or the fallback selects:

  - `NULL`: the closed form of the family where it has one, and
    [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
    otherwise;

  - [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html):
    the outer product of the observed scores, which costs one gradient;

  - [`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html):
    quadrature (or summation over the support), with
    `identity = "hessian"` or `identity = "score"`;

  - [`distributions7::expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.html):
    Monte Carlo.

  An object is used also for a family that has the closed form.
  [distributions7::expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.html)
  compares the methods. Any other value signals an error.

- repair_floor:

  A single positive number. A penalized observed information whose
  smallest eigenvalue is below `repair_floor` times the largest absolute
  eigenvalue is read as indefinite by
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md).
  The step on it then floors the eigenvalues (the floor of the solver,
  see
  [`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md)),
  the expected information takes the step only if the search on the
  repaired matrix fails, and a run on it does not take the extra step
  that follows convergence. The default is `1e-8`. The argument is read
  by
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  only.

## Value

An S7 object of class `AutoInfo`, which inherits from
[distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html),
with the properties `method` and `repair_floor`.

## See also

[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
and
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
for the objects that name a matrix;
[information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)
for the comparison;
[distributions7::expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.html)
for the methods;
[AutoInfo](https://statmodels7.github.io/statmodels7/reference/AutoInfo.md)
for the class;
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
for the functions that take the object.

## Examples

``` r
auto_info()
#> Information chosen against the family (closed form where the family has one, otherwise the outer product of the scores)
auto_info(method = distributions7::expect_integrate(identity = "score"),
          repair_floor = 1e-10)
#> Information chosen against the family (expectation by integration against the density (identity = score))

set.seed(1)
dd <- data.frame(x = runif(120, -2, 2))
dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
G <- distributions7::gaussian1_distrib()
fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
               inner_information = auto_info(),
               outer_criterion = reml(information = auto_info()))
fit@methods$smooth
#> iwls: qr, expected information (the observed after 10 iterations)
#>   maxit 1000, tol 1e-06
```
