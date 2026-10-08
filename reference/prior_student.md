# The Scale Matrix and Degrees of Freedom of a Student t Prior

The Scale Matrix and Degrees of Freedom of a Student t Prior

## Usage

``` r
prior_student(pen, th, d)
```

## Arguments

- pen:

  A penalty.

- th:

  Its hyperparameters, a named list.

- d:

  The number of coefficients a group carries.

## Value

A list with `Sigma` and `nu` where the penalty is a Student t prior with
one mixing variable a group (a multivariate t, or a univariate one on a
single coefficient), and `NULL` otherwise.

## Examples

``` r
tp <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
pen <- penalties7::distrib_penalty(tp, n_coef = 4)
statmodels7:::prior_student(pen, list(sigma = 2, nu = 3), 1)
#> $Sigma
#>      [,1]
#> [1,]    4
#> 
#> $nu
#> [1] 3
#> 
```
