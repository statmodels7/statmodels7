# Settle a Marginal Criterion's Information Against the Family

Replaces `hessian = "auto"` by the information the criterion uses for
`distrib`: `"expected"` for a Student t, a skew t or a Cauchy response,
wrapped or not, and `"observed"` for every other family. Any other value
is returned unchanged.

## Usage

``` r
outer_resolve(method, distrib)
```

## Arguments

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- distrib:

  The response family.

## Value

`method`, with `hessian` settled.

## Details

The log-density of these families is not concave in the location, so
their observed information is not positive definite at an outlier, and
the Laplace determinant built on it fails there. Measured on
[`MASS::GAGurine`](https://rdrr.io/pkg/MASS/man/GAGurine.html) with
`GAG ~ Age | sigma ~ Age | nu ~ Age` and a `student_t1` response: on the
observed information the outer search ran 91.5 s and did not converge,
on the expected information it converged in 7.9 s at a larger
log-likelihood (-923.891 against -924.106). Without the equation in `nu`
the two converge to the same point, the expected information costing 8.2
s against 2.6 s (Giovanni, 2026-10-07: the expected information for the
t).
