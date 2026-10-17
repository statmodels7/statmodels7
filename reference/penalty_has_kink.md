# Does a Penalty Have a Kink?

Returns `TRUE` when the penalty is not differentiable at some point that
a coefficient can take. A block whose penalty gives `TRUE` is estimated
outside the jointly fitted system by
[`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md):
a compiled coordinate descent where the penalty describes its proximal
operator as a table, and an accelerated proximal gradient iteration
otherwise.

## Usage

``` r
penalty_has_kink(pen, what = "a penalty")
```

## Arguments

- pen:

  A penalties7 penalty.

- what:

  How to name the penalty in the error message.

## Value

A single logical.

## Details

[`penalties7::penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.html)
is read at a probe value of the hyperparameters, which
[`penalty_theta_start()`](https://statmodels7.github.io/statmodels7/reference/penalty_theta_start.md)
gives. The probe is the midpoint of the interval where both bounds are
finite, one unit inside the bound where only one is finite, and zero
where both bounds are infinite. modelterms7 uses the same rule. Whether
a kink exists is a property of the family and not of a point, so any
admissible probe gives the same answer.

If
[`penalties7::penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.html)
signals an error, `penalty_has_kink()` signals an error that names the
penalty, and the penalty is not treated as smooth. If the failure were
read as an answer, a kinked penalty, such as `scad()` or `mcp()`, would
be sent to the scheme for smooth penalties, which fits it by the
curvature of a function that has none.
