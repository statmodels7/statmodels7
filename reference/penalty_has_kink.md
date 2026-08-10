# Does a Penalty Have a Kink?

`TRUE` when the penalty is not differentiable somewhere a coefficient
can be, which is what puts its block outside the jointly fitted system.

## Usage

``` r
penalty_has_kink(pen, what = "a penalty")
```

## Arguments

- pen:

  A penalties7 penalty.

- what:

  How to name the penalty if it cannot answer.

## Value

A single logical.

## Details

[`penalty_kinks`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.html)
is read at a probe value of the hyperparameters – the midpoint of their
bounds, the rule modelterms7 already uses – because whether a kink
exists is a property of the family and not of a point.

A penalty that stops when asked is reported rather than treated as
smooth. Reading the failure as an answer sends the term to the scheme
for the opposite property: `scad()` and `mcp()` were fitted by the
curvature of a function that has none, reporting an effective 19.00
degrees of freedom out of 20 on a design of pure noise, which is no
selection at all.
