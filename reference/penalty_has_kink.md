# Does a Penalty Have a Kink?

`TRUE` when the penalty is not differentiable somewhere a coefficient
can be, which is what puts its block outside the jointly fitted system.

## Usage

``` r
penalty_has_kink(pen)
```

## Arguments

- pen:

  A penalties7 penalty.

## Value

A single logical.

## Details

[`penalty_kinks`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.html)
is read at a probe value of the hyperparameters – the midpoint of their
bounds, the rule modelterms7 already uses – because whether a kink
exists is a property of the family and not of a point.
