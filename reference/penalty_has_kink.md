# Does a Penalty Have a Kink?

`TRUE` when the penalty is not differentiable somewhere a coefficient
can be. A block whose penalty answers `TRUE` is estimated outside the
jointly fitted system, by a coordinate descent of its own.

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

[`penalties7::penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.html)
is read at a probe value of the hyperparameters, the midpoint of their
bounds, which is the rule modelterms7 already uses. Whether a kink
exists is a property of the family, not of a point, so any admissible
probe answers.

A penalty that stops when asked is reported, never treated as smooth.
Reading the failure as an answer sends the term to the scheme for the
opposite property: `scad()` and `mcp()` were fitted by the curvature of
a function that has none, reporting an effective 19.00 degrees of
freedom out of 20 on a design of pure noise, which is no selection at
all.
