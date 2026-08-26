# Which Hyperparameters Set the Size of the Kink

The names whose value moves
[`kink_scale()`](https://statmodels7.github.io/statmodels7/reference/kink_scale.md),
which are the ones a path over the penalty has to vary.

## Usage

``` r
kink_hypers(pen, theta, unbounded = TRUE)
```

## Arguments

- pen:

  A penalties7 penalty.

- theta:

  The hyperparameters in force.

- unbounded:

  Whether to keep only those with an infinite upper bound.

## Value

A character vector, possibly empty.

## Details

The question is put to the penalty and never answered from a list of
families: the shape parameters of SCAD and MCP leave the subdifferential
at zero unchanged and govern how fast the penalty flattens further out,
while \\\lambda\\ and the elastic net's \\\alpha\\ both scale it.

`unbounded` restricts the answer to the hyperparameters with no upper
bound, which is the default choice of what to select. The reference
implementations do the same by convention, glmnet holding \\\alpha\\
fixed and ncvreg holding \\\gamma\\, and a bounded shape is better swept
by hand over a few values than searched.

## See also

[`kink_scale()`](https://statmodels7.github.io/statmodels7/reference/kink_scale.md)
