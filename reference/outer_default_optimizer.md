# Which Optimizer the Outer Search Uses When the Caller Names None

The choice is made from what the criterion can supply: its exact
Hessian, its exact gradient, or neither.

## Usage

``` r
outer_default_optimizer(exact, exact2)
```

## Arguments

- exact:

  Whether the criterion has an exact gradient.

- exact2:

  Whether it has an exact Hessian as well.

## Value

An optimizers7 optimizer.

## See also

[`outer_fit`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md),
[`outer_gradient_ok`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
