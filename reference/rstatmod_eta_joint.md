# Simulate a Model Carrying a Joint Filter

Runs the joint recursion forward drawing the response at each step, the
branch of
[`rstatmod_eta()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_eta.md)
for a joint filter.

## Usage

``` r
rstatmod_eta_joint(spec, tm, eta, theta, psi)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- tm:

  The joint filter.

- eta, theta:

  The static predictors and parameters, named lists.

- psi:

  The filter's parameters on the parameter scale.

## Value

The list that
[`rstatmod_eta()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_eta.md)
returns.

## Details

At each row the response is drawn at the filtered predictors that the
recursion has reached there and the static predictors of the other
equations, and the driving vector is read at the drawn value, the same
vector, scaled as in the fit, that drove the fitted filter.
