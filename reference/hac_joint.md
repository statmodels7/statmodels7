# The HAC Variance Over the Joint Unknowns

Computes \\V_b B V_b\\ over every unknown of the joint information, for
the table of a structural term's own parameters.

## Usage

``` r
hac_joint(fit, design, Vb, sw)
```

## Arguments

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- design:

  The design of the fit.

- Vb:

  The inverse penalized information over the unknowns of the joint
  information.

- sw:

  The
  [HacSandwich](https://statmodels7.github.io/statmodels7/reference/HacSandwich.md)
  object.

## Value

A matrix of the dimension of `Vb`.
