# Does a Penalty Supply What a Marginal Criterion Needs?

Asks the penalty for the derivatives, at a probe value of its
hyperparameters, and reports whether it answered.

## Usage

``` r
penalty_answers(pen, order = 1L)
```

## Arguments

- pen:

  A penalties7 penalty.

- order:

  `1` or `2`.

## Value

A single logical.

## Details

The order-2 route additionally requires the penalty to be quadratic in
the coefficients
([`beta_quadratic`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.html)),
since otherwise the third and fourth derivatives of the penalty in
\\\beta\\ enter the criterion's own second derivative and penalties7
does not carry them.
