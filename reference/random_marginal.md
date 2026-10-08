# A Marginal Prediction: the Population Average Over New Groups

Integrates the quantity asked for over the prior of every random-effect
term read at `"marginal"`, at the fixed part of each equation's
predictor.

## Usage

``` r
random_marginal(
  spec,
  design,
  ep,
  mt,
  nodes,
  what,
  eta_at = NULL,
  deriv = FALSE,
  eta_node = NULL
)
```

## Arguments

- spec, design:

  The specification at the new rows and its design, with the marginal
  terms' columns set to zero.

- ep:

  The predictors at the fixed part, as
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  returns them.

- mt:

  The marginal terms, rows of
  [`random_modes()`](https://statmodels7.github.io/statmodels7/reference/random_modes.md).

- nodes:

  What
  [`random_nodes()`](https://statmodels7.github.io/statmodels7/reference/random_nodes.md)
  returns, in the order of `mt`.

- what:

  What was asked for.

- eta_at:

  A named list of predictors to integrate at in place of `ep$eta`, for
  the ends of an interval; `NULL` for `ep$eta`.

- deriv:

  `TRUE` to average the derivative of the inverse link instead of the
  inverse link, which is what the delta method of a marginal parameter
  multiplies the predictor's standard error by. Read only for a
  parameter.

- eta_node:

  For random effects written inside a subformula, a function of the node
  index returning the predictors at that node, where the term holding
  them is evaluated at the node's values. It takes the place of `eta0`,
  and the effects written in an equation are added to it as usual. The
  check of a predictor's domain reads `eta0` and is not run.

## Value

A numeric vector, or a named list of them for `"parameter"`.

## Details

At each node the predictor of equation \\p\\ is \\\eta\_{0p} + \sum_t
Z_t b_t\\, with \\Z_t\\ the term's within-group design
([`modelterms7::term_within()`](https://statmodels7.github.io/modelterms7/reference/term_within.html)).
A parameter is averaged, \\E_b\[h^{-1}(\eta_0 + Zb)\]\\; the mean is
averaged too, by the law of total expectation; the variance adds the
variance of the conditional mean, \\E_b\[\mathrm{Var}(Y\mid b)\] +
\mathrm{Var}\_b(E\[Y\mid b\])\\, which is where a random effect on the
scale enters. A predictor, `"link"` or `"link:p"`, is \\\eta_0\\ itself,
the effects having mean zero.
