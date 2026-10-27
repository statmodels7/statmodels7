# Rows of the Derivative of a Block

Computes the derivative of order `order` with respect to `var` of every
row of the block of a term at the rows of `grid`, by one central
finite-difference stencil of accuracy four on the rows that
[`modelterms7::term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.html)
returns.

## Usage

``` r
design_derivative(tm, grid, var, order, lo, hi)
```

## Arguments

- tm:

  A built
  [modelterms7::SmoothTerm](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.html).

- grid:

  A data frame with every variable the term reads.

- var:

  The variable that is differentiated.

- order:

  The order of the derivative.

- lo, hi:

  The range of the variable.

## Value

A matrix with one row for each row of `grid` and one column for each
coefficient of the term.

## Details

The step is \\h = 0.05\\\mathrm{fd\\step}(1, m, 4)\\(b - a)\\, where
\\m\\ is the order and \\\[a, b\]\\ the range of the variable. The
center of the stencil is moved inside \\\[a + rh, b - rh\]\\, where
\\r\\ is the reach of the stencil, and the weights of the shifted nodes
are computed for the point itself, so that the nodes stay inside the
range over which the basis is defined.
