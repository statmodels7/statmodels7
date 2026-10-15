# The Blocks of One Distribution Parameter

Groups the terms of a distribution parameter into the blocks that a
summary prints: the parametric terms together, and one block for each
other term (smooth, random effect, selection, other penalized term,
break-point term or structural term).

## Usage

``` r
summary_blocks(
  fit,
  spec,
  design,
  p,
  ci,
  level = 0.95,
  V = NULL,
  st = NULL,
  at_edge = character(0),
  test = "wald"
)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec:

  The specification.

- design:

  The design.

- p:

  The distribution parameter.

- ci:

  The flat interval table, as
  [`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
  returns it with the statistic and the p-value added.

- level:

  The confidence level the intervals are built at.

- V:

  The variance matrix over the stacked coefficients, or `NULL`. It is
  needed only by a term reported through
  [`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html),
  whose quantities are functions of several coefficients at once and
  whose standard errors therefore come from the delta method and not
  from a diagonal entry.

- st:

  The structural table, or `NULL`. A structural term has no design
  columns, so its block is built from the parameters that the term
  reports and not from a block of the design, and its hyperparameter is
  reported there instead of in a block of its own that would carry
  nothing else.

- at_edge:

  The `boundary_key` of
  [`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md):
  the hyperparameter coordinates that have run to an edge of their range
  with their own gradient already met. Such a coordinate reports no
  standard error and no interval, and neither does any quantity that
  depends on it; the variance matrix around it is computed and reported
  as it always was.

- test:

  Which test the statistic column reports. Under a test other than
  `"wald"`, a quantity of a break-point term that is one coefficient
  carries the restricted statistic of that coefficient, and any other
  quantity carries none.

## Value

A list of block records, each with `kind`, `label`, `term`, `n_coef`,
`edf`, `n_zero`, `table` and `components`, and, where they apply,
`classes` and `note` (the covariance blocks and the note that locates
them), `n_na` (the aliased coefficients of a selection) and `spread`
(the summary of the predicted effects of a random effect). A term
written in parameters of its own that develops one of them over
covariates reports that parameter as a compartment of its own, carrying
its hyperparameter and the rows of its sub-terms, and `table` keeps only
what is left. Such a term with a kinked sub-term also carries
`n_kinked`, the number of coordinates under the kink, which `n_zero` and
the count of selected coefficients read. The element `head` is always
`NULL`.
