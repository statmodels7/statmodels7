# Which Coefficients Belong to a Block That Is Not a Jacobian

The positions of every coefficient of a term whose design block is a
working linearization with a frozen weight.

## Usage

``` r
frozen_block(spec, lab)
```

## Arguments

- spec:

  The fitted specification.

- lab:

  The coefficient labels.

## Value

A logical vector over the coefficients.

## Details

A block whose term answers `FALSE` to
[`modelterms7::term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.html)
is a working linearization with a weight held at the previous iterate,
so the curvature that the block carries is that of the working model and
not that of the model. The standard errors computed from it would be
zero, which would be reported as a measurement, and those coefficients
are left missing.

The check reads the term through
[`modelterms7::term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.html)
and never through its class, so a construction whose block is a Jacobian
keeps its inference: a continuous
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
and a discontinuous one smoothed by an
[`numericals7::abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.html),
both answer `TRUE`.
