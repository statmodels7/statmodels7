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

A discontinuous break-point term is fitted through a block whose weight
is held at the previous iterate, so the curvature that block carries is
the working model's and not the model's. Measured on a jump at 400
observations against a bootstrap of 200 resamples: the working
information gives the change of level and the auxiliary coordinate a
standard error of EXACTLY ZERO, where the resamples give 0.063 and
0.540, and the position read off them 1.8e-05 against 0.090. A zero
looks like a number and is worse than a gap, so those coefficients are
left missing.

The question is asked of the term through
[`term_jacobian_block`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.html)
rather than of its class, so a construction whose block IS a Jacobian
keeps its inference: a continuous
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.html),
and a discontinuous one smoothed by an
[`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.html),
both answer yes.
