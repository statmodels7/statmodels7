# The Joint Second Derivative of the Penalized Information, Joint Filter

Returns, for a model carrying a joint filter, the quantity that
[`structural_chain_extra2()`](https://statmodels7.github.io/statmodels7/reference/structural_chain_extra2.md)
returns for one filter: \\\mathrm{tr}(M\\\partial^2 K/\partial u^2\[v,
w\])\\ less its design part.

## Usage

``` r
joint_chain_extra2(spec, design, jd, M, st, v, w)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- jd:

  The joint rows.

- M:

  The matrix against which the trace is taken.

- st:

  The shared quantities, from
  [`joint_grad_parts()`](https://statmodels7.github.io/statmodels7/reference/joint_grad_parts.md)
  at order 4.

- v, w:

  The two directions, over the estimated coordinates.

## Value

A single number.

## Details

The terms are those of the single-filter function, with every quantity
that reads the filter's equation summed over the filtered parameters.
The second derivative of the score of parameter \\k\\ along the two
directions is \$\$\sum\_{a,b}\ell\_{kab}(V_a\cdot v)(V_b\cdot w) +
\sum_j \ell\_{kj}\\ v^\top\Phi_j w,\$\$ in which only the filtered
equations have a second derivative \\\Phi_j\\. The movements of the
filtered rows along the two directions pair every two filtered
parameters through \\\ell\_{jj'}\\.
