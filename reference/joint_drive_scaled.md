# The Scaled Driving Vector's Derivatives in Every Predictor

Computes the entries of
[`joint_drive_data()`](https://statmodels7.github.io/statmodels7/reference/joint_drive_data.md)
for a joint filter whose `scaling` is 1/2 or 1.

## Usage

``` r
joint_drive_scaled(spec, tm, theta, gl, H, D3, D4, D5, order)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- tm:

  The joint filter.

- theta:

  The per-observation parameters at the fitted predictors.

- gl, H, D3, D4, D5:

  The family's link-scale derivatives.

- order:

  As for
  [`joint_drive_data()`](https://statmodels7.github.io/statmodels7/reference/joint_drive_data.md).

## Value

The list that
[`joint_drive_data()`](https://statmodels7.github.io/statmodels7/reference/joint_drive_data.md)
returns.

## Details

With \\u_l = \sum_c S\_{lc}\\ g_c\\ and \\g_c\\ the link-scale gradient,
Leibniz's rule over the positions of a multiset \\T\\ of predictors
gives \$\$U\_{l;T} = \sum\_{A \subseteq T} \sum_c S\_{lc;T \setminus
A}\\ g\_{c;A},\$\$ where \\g\_{c;A}\\ is the family's derivative of
order \\1 + \|A\|\\ with indices \\c\\ and \\A\\, and the derivatives of
\\S\\ come from
[`joint_scale_derivs()`](https://statmodels7.github.io/statmodels7/reference/joint_scale_derivs.md).
The entries are symmetric in \\T\\, so each multiset is computed once
and copied to its permutations.
