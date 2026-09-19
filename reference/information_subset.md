# The Information Over Some of the Coefficients

The rows and columns of
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
at the given positions, assembled from those columns of the designs
alone.

## Usage

``` r
information_subset(spec, design, ev, expected, approx, H, index)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, already refreshed at the coefficients.

- ev:

  What
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  returns at the coefficients.

- expected, approx:

  As in
  [`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md).

- H:

  The second-derivative components, or `NULL` to ask for them.

- index:

  Integer positions in the stacked coefficients.

## Value

A symmetric matrix over `index`, in that order: a Matrix object when any
equation's design is sparse, a base matrix otherwise.

## Details

Block \\(a, b)\\ of the result is \\X\_{a,J_a}'\\\mathrm{diag}(w\\
h\_{ab})\\X\_{b,J_b}\\ over the columns \\J_a\\ and \\J_b\\ the
positions select in each equation, which is exactly the corresponding
submatrix of the whole information and costs a cross product over those
columns instead of over all of them. A mixture over regimes sums the
components' blocks as the whole matrix does.

## See also

[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
