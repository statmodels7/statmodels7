# The Notes That a Smoothed Break-Point Term Adds to a Summary

Returns a note for each smoothed break-point term that names the
smoother and the width that was resolved at build, and a second note for
each random development of a break-point under a Gaussian prior, when
the smoother declares a scale correction, that gives the corrected scale
beside the apparent one.

## Usage

``` r
smoothed_notes(spec, object)
```

## Arguments

- spec:

  The fitted specification, whose terms are those that the fit left.

- object:

  The fit.

## Value

A character vector, possibly empty.

## Details

The correction is declared by the smoother. The probit satisfies the
exact convolution identity \\\tau^2\_{\mathrm{apparent}} = \tau^2 +
h^2\\, so the corrected scale is \\\sqrt{\tau^2 - h^2}\\. A smoother
that declares no correction (the hyperbolic and the quintic) receives
only the first note. The apparent scale is the standard deviation
parameter of a Gaussian prior written by its scale, or
\\1/\sqrt{\lambda}\\ for a quadratic penalty with the single
hyperparameter \\\lambda\\. A development with any other prior receives
no scale note.

## See also

[`numericals7::abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.html)
