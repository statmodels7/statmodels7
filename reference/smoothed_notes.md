# The Notes a Smoothed Break-Point Term Adds to a Summary

One note per smoothed term, naming the smoother and the width the build
resolved – the width of the transition, which is part of the model and
not a detail – and, where a break-point carries a random development
under a Gaussian precision and the smoother declares a scale correction,
the corrected scale beside the apparent one.

## Usage

``` r
smoothed_notes(spec, object)
```

## Arguments

- spec:

  The fitted specification, whose terms are the ones the fit left.

- object:

  The fit.

## Value

A character vector, possibly empty.

## Details

The correction is the smoother's own: the probit satisfies the exact
convolution identity \\\tau^2\_{\mathrm{apparent}} = \tau^2 + h^2\\, so
the corrected scale is \\\sqrt{\tau^2 - h^2}\\; a smoother declaring
none (the hyperbolic, the quintic) gets the apparent scale alone, which
the random effect's own block already reports. The apparent scale is
read off the ridge precision as \\1/\sqrt{\lambda}\\, which is only a
scale where the penalty is the quadratic branch with that
hyperparameter; any other development is left without the note rather
than given a number of the wrong meaning.

## See also

[`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.html)
