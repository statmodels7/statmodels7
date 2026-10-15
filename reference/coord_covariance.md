# Which Way of Holding the Gradient Is Cheaper

Chooses how the compiled sweeps keep the gradient: `TRUE` for the
covariance form, which holds the gradient itself and caches columns of
\\X'WX\\, and `FALSE` for the running residual. The test is
`m <= 32 && n > 8 * m`.

## Usage

``` r
coord_covariance(n, m)
```

## Arguments

- n:

  The number of observations.

- m:

  How many coordinates the strong rule left to visit.

## Value

A single logical.

## Details

The covariance form replaces an \\O(n)\\ read of the gradient with an
\\O(m)\\ one, and pays for it by building a column of \\X'WX\\ at
\\O(nm)\\ the first time the coordinate is nonzero. It is used only
while \\m\\ is small next to \\n\\: \\m \le 32\\ and \\n \> 8m\\. When
many coordinates are live, the Gram columns cost more than the residual
passes that they replace, so the residual form is used.

## See also

[`coord_call()`](https://statmodels7.github.io/statmodels7/reference/coord_call.md),
which passes the answer to the kernel,
[`coord_screen()`](https://statmodels7.github.io/statmodels7/reference/coord_screen.md)
for `m`.
