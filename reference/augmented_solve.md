# Solve a Scoring Step From the Square-Root Design

Decomposes the augmented matrix \\\[R;\\ C\]\\, whose cross-product is
the penalized information, and returns the increment solving \\(R'R +
C'C)\delta = u\\.

## Usage

``` r
augmented_solve(R, C, u, how, threads = 1L)
```

## Arguments

- R:

  The square-root design.

- C:

  The penalty's factor.

- u:

  The right-hand side.

- how:

  Either `"qr"` or `"svd"`.

- threads:

  How many threads the triangular factor may use.

## Value

A list with `delta` and `rank`.
