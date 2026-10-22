# The Scaling Matrix of a Joint Filter and Its First Derivative

Returns \\S = \mathcal{I}^{-d}\\ for a symmetric positive definite
block, with \\d\\ in 0, 1/2 or 1, and its derivative along each of a set
of derivatives of the block.

## Usage

``` r
joint_scale_matrix(I, dI, d)
```

## Arguments

- I:

  The block.

- dI:

  A list of its derivatives.

- d:

  The power.

## Value

A list with `S` and `dS`, a list like `dI`.

## Details

For \\d = 1\\, \\\dot S = -S\dot{\mathcal{I}}S\\. For \\d = 1/2\\, \\X =
\mathcal{I}^{1/2}\\ is the symmetric root, whose derivative solves the
Sylvester equation \\X\dot X + \dot X X = \dot{\mathcal{I}}\\; in the
eigenbasis of \\\mathcal{I}\\ that is a division by \\\sqrt{\lambda_i} +
\sqrt{\lambda_j}\\, and then \\\dot S = -S\dot X S\\.
