# Iterated Weighted Least Squares

The scoring step, written out, with the curvature and the decomposition
left to the caller.

## Usage

``` r
iwls(
  hessian = c("expected", "observed"),
  approx = c("bartlett", "integrate", "mc", "opg"),
  decomposition = c("qr", "svd", "chol", "chol_crossprod"),
  maxit = 100L,
  tol = 1e-06,
  step_halving = 30L
)

# S3 method for class 'Iwls'
print(x, ...)
```

## Arguments

- hessian:

  Either `"expected"` – Fisher scoring – or `"observed"`, which is
  Newton.

- approx:

  The approximation of the expected information, read only where the
  family has no closed form.

- decomposition:

  How the step is solved: `"qr"`, `"svd"`, `"chol"` or
  `"chol_crossprod"`.

- maxit:

  The iteration budget.

- tol:

  The stopping tolerance on the scaled score.

- step_halving:

  The number of halvings allowed before a step is abandoned.

- x:

  An
  [`Iwls`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

- ...:

  Unused.

## Value

An object of class
[`Iwls`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md).

## Details

The name is `iwls` and not `irls`: the re-weighting is already implied
by the iteration, so the second `r` says nothing that *iterated* has not
said.

**The curvature.** `hessian = "expected"` inverts the expected
information, which is Fisher scoring; `"observed"` inverts the observed
Hessian, which is Newton. `approx` is handed to distributions7 and is
read only where the family has no closed expected information; asking
for one where it would be ignored is refused, since an argument accepted
and ignored is worse than one that errors.

**The decomposition.** The step solves \\(X'WX + S)\delta = X'g -
S\beta\\, and how depends on this argument.

- `"qr"`:

  a QR of the augmented design \\\[\sqrt{W}X;\\ \mathrm{chol}(S)\]\\.
  The default, and not as a matter of taste: forming \\X'X\\ squares a
  conditioning that the break-point terms bound to \\\epsilon^{-1/2}\\
  by construction, and it is why a closed-form ridge is beaten at \\p =
  200\\ by a method that never forms it.

- `"svd"`:

  a singular value decomposition of the same augmented design, which
  reports the numerical rank rather than failing on a deficient block.

- `"chol"`:

  a Cholesky factor of the penalized information, assembled without the
  augmentation.

- `"chol_crossprod"`:

  the same, forming \\X'WX\\ explicitly. The fastest per iteration and
  the worst conditioned; offered because the choice is the user's.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
iwls()
#> iwls: expected information, qr
#>   maxit 100, tol 1e-06
iwls(hessian = "observed", decomposition = "svd")
#> iwls: observed information, svd
#>   maxit 100, tol 1e-06
```
