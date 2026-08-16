# Reject an Offset Buried Inside a Term

Stops where [`offset()`](https://rdrr.io/r/stats/offset.html) appears
anywhere other than as a term of the equation itself.

## Usage

``` r
reject_nested_offsets(eq, param)
```

## Arguments

- eq:

  A one-sided formula, with the equation's own offsets already taken
  out.

- param:

  The parameter the equation belongs to, for the message.

## Value

`NULL`, invisibly; an error where one is found.

## Details

[`split_offsets`](https://statmodels7.github.io/statmodels7/reference/split_offsets.md)
takes the top-level additive terms, which is where R recognizes an
offset. One written inside another term's formula –
`ridge(~ z + offset(o))`, `random(~ 1 + offset(o) | g)`,
`nl(a ~ 0 + ridge(~ g + offset(o)))` – reaches that term's own
`model.matrix`, where it is dropped exactly as it used to be dropped at
the equation level: measured, the fit ran, the block had the columns of
the model WITHOUT it, and the intercept came back at 1.33 against the
-5.01 the offset gives, a factor of 566.

It is REFUSED rather than routed up, and the reason is that the meaning
differs by where it sits. In a penalized term's formula an offset would
be a contribution to the equation's predictor, which is what writing it
at the equation level already says. In a SUBFORMULA it would be a
contribution to that parameter's own chart, which is a different
quantity. Rather than pick one reading and apply it everywhere, the
construct is rejected with the place it belongs named, which costs a
user one edit and cannot silently fit the wrong model.

## See also

[`split_offsets`](https://statmodels7.github.io/statmodels7/reference/split_offsets.md)
