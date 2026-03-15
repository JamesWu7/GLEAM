# Canonical spatial signature testing interface

Canonical spatial signature testing interface

## Usage

``` r
test_signature_spatial(
  score,
  region,
  group = NULL,
  sample = NULL,
  method = c("wilcox", "t"),
  adjust_method = "BH",
  level = c("region", "sample_region")
)
```

## Arguments

- score:

  `gleam_score` object.

- region:

  Region/domain metadata column or vector.

- group:

  Optional group for within-region contrasts.

- sample:

  Optional sample variable.

- method:

  Statistical method.

- adjust_method:

  Multiple testing adjustment.

- level:

  Region level (`region` or `sample_region`).

## Value

`gleam_test` object.
