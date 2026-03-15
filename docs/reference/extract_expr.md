# Extract expression matrix from input

Extract expression matrix from input

## Usage

``` r
extract_expr(
  object = NULL,
  expr = NULL,
  assay = NULL,
  layer = NULL,
  slot = NULL,
  seurat = TRUE
)
```

## Arguments

- object:

  Seurat object.

- expr:

  Matrix input when `seurat = FALSE`.

- assay:

  Assay name.

- layer:

  Layer name.

- slot:

  Legacy slot.

- seurat:

  Input mode.

## Value

Matrix or `dgCMatrix`.
