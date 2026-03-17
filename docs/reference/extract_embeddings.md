# Extract embeddings from generic input

Extract embeddings from generic input

## Usage

``` r
extract_embeddings(
  object = NULL,
  embeddings = NULL,
  reduction = "umap",
  seurat = TRUE
)
```

## Arguments

- object:

  Seurat object.

- embeddings:

  Numeric matrix, optional direct input.

- reduction:

  Reduction name when extracting from Seurat.

- seurat:

  Input mode.

## Value

Embedding matrix.
