# Get geneset collection

Get geneset collection

## Usage

``` r
get_geneset(
  geneset,
  source = c("auto", "builtin", "list", "gmt", "data.frame", "msigdb", "go", "kegg",
    "reactome"),
  species = "Homo sapiens",
  collection = "H",
  subcollection = NULL,
  ontology = c("BP", "MF", "CC"),
  version = NA_character_
)
```

## Arguments

- geneset:

  Geneset input.

- source:

  Geneset source. One of `auto`, `builtin`, `list`, `gmt`, `data.frame`,
  `msigdb`, `go`, `kegg`, `reactome`.

- species:

  Species label used by external geneset sources.

- collection:

  MSigDB collection (for source = `msigdb`).

- subcollection:

  MSigDB subcollection (for source = `msigdb`).

- ontology:

  GO ontology (`BP`, `MF`, `CC`) for source = `go`.

- version:

  Optional source version string.

## Value

Named list of pathways to genes with source metadata attached.
