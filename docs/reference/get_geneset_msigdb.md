# Get MSigDB genesets via msigdbr

Get MSigDB genesets via msigdbr

## Usage

``` r
get_geneset_msigdb(
  species = "Homo sapiens",
  collection = "H",
  subcollection = NULL
)
```

## Arguments

- species:

  Species label accepted by `msigdbr`.

- collection:

  MSigDB collection, e.g. `H`, `C2`, `C5`.

- subcollection:

  Optional subcollection, e.g. `CP:KEGG`.

## Value

Named list genesets.
