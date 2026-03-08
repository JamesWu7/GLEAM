# GLEAM supported geneset sources

## Official built-in species focus

GLEAM built-in geneset workflows officially target:

- `human`
- `mouse`

Unsupported species for built-in sources fail clearly and recommend
custom genesets.

## Source overview

``` r

list_geneset_sources()
#>       source requires_package supported_builtin_species built_in_available
#> 1    builtin             <NA>               human,mouse               TRUE
#> 2       list             <NA>               any(custom)              FALSE
#> 3        gmt             <NA>               any(custom)              FALSE
#> 4 data.frame             <NA>               any(custom)              FALSE
#> 5     msigdb          msigdbr               human,mouse              FALSE
#> 6         go          msigdbr               human,mouse              FALSE
#> 7       kegg          msigdbr               human,mouse              FALSE
#> 8   reactome          msigdbr               human,mouse              FALSE
#>   internet_required
#> 1             FALSE
#> 2             FALSE
#> 3             FALSE
#> 4             FALSE
#> 5             FALSE
#> 6             FALSE
#> 7             FALSE
#> 8             FALSE
```

## Source examples

``` r

# Built-in hallmark-like set
gs_builtin <- get_geneset("hallmark", source = "builtin", species = "human")

# User-defined list
gs_list <- get_geneset(list(PATH_A = c("G1", "G2")), source = "list")

# GMT file
gs_gmt <- get_geneset("path/to/geneset.gmt", source = "gmt")

# data.frame input
gs_df <- get_geneset(data.frame(pathway = c("A", "A"), gene = c("G1", "G2")), source = "data.frame")

# Optional curated sources (msigdbr backend)
gs_go <- get_geneset(NULL, source = "go", species = "human", ontology = "BP")
gs_kegg <- get_geneset(NULL, source = "kegg", species = "human")
gs_reactome <- get_geneset(NULL, source = "reactome", species = "human")
gs_msig <- get_geneset(NULL, source = "msigdb", species = "human", collection = "H")
```
