# Export score table

Export score table

## Usage

``` r
export_scores(x, file, format = c("csv", "tsv"), include_meta = TRUE)
```

## Arguments

- x:

  `gleam_score` object or data.frame.

- file:

  Output file path.

- format:

  Output format: `csv` or `tsv`.

- include_meta:

  Whether to include metadata columns when input is score object.

## Value

Invisibly the exported data.frame.
