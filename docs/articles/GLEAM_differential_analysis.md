# GLEAM differential analysis guide

## Comparison levels

[`test_signature()`](https://JamesWu7.github.io/GLEAM/reference/test_signature.md)
supports multiple levels:

- `cell` (exploratory)
- `sample`
- `celltype`
- `sample_celltype`
- `pseudobulk`
- `region`
- `sample_region`
- `trajectory`

## Common patterns

``` r

res_cell <- test_signature(sc, group = group_col, level = "cell", method = "wilcox")
res_sample <- test_signature(sc, group = group_col, sample = "sample", level = "sample", method = "wilcox")
res_celltype <- test_signature(sc, celltype = celltype_col, group = group_col, level = "celltype", method = "wilcox")
res_sample_celltype <- test_signature(sc, group = group_col, sample = "sample", celltype = celltype_col, level = "sample_celltype")
res_pb <- test_signature(sc, group = group_col, sample = "sample", celltype = celltype_col, level = "pseudobulk")

res_wct <- test_signature(
  score = sc,
  group = group_col,
  sample = "sample",
  celltype = celltype_col,
  target_celltype = "CD8_T",
  level = "sample"
)
```

## Differential result visualizations

``` r

plot_volcano(res_pb)
plot_celltype_compare(res_celltype)
plot_group_in_celltype(res_wct)
```
