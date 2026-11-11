# Spatial weights for a regular lattice

Builds the contiguity weights matrix of an `nrow` x `ncol` regular grid,
the design used in the Monte Carlo study of the ZISAR paper. Cells are
numbered column by column, so the cell in row \\i\\ and column \\j\\ is
unit \\i + (j - 1) \times\\ `nrow`, the order of
`expand.grid(row = 1:nrow, col = 1:ncol)`.

## Usage

``` r
lattice_weights(
  nrow,
  ncol = nrow,
  type = c("queen", "rook"),
  style = c("W", "B")
)
```

## Arguments

- nrow, ncol:

  Number of rows and columns of the grid.

- type:

  `"queen"` (cells sharing an edge or a corner are neighbours) or
  `"rook"` (edge only).

- style:

  `"W"` for a row-standardised matrix (rows sum to one), `"B"` for the
  binary contiguity matrix.

## Value

A sparse `nrow * ncol` square matrix of class `dgCMatrix`.

## Examples

``` r
W <- lattice_weights(4)
dim(W)
#> [1] 16 16
Matrix::rowSums(W)
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
```
