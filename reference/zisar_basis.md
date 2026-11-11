# Basis of a reduced-rank inflation field

Builds the knot lattice and basis functions of the reduced-rank field
\\g = B\alpha\\ used by `zisar(..., zero_spatial = "reduced")`.

## Usage

``` r
zisar_basis(coords, h, r_fac = 1.5)

# S3 method for class 'zisar_basis'
print(x, ...)
```

## Arguments

- coords:

  Two-column matrix or data frame of unit coordinates (for example
  centroids), in the order of the rows of `W`. Use projected
  coordinates, in the same unit as `h`.

- h:

  Spacing of the knot lattice.

- r_fac:

  Radius of the basis functions, as a multiple of `h`.

- x:

  A `"zisar_basis"` object.

- ...:

  Not used.

## Value

An object of class `"zisar_basis"`: a list with the basis `B` (\\n
\times K\\, sparse), the knot weights `Wk` (\\K \times K\\, sparse), the
eigenvalues of `Wk` (`eig_re`, `eig_im`), the knot coordinates `knots`,
the number of knots `K`, the spacing `h` and the radius `r`.

## Details

Knots lie on a square lattice of spacing `h` laid over the units, and
are kept where a unit lies within `h` of the knot. Column \\k\\ of \\B\\
is the bisquare function \\(1 - d^2/r^2)^2\\ of the distance \\d\\ from
knot \\k\\, zero beyond the radius \\r\\ = `r_fac * h`, and the rows of
\\B\\ are normalised to sum to one, so that \\g_i\\ is a weighted mean
of the knot values around unit \\i\\ (the fixed-rank representation of
Cressie and Johannesson, 2008). The knot weights \\W_K\\ are rook
contiguity on the lattice, row-standardised, with isolated knots linked
to their nearest knot.

The spacing sets the resolution of the field: variation at scales below
`h` is smoothed out. Each knot value should be informed by many units
(the ZISAR paper uses knots 10 km apart for English neighbourhoods,
about 150 units within reach of each knot); with few units per knot the
Laplace approximation degrades towards that of the unit-level field.

## References

Cressie, N. and Johannesson, G. (2008). Fixed rank kriging for very
large spatial data sets. *Journal of the Royal Statistical Society:
Series B*, 70(1), 209–226.
[doi:10.1111/j.1467-9868.2007.00633.x](https://doi.org/10.1111/j.1467-9868.2007.00633.x)

## See also

[`zisar()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar.md)

## Examples

``` r
# 20 x 20 lattice of units, knots every 4 units
coords <- expand.grid(x = 1:20, y = 1:20)
b <- zisar_basis(coords, h = 4)
b
#> Reduced-rank basis: 36 knots, spacing 4, radius 6, 400 units (6.2 knots per unit)
# Each unit's field value is a weighted mean of nearby knot values
range(Matrix::rowSums(b$B))
#> [1] 1 1
```
