#' Basis of a reduced-rank inflation field
#'
#' Builds the knot lattice and basis functions of the reduced-rank field
#' \eqn{g = B\alpha} used by `zisar(..., zero_spatial = "reduced")`.
#'
#' Knots lie on a square lattice of spacing `h` laid over the units, and are
#' kept where a unit lies within `h` of the knot. Column \eqn{k} of \eqn{B} is
#' the bisquare function \eqn{(1 - d^2/r^2)^2} of the distance \eqn{d} from
#' knot \eqn{k}, zero beyond the radius \eqn{r} = `r_fac * h`, and the rows of
#' \eqn{B} are normalised to sum to one, so that \eqn{g_i} is a weighted mean
#' of the knot values around unit \eqn{i} (the fixed-rank representation of
#' Cressie and Johannesson, 2008). The knot weights \eqn{W_K} are rook
#' contiguity on the lattice, row-standardised, with isolated knots linked to
#' their nearest knot.
#'
#' The spacing sets the resolution of the field: variation at scales below
#' `h` is smoothed out. Each knot value should be informed by many units
#' (the ZISAR paper uses knots 10 km apart for English neighbourhoods, about
#' 150 units within reach of each knot); with few units per knot the Laplace
#' approximation degrades towards that of the unit-level field.
#'
#' @param coords Two-column matrix or data frame of unit coordinates (for
#'   example centroids), in the order of the rows of `W`. Use projected
#'   coordinates, in the same unit as `h`.
#' @param h Spacing of the knot lattice.
#' @param r_fac Radius of the basis functions, as a multiple of `h`.
#'
#' @return An object of class `"zisar_basis"`: a list with the basis `B`
#'   (\eqn{n \times K}, sparse), the knot weights `Wk` (\eqn{K \times K},
#'   sparse), the eigenvalues of `Wk` (`eig_re`, `eig_im`), the knot
#'   coordinates `knots`, the number of knots `K`, the spacing `h` and the
#'   radius `r`.
#'
#' @references
#' Cressie, N. and Johannesson, G. (2008). Fixed rank kriging for very large
#' spatial data sets. *Journal of the Royal Statistical Society: Series B*,
#' 70(1), 209--226. \doi{10.1111/j.1467-9868.2007.00633.x}
#' @seealso [zisar()]
#' @export
#' @examples
#' # 20 x 20 lattice of units, knots every 4 units
#' coords <- expand.grid(x = 1:20, y = 1:20)
#' b <- zisar_basis(coords, h = 4)
#' b
#' # Each unit's field value is a weighted mean of nearby knot values
#' range(Matrix::rowSums(b$B))
zisar_basis <- function(coords, h, r_fac = 1.5) {
  coords <- as.matrix(coords)
  if (!is.numeric(coords) || ncol(coords) != 2L || any(!is.finite(coords)))
    stop("`coords` must be a finite numeric matrix with two columns.",
         call. = FALSE)
  if (length(h) != 1L || !is.finite(h) || h <= 0)
    stop("`h` must be a positive number.", call. = FALSE)
  if (length(r_fac) != 1L || !is.finite(r_fac) || r_fac < 1)
    stop("`r_fac` must be a number of at least 1.", call. = FALSE)
  east <- coords[, 1]
  north <- coords[, 2]

  x0 <- floor(min(east) / h) * h
  y0 <- floor(min(north) / h) * h
  nx <- ceiling((max(east) - x0) / h) + 1
  ny <- ceiling((max(north) - y0) / h) + 1
  kn <- expand.grid(i = seq_len(nx) - 1, j = seq_len(ny) - 1)
  kn$x <- x0 + kn$i * h
  kn$y <- y0 + kn$j * h
  # Units indexed by their nearest lattice cell
  cells <- split(seq_along(east), paste(round((east - x0) / h),
                                        round((north - y0) / h)))
  units_near <- function(i, j, m) {
    keys <- as.vector(outer(i + (-m:m), j + (-m:m), paste))
    unlist(cells[keys], use.names = FALSE)
  }
  keep <- vapply(seq_len(nrow(kn)), function(k) {
    u <- units_near(kn$i[k], kn$j[k], 1)
    length(u) > 0 &&
      min((east[u] - kn$x[k])^2 + (north[u] - kn$y[k])^2) <= h^2
  }, logical(1))
  kn <- kn[keep, ]
  K <- nrow(kn)
  r <- r_fac * h
  m_r <- ceiling(r_fac + 0.5)
  trip <- do.call(rbind, lapply(seq_len(K), function(k) {
    u <- units_near(kn$i[k], kn$j[k], m_r)
    d2 <- (east[u] - kn$x[k])^2 + (north[u] - kn$y[k])^2
    s <- d2 < r^2
    if (!any(s)) return(NULL)
    cbind(u[s], k, (1 - d2[s] / r^2)^2)
  }))
  B <- Matrix::sparseMatrix(i = trip[, 1], j = trip[, 2], x = trip[, 3],
                            dims = c(length(east), K))
  if (any(Matrix::rowSums(B) == 0))
    stop("Some units are not covered by the basis: increase `r_fac`.",
         call. = FALSE)
  B <- Matrix::Diagonal(x = 1 / Matrix::rowSums(B)) %*% B

  key <- paste(kn$i, kn$j)
  nbs <- lapply(seq_len(K), function(k) {
    v <- match(paste(kn$i[k] + c(1, -1, 0, 0), kn$j[k] + c(0, 0, 1, -1)), key)
    v[!is.na(v)]
  })
  for (k in which(lengths(nbs) == 0)) {
    d <- (kn$x - kn$x[k])^2 + (kn$y - kn$y[k])^2
    d[k] <- Inf
    j <- which.min(d)
    nbs[[k]] <- j
    nbs[[j]] <- c(nbs[[j]], k)
  }
  A <- Matrix::sparseMatrix(i = rep(seq_len(K), lengths(nbs)),
                            j = unlist(nbs), x = 1, dims = c(K, K))
  A <- (A + Matrix::t(A) > 0) * 1
  Wk <- Matrix::Diagonal(x = 1 / Matrix::rowSums(A)) %*% A
  ev <- eigen(as.matrix(Wk), only.values = TRUE)$values
  rownames(kn) <- NULL
  structure(
    list(B = as(as(B, "CsparseMatrix"), "TsparseMatrix"),
         Wk = as(as(Wk, "CsparseMatrix"), "TsparseMatrix"),
         eig_re = Re(ev), eig_im = Im(ev), knots = kn[, c("x", "y")],
         K = K, h = h, r = r),
    class = "zisar_basis")
}

#' @rdname zisar_basis
#' @param x A `"zisar_basis"` object.
#' @param ... Not used.
#' @export
print.zisar_basis <- function(x, ...) {
  cat(sprintf(paste0("Reduced-rank basis: %d knots, spacing %s, radius %s, ",
                     "%d units (%.1f knots per unit)\n"),
              x$K, format(x$h), format(x$r), nrow(x$B),
              Matrix::nnzero(x$B) / nrow(x$B)))
  invisible(x)
}
