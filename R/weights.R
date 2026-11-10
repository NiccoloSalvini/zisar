#' Spatial weights for a regular lattice
#'
#' Builds the contiguity weights matrix of an `nrow` x `ncol` regular grid,
#' the design used in the Monte Carlo study of the ZISAR paper. Cells are
#' numbered column by column, so the cell in row \eqn{i} and column \eqn{j}
#' is unit \eqn{i + (j - 1) \times} `nrow`, the order of
#' `expand.grid(row = 1:nrow, col = 1:ncol)`.
#'
#' @param nrow,ncol Number of rows and columns of the grid.
#' @param type `"queen"` (cells sharing an edge or a corner are neighbours)
#'   or `"rook"` (edge only).
#' @param style `"W"` for a row-standardised matrix (rows sum to one), `"B"`
#'   for the binary contiguity matrix.
#'
#' @return A sparse `nrow * ncol` square matrix of class `dgCMatrix`.
#' @export
#' @examples
#' W <- lattice_weights(4)
#' dim(W)
#' Matrix::rowSums(W)
lattice_weights <- function(nrow, ncol = nrow, type = c("queen", "rook"),
                            style = c("W", "B")) {
  type <- match.arg(type)
  style <- match.arg(style)
  stopifnot(nrow >= 1, ncol >= 1, nrow * ncol >= 2)
  n <- nrow * ncol
  id <- matrix(seq_len(n), nrow, ncol)

  # Half of the neighbour offsets; the other half comes from symmetry
  offs <- if (type == "rook") list(c(0, 1), c(1, 0)) else
    list(c(0, 1), c(1, 0), c(1, 1), c(1, -1))
  i <- j <- integer(0)
  for (o in offs) {
    rr <- seq_len(nrow)[seq_len(nrow) + o[1] <= nrow]
    cc <- seq_len(ncol)[seq_len(ncol) + o[2] >= 1 & seq_len(ncol) + o[2] <= ncol]
    if (!length(rr) || !length(cc)) next
    from <- as.vector(id[rr, cc, drop = FALSE])
    to <- as.vector(id[rr + o[1], cc + o[2], drop = FALSE])
    i <- c(i, from, to)
    j <- c(j, to, from)
  }
  B <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(n, n))
  if (style == "B") return(B)
  rs <- Matrix::rowSums(B)
  rs[rs == 0] <- 1
  as(Matrix::Diagonal(x = 1 / rs) %*% B, "CsparseMatrix")
}

# Coerce W (base matrix, Matrix, or spdep listw) to a checked dgCMatrix
.as_weights <- function(W, n) {
  if (inherits(W, "listw")) W <- .listw_to_sparse(W)
  if (!(is.matrix(W) || inherits(W, "Matrix")))
    stop("`W` must be a matrix, a sparse 'Matrix' or an spdep 'listw' object.",
         call. = FALSE)
  # Matrix::Matrix() first: it also loads the coercion methods of 'Matrix',
  # which as() needs in a session where 'Matrix' is not loaded yet
  if (is.matrix(W)) W <- Matrix::Matrix(W, sparse = TRUE)
  W <- as(as(as(W, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  if (nrow(W) != ncol(W))
    stop("`W` must be a square matrix.", call. = FALSE)
  if (nrow(W) != n)
    stop(sprintf("`W` is %d x %d but there are %d observations.",
                 nrow(W), ncol(W), n), call. = FALSE)
  if (anyNA(W@x) || any(!is.finite(W@x)))
    stop("`W` contains missing or infinite values.", call. = FALSE)
  if (any(Matrix::diag(W) != 0))
    stop("`W` must have a zero diagonal: no unit is its own neighbour.",
         call. = FALSE)
  rs <- Matrix::rowSums(W)
  if (any(abs(rs[rs != 0] - 1) > 1e-8))
    warning("`W` is not row-standardised. `rho` is restricted to (-1, 1), ",
            "the stationary range only when W is row-standardised.",
            call. = FALSE)
  W
}

# spdep listw -> dgCMatrix, without depending on spdep. A unit with no
# neighbours has neighbours[[i]] == 0L and a NULL weights entry.
.listw_to_sparse <- function(lw) {
  nb <- lw$neighbours
  n <- length(nb)
  nbrs <- lapply(nb, function(v) v[v > 0])
  i <- rep(seq_len(n), lengths(nbrs))
  Matrix::sparseMatrix(i = i, j = unlist(nbrs), x = unlist(lw$weights),
                       dims = c(n, n))
}
