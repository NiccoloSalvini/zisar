test_that("lattice_weights() builds queen and rook contiguity", {
  B <- lattice_weights(4, type = "queen", style = "B")
  expect_equal(dim(B), c(16, 16))
  expect_true(Matrix::isSymmetric(B))
  expect_equal(sum(Matrix::diag(B)), 0)
  # Interior cell (2, 2) = unit 6 has 8 queen neighbours, a corner 3
  expect_equal(sum(B[6, ]), 8)
  expect_equal(sum(B[1, ]), 3)

  R <- lattice_weights(4, type = "rook", style = "B")
  expect_equal(sum(R[6, ]), 4)
  expect_equal(sum(R[1, ]), 2)

  W <- lattice_weights(3, 5)
  expect_equal(dim(W), c(15, 15))
  expect_equal(Matrix::rowSums(W), rep(1, 15))
})

test_that("lattice_weights() numbers cells column by column", {
  B <- lattice_weights(3, 4, type = "rook", style = "B")
  # Unit 1 = (1, 1): neighbours (2, 1) = unit 2 and (1, 2) = unit 4
  expect_equal(which(B[1, ] == 1), c(2, 4))
})

test_that("weights are validated", {
  W <- as.matrix(lattice_weights(3))
  expect_error(zisar:::.as_weights(W[, -1], 9), "square")
  expect_error(zisar:::.as_weights(W, 10), "9 x 9")
  diag(W)[1] <- 0.5
  expect_error(zisar:::.as_weights(W, 9), "zero diagonal")
  expect_error(zisar:::.as_weights(data.frame(a = 1), 1), "must be a matrix")
  expect_warning(
    zisar:::.as_weights(lattice_weights(3, style = "B"), 9),
    "not row-standardised"
  )
  expect_s4_class(zisar:::.as_weights(as.matrix(lattice_weights(3)), 9),
                  "dgCMatrix")
})

test_that("spdep listw objects are converted, islands included", {
  skip_if_not_installed("spdep")
  lw <- spdep::nb2listw(spdep::cell2nb(4, 4, type = "queen"), style = "W")
  W <- zisar:::.as_weights(lw, 16)
  expect_equal(as.matrix(W), spdep::listw2mat(lw), ignore_attr = TRUE)

  nb <- spdep::cell2nb(3, 3)
  nb[[9]] <- 0L
  nb[[6]] <- setdiff(nb[[6]], 9L)
  nb[[8]] <- setdiff(nb[[8]], 9L)
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  W <- zisar:::.listw_to_sparse(lw)
  expect_equal(sum(W[9, ]), 0)
  expect_equal(as.matrix(W), spdep::listw2mat(lw), ignore_attr = TRUE)
})
