test_that("zisar_basis() builds a normalised basis and knot weights", {
  coords <- expand.grid(x = 1:20, y = 1:20)
  b <- zisar_basis(coords, h = 4)
  expect_s3_class(b, "zisar_basis")
  expect_equal(dim(b$B), c(400, b$K))
  expect_equal(Matrix::rowSums(b$B), rep(1, 400))
  expect_true(all(b$B@x >= 0))
  expect_equal(dim(b$Wk), c(b$K, b$K))
  expect_equal(Matrix::rowSums(b$Wk), rep(1, b$K))
  expect_equal(Matrix::diag(b$Wk), rep(0, b$K))
  expect_equal(max(b$eig_re), 1, tolerance = 1e-8)
  expect_equal(nrow(b$knots), b$K)
  expect_equal(b$r, 6)
  # Every kept knot has a unit within h
  d <- sapply(seq_len(b$K), function(k)
    min((coords$x - b$knots$x[k])^2 + (coords$y - b$knots$y[k])^2))
  expect_true(all(d <= 16))
  expect_output(print(b), "knots, spacing 4")
})

test_that("zisar_basis() links isolated knots and checks its arguments", {
  # Two clusters far apart: their knots are linked to the nearest knot
  coords <- rbind(expand.grid(x = 0:2, y = 0:2),
                  expand.grid(x = 30:32, y = 30:32))
  b <- zisar_basis(coords, h = 2)
  expect_true(all(Matrix::rowSums(b$Wk > 0) >= 1))
  expect_error(zisar_basis(coords[, 1, drop = FALSE], h = 2), "two columns")
  expect_error(zisar_basis(coords, h = -1), "positive")
  expect_error(zisar_basis(coords, h = 2, r_fac = 0.5), "at least 1")
})
