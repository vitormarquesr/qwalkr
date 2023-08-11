test_that("unique values are tagged correctly", {
  expect_equal(unique_eigvals(c(1, 2.1, 2.1, 3, 4.5)), c(TRUE, TRUE, FALSE, TRUE, TRUE))
  expect_equal(unique_eigvals(c(1, 1, 3, 2, sqrt(2)^2)), c(TRUE, FALSE, TRUE, TRUE, FALSE))
  expect_equal(unique_eigvals(c(1)), c(TRUE))
})

test_that("multiplicity calculation from tagged values work", {
  expect_equal(mult_eigvals(unique_eigvals(c(1,2,2,3))), c(1,2,1))
  expect_equal(mult_eigvals(unique_eigvals(c(1,1,2,3,3,3))), c(2,1,3))
  expect_equal(mult_eigvals(unique_eigvals(c(1))), 1)
})

test_that("extraction of eigenspace index works", {
  expect_equal(index_eigspace(c(1,1,2,1,3), 3), c(3,4))
  expect_equal(index_eigspace(c(1,1,2,1,3), 5), c(6,7,8))
  expect_equal(index_eigspace(c(1,1), 2), c(2))
  expect_equal(index_eigspace(c(1), 1), 1)
})


