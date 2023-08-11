test_that("inference on eigenvalue multiplicity works", {
  Ak33 <- matrix(c(0,0,0,1,1,1,
                           0,0,0,1,1,1,
                           0,0,0,1,1,1,
                           1,1,1,0,0,0,
                           1,1,1,0,0,0,
                           1,1,1,0,0,0), nrow=6)

  expect_equal(spectral(Ak33, multiplicity=TRUE)$multiplicity, c(1,4,1))
  expect_equal(spectral(Ak33, multiplicity=FALSE)$multiplicity, rep(1, times=6))
  expect_equal(spectral(Ak33, multiplicity=TRUE, tol=4)$multiplicity, 6)
  expect_equal(spectral(matrix(c(1)))$multiplicity, 1)
})

test_that("extraction of projectors work", {
  K3 <- spectral(matrix(c(0,1,1,1,0,1, 1,1,0), nrow=3), multiplicity = TRUE)
  sum_proj <- extractPROJ(K3, 1) + extractPROJ(K3, 2)
  I3 <- diag(3)

  expect_true(all(equal(sum_proj, I3)))
})

test_that("schur product of projectors work", {
  K3 <- spectral(matrix(c(0,1,1,1,0,1, 1, 1, 0), nrow=3), multiplicity = TRUE)
  E1_squared <- extractSCHUR(K3, id1=1)

  expect_true(all(equal(E1_squared, 1/9)))
})





