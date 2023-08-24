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
  sum_proj <- get_eigproj(K3, 1) + get_eigproj(K3, 2)
  I3 <- diag(3)

  expect_true(all(equal(sum_proj, I3)))

  A_h <- spectral(matrix(c(0,1-1i,1-2i,1+1i,0,1+3i, 1+2i,1-3i,0), nrow=3), multiplicity = TRUE)
  sum_proj_h <- get_eigproj(A_h, 1) + get_eigproj(A_h, 2) + get_eigproj(A_h, 3)
  I3 <- diag(3)

  expect_true(all(equal(sum_proj_h, I3)))

})

test_that("schur product of projectors work", {
  K3 <- spectral(matrix(c(0,1,1,1,0,1, 1, 1, 0), nrow=3), multiplicity = TRUE)
  E1_squared <- extractSCHUR(K3, id1=1)

  expect_true(all(equal(E1_squared, 1/9)))
})

test_that("matrix functions work", {
  H <- matrix(c(0,1,1,1,0,1,1,1,0), nrow=3)
  decomp <- spectral(H)

  expect_equal(evalMFUN(decomp, FUN = function(x) x^3), H %*% H %*% H)
  expect_equal(evalMFUN(decomp, FUN = function(x) 1/x), solve(H))

  H_cmplx <- matrix(c(0,1+3i,1-2i,1-3i,0,3,1+2i,3,5), nrow=3)
  decomp_cmplx <- spectral(H_cmplx)
  expect_equal(evalMFUN(decomp_cmplx, function(x,y) x^y, 2), H_cmplx %*% H_cmplx)

})




