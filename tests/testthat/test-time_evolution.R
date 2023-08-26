test_that("unitary evolution works", {
  walk_P3 <- ctqwalk(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))
  U_pst <- unitary_matrix(walk_P3, pi/sqrt(2))
  op_diag <- U_pst[row(U_pst) == (4 - col(U_pst))]

  expect_equal(op_diag, c(-1+ 0i, -1+0i, -1+0i))


  walk_P2 <- ctqwalk(rbind(c(0,1), c(1,0)))
  U_mix <- unitary_matrix(walk_P2, pi/4)

  expect_equal(U_mix, (1/sqrt(2))*rbind(c(1, 1i), c(1i, 1)))

})
