test_that("standard average mixing matrix works", {
  walk_P3 <- ctqwalk(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))
  M <- avg_matrix(walk_P3)
  M_exp <- (1/8)*rbind(c(3, 2, 3),
                       c(2, 4, 2),
                       c(3, 2, 3))
  expect_true(all(equal(M, M_exp)))

  walk_P4 <- ctqwalk(rbind(c(0, 1, 0, 0),
                           c(1, 0, 1, 0),
                           c(0, 1, 0, 1),
                           c(0, 0, 1, 0)))
  M <- avg_matrix(walk_P4)
  M_exp <- (1/10)*rbind(c(3, 2, 2, 3),
                       c(2, 3, 3, 2),
                       c(2, 3, 3, 2),
                       c(3, 2, 2, 3))
  expect_true(all(equal(M, M_exp)))

})

test_that("generalized average mixing matrix works", {
  walk_K3 <- ctqwalk(matrix(c(0,1,1,1,0,1,1,1,0), nrow=3))

  expect_equal(gavg_matrix(walk_K3, c(2*pi/9)), matrix(1/3, nrow=3, ncol=3))

  walk_P3 <- ctqwalk(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))

  theta1 <- acos(-1/3)/(2*sqrt(2))
  theta2 <- pi/sqrt(2) - theta1

  expect_equal(gavg_matrix(walk_P3, c(theta1, theta2)), matrix(1/3, nrow=3, ncol=3))

})



