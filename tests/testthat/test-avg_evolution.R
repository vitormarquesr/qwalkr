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
