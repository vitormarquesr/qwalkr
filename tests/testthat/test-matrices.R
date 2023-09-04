test_that("cartesian product works", {
  P2 <- rbind(c(0, 1),
              c(1, 0))
  P3 <- rbind(c(0, 1, 0),
              c(1, 0, 1),
              c(0, 1, 0))
  expect_equal(cartesian(P2), rbind(c(0, 1, 1, 0),
                                    c(1, 0, 0, 1),
                                    c(1, 0, 0, 1),
                                    c(0, 1, 1, 0)))


  expect_equal(cartesian(P2, P3), rbind(c(0, 1, 0, 1, 0, 0),
                                        c(1, 0, 1, 0, 1, 0),
                                        c(0, 1, 0, 0, 0, 1),
                                        c(1, 0, 0, 0, 1, 0),
                                        c(0, 1, 0, 1, 0, 1),
                                        c(0, 0, 1, 0, 1, 0)))
})

