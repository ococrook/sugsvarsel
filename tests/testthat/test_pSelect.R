context("pSelect")


test_that("pSelect Output",{

  Noise <- matrix(rnorm(250*4,0,1), nrow = 250, ncol = 4)

  cc <- pSelect(Noise, iter = 2, p = 3, numSelect = 1, Model = "ML")
  expect_equal(ncol(cc), ncol(Noise))
  expect_equal(sum(cc %in% c(0,1)), ncol(Noise))
  expect_error(pSelect(Noise, iter = 1, p = 7, numSelect = 1, Model = "ML"))
}
)
