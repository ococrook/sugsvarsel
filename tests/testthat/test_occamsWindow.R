context("occamsWindow")


test_that("occamsWindow Output",{

  Noise <- matrix(rnorm(250*4,0,1), nrow = 250, ncol = 4)
  intfeatures <- rep(1,4)

  cc <- sugsvarsel(Noise, 1, 2, intfeatures, 1, Model="ML")
  dd <- occamsWindow(cc$ML, 1)
  expect_true(all(dd>=0))
  expect_true(all(occamsWindow(c(0, 5, 0.5), 1)>=0))

}
)
