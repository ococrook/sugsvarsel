context("sugsvarsel")


test_that("sugsvarsel Output",{


  Noise <- matrix(rnorm(250*4,0,1), nrow = 250, ncol = 4)
  intfeatures <- rep(1, ncol(Noise))

  cc <- sugsvarsel(Noise, 1, 2, intfeatures, 1, Model="ML")
  expect_equal(ncol(cc$member), nrow(Noise))
  expect_true(all(cc$features %in% c(0,1)))
  expect_equal(max(cc$ordering),nrow(Noise))
}
)
