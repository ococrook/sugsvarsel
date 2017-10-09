context("runSUGS")

test_that("SUGS Output",{

  Noise <- matrix(rnorm(250*4,0,1), nrow = 250, ncol = 4)
  cc <- runSugs(iter = 2, Noise, Model = "PML")
  expect_equal(cc$member[,1],c(1,1))
  expect_equal(max(cc$member),max(cc$cluster))
  expect_equal(length(cc$LPML),2)
  expect_equal(max(cc$ordering),nrow(Noise))
  expect_equal(min(cc$ordering),1)
}
)



