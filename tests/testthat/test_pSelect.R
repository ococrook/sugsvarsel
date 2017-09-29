library("sugsvaRsel")
context("pSelect")


test_that("pSelect Output",{

  data(Isomix)
  Noise <- matrix(rnorm(250*4,0,1), nrow = 250, ncol = 4)
  IsoNoise<-cbind(Isomix$samplesIso,Noise)

  cc <- pSelect(IsoNoise, iter = 2, p = 3, numSelect = 1, Model = "ML")
  expect_equal(ncol(cc), ncol(IsoNoise))
  expect_equal(sum(cc %in% c(0,1)), ncol(IsoNoise))
  expect_error(pSelect(IsoNoise, iter = 1, p = 7, numSelect = 1, Model = "ML"))
}
)
