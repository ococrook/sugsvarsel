library("sugsvaRsel")
context("sugsvarsel")


test_that("sugsvarsel Output",{

  data(Isomix)
  Noise <- matrix(rnorm(250*4,0,1), nrow = 250, ncol = 4)
  IsoNoise <- cbind(Isomix$samplesIso,Noise)
  intfeatures <- rep(1, ncol(IsoNoise))

  cc <- sugsvarsel(IsoNoise, 1, 2, intfeatures, 1, Model="ML")
  expect_equal(ncol(cc$member), nrow(IsoNoise))
  expect_true(all(cc$features %in% c(0,1)))
  expect_equal(max(cc$ordering),nrow(IsoNoise))
}
)
