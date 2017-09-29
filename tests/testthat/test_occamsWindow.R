library("sugsvaRsel")
context("occamsWindow")


test_that("occamsWindow Output",{

  data(Isomix)
  Noise <- matrix(rnorm(250*4,0,1), nrow = 250, ncol = 4)
  IsoNoise<-cbind(Isomix$samplesIso,Noise)
  intfeatures <- rep(1,6)

  cc <- sugsvarsel(IsoNoise, 1, 2, intfeatures, 1, Model="ML")
  dd <- occamsWindow(cc$ML, 1)
  expect_true(all(dd>=0))
  expect_true(all(occamsWindow(c(0, 5, 0.5), 1)>=0))

}
)
