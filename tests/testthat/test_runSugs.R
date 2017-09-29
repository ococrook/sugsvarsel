library("sugsvaRsel")
context("runSUGS")

test_that("SUGS Output",{

  data(Isomix)
  cc<-runSugs(iter = 2, Isomix$samplesIso, Model = "PML")
  expect_equal(cc$member[,1],c(1,1))
  expect_equal(max(cc$member),max(cc$cluster))
  expect_equal(length(cc$LPML),2)
  expect_equal(max(cc$ordering),nrow(Isomix$samplesIso))
  expect_equal(min(cc$ordering),1)
}
)



