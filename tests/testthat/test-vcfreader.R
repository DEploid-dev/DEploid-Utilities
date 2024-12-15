test_that("Extracted coverage", {
  PG0390CoverageTxt <- extractCoverageFromTxt(refFileName, altFileName)
  expect_that(PG0390CoverageTxt, is_a("data.frame"))
  expect_that(PG0390CoverageVcf, is_a("data.frame"))
  expect_equal(PG0390CoverageTxt, PG0390CoverageVcf)
})


test_that("Extracted plaf", {
  expect_that(plaf, is_a("numeric"))
})


test_that("computeObsWSAF", {
  expect_equal(computeObsWSAF(0, 0), 0)
  expect_equal(computeObsWSAF(0, 100), 0)
  expect_equal(computeObsWSAF(1, 99), 0.01)
  expect_equal(computeObsWSAF(99, 1), 0.99)
  expect_equal(computeObsWSAF(50, 50), 0.5)
  expect_equal(computeObsWSAF(50, 100), 0.3333333333333)
})


test_that("WSAF Related", {
  obsWSAF <- computeObsWSAF(
    PG0390CoverageVcf$altCount,
    PG0390CoverageVcf$refCount
  )
  potentialOutliers <- c(5, 12, 25, 30, 35, 50)

  png(paste0(tempdir(), "/histWSAF.png"))
  expect_that(histWSAF(obsWSAF), is_a("histogram"))
  dev.off()

  ####
  png(paste0(tempdir(), "/WSAFvsPLAF.png"))
  expect_null(plotWSAFvsPLAF(plaf, obsWSAF))
  dev.off()
})


test_that("plotAltVsRef", {
  png(paste0(tempdir(), "/AltVsRef.png"))
  expect_null(plotAltVsRef(
    PG0390CoverageVcf$refCount,
    PG0390CoverageVcf$altCount
  ))
  dev.off()
})


test_that("plotAltVsRefWithOutliers", {
  potentialOutliers <- c(1, 10, 20, 30, 40)
  png(paste0(tempdir(), "/AltVsRefOutlier.png"))
  expect_null(plotAltVsRef(PG0390CoverageVcf$refCount,
                           PG0390CoverageVcf$altCount,
                           potentialOutliers = potentialOutliers
  ))
  dev.off()
})
