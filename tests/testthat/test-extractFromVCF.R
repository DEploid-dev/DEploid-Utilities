#' @title Extract read counts from VCF
#'
#' @description Extract read counts from VCF file of a single sample.
#'
#' @note The VCF file should only contain one sample. If more samples present in
#'  the VCF, it only returns coverage for of the first sample.
#'
#' @param vcfFileName Path of the VCF file.
#'
#' @param ADFieldIndex Index of the AD field of the sample field. For example,
#'  if the format is "GT:AD:DP:GQ:PL", the AD index is 2 (by default).
#'
#' @return A data.frame contains four columns: chromosomes, positions, reference
#'  allele count, alternative allele count.
#'
#' @examples
#' vcfFile <- system.file("extdata", "PG0390-C.test.vcf.gz", package = "DEploid")
#' PG0390 <- extractCoverageFromVcf(vcfFile)
#'
old_extractCoverageFromVcf <- function(vcfFileName, ADFieldIndex = 2) {
  # Assume that AD is the second field
  h <- function(w) {
    if (any(grepl("gzfile connection", w))) {
      invokeRestart("muffleWarning")
    }
  }

  gzf <- gzfile(vcfFileName, open = "rb")
  numberOfHeaderLines <- 0
  line <- withCallingHandlers(readLines(gzf, n = 1), warning = h)
  while (length(line) > 0) {
    if (grepl("##", line)) {
      numberOfHeaderLines <- numberOfHeaderLines + 1
    } else {
      break
    }
    line <- withCallingHandlers(readLines(gzf, n = 1), warning = h)
  }
  close(gzf)

  vcf <- read.table(gzfile(vcfFileName),
    skip = numberOfHeaderLines,
    header = T, comment.char = "", stringsAsFactors = FALSE,
    check.names = FALSE
  )

  sampleName <- names(vcf)[10]

  tmp <- vcf[[sampleName]]
  field <- strsplit(as.character(tmp), ":")

  tmpCovStr <- unlist(lapply(field, `[[`, ADFieldIndex))
  tmpCov <- strsplit(as.character(tmpCovStr), ",")

  refCount <- as.numeric(unlist(lapply(tmpCov, `[[`, 1)))
  altCount <- as.numeric(unlist(lapply(tmpCov, `[[`, 2)))

  return(data.frame(
    CHROM = vcf[, 1],
    POS = vcf[, 2],
    refCount = refCount,
    altCount = altCount
  ))
}


test_that("test R vs cpp vcf extract", {
  PG0390 <- old_extractCoverageFromVcf(vcfFileName)
  testthat::expect_equal(names(PG0390CoverageVcf), names(PG0390))
  testthat::expect_equal(PG0390CoverageVcf$CHROM, PG0390$CHROM)
  testthat::expect_equal(PG0390CoverageVcf$POS, PG0390$POS)
  testthat::expect_equal(PG0390CoverageVcf$refCount, PG0390$refCount)
  testthat::expect_equal(PG0390CoverageVcf$altCount, PG0390$altCount)
})
