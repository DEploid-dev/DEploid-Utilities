vcfFileName <- system.file("extdata", "PG0390-C.test.vcf.gz",
  package = "DEploid.utils"
)
plafFileName <- system.file("extdata", "labStrains.test.PLAF.txt",
  package = "DEploid.utils"
)
panelFileName <- system.file("extdata", "labStrains.test.panel.txt",
  package = "DEploid.utils"
)
refFileName <- system.file("extdata", "PG0390-C.test.ref", package = "DEploid.utils")
altFileName <- system.file("extdata", "PG0390-C.test.alt", package = "DEploid.utils")

PG0390CoverageVcf <- extractCoverageFromVcf(vcfFileName, "PG0390-C")
plaf <- extractPLAF(plafFileName)
