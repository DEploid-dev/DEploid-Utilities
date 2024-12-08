if (require("testthat")) {
  test_check("DEploid.utils")
} else {
  warning("testthat not available. Skipping unittests!")
}
