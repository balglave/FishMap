unlink("pkgdown/assets", recursive = TRUE)
dir.create("pkgdown/assets", recursive = TRUE)

if (!requireNamespace("covr", quietly = TRUE)) install.packages("covr")
if (!requireNamespace("testdown", quietly = TRUE)) install.packages("testdown", repos = c("thinkropen" = "https://thinkr-open.r-universe.dev"))

# Create covr book and add it to pkgdown ----
covr::report(file = "pkgdown/assets/coverage.html", browse = FALSE)
file.copy("dev/codecoverage_explanation.md", "pkgdown/assets/codecoverage_explanation.md")

# Create testdown book and add it to pkgdown ----
testdown::test_down(open = FALSE)
file.copy("tests/testdown", "pkgdown/assets", recursive = TRUE)
unlink("tests/testdown", recursive = TRUE)

pkgdown::build_site()
