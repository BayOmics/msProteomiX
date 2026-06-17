# ==============================================================================
# msProteomiX — Basic Smoke Tests
# ==============================================================================
# These tests verify that core functions are loadable and return expected types.
# No test data needed — uses built-in defaults only.

test_that("package loads successfully", {
  expect_true(requireNamespace("msProteomiX", quietly = TRUE))
})

test_that("calc_pI_seq returns numeric", {
  result <- calc_pI_seq(c("ACDEFGHIK", "PEPTIDE"))
  expect_type(result, "double")
  expect_length(result, 2)
  expect_true(all(result > 0 & result < 14))
})

test_that("calc_gravy returns numeric", {
  result <- calc_gravy(c("ACDEFGHIK", "PEPTIDE"))
  expect_type(result, "double")
  expect_length(result, 2)
})

test_that("calc_missed_cleavage handles edge cases", {
  # Tryptic peptide with no missed cleavage
  expect_equal(calc_missed_cleavage("PEPTIDEK"), 0L)
  # One missed cleavage (K in middle)
  expect_equal(calc_missed_cleavage("PEPKIDER"), 1L)
  # NA input
  expect_true(is.na(calc_missed_cleavage(NA_character_)))
})
