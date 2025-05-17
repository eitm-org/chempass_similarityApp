# File: tests/testthat/test-sanitize_names.R


test_that("sanitize_names handles complex compound names", {
  expect_equal(sanitize_names("N-(3-chlorophenyl)-1H-indole"), "N_3_chlorophenyl_1H_indole")
  expect_equal(sanitize_names("4'-hydroxyacetophenone"), "4_hydroxyacetophenone")
  expect_equal(sanitize_names("(+)-catechin"), "catechin")
  expect_equal(sanitize_names("α-methyltryptamine"), "methyltryptamine")  # if non-ASCII stripped
  expect_equal(sanitize_names("O,O-dimethyl-O-(3-methylphenyl) phosphorothioate"), "O_O_dimethyl_O_3_methylphenyl_phosphorothioate")
  expect_equal(sanitize_names("γ-aminobutyric acid"), "aminobutyric_acid")
  expect_equal(sanitize_names("(+)-catechin"), "catechin")
})

