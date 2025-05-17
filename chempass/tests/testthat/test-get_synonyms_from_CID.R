# get_synonyms_from_CID function

test_that("get_synonyms_from_CID works with valid CIDs", {
  df <- get_properties_from_CIDs("2244,2519")  # Aspirin and Caffeine
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

test_that("get_synonyms_from_CID gives errors with incorrect CID", {
  df <- get_properties_from_CIDs("00snc")  
  expect_equal(df, "error")
})

test_that("get_synonyms_from_CID returns error for empty input", {
  result <- get_synonyms_from_CID("")
  expect_equal(result, "error")
})

test_that("get_synonyms_from_CID returns <= 10 synonyms", {
  result <- get_synonyms_from_CID("1983")  # Caffeine
  expect_type(result, "character")
  expect_lte(length(strsplit(result, ";")[[1]]), 10)
})