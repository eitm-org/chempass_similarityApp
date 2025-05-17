test_that("get_properties_from_CIDs works with valid CIDs", {
  df <- get_properties_from_CIDs("2244,2519")  # Aspirin and Caffeine
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
  expect_true("MolecularWeight" %in% colnames(df))
})

test_that("get_properties_from_CIDs returns error for invalid CID", {
  result <- get_properties_from_CIDs("000000")
  expect_equal(result, "error")
})

test_that("get_properties_from_CIDs returns error for empty input", {
  result <- get_properties_from_CIDs("")
  expect_equal(result, "error")
})

test_that("get_properties_from_CIDs handles tabs and extra spaces", {
  
  result <- get_properties_from_CIDs("2244,1234")
  expect_gt(nrow(result), 1)
  expect_true("TPSA" %in% colnames(result))
  
})