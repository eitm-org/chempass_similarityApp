# test input_cid_reformating

test_that("input_cid_reformating returns title for valid CID", {
  input_df <- data.frame(V1 = c("2244"), stringsAsFactors = FALSE)  # Aspirin
  result <- input_cid_reformating(input_df)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true("TITLE" %in% colnames(result))
  expect_false(result$TITLE[1] == "error")
  expect_true(grepl("aspirin", tolower(result$TITLE[1])))
})

test_that("input_cid_reformating returns 'error' for invalid CID", {
  input_df <- data.frame(V1 = c("999999999"), stringsAsFactors = FALSE)  # Invalid CID
  result <- input_cid_reformating(input_df)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)  # filtered out after group_by
})

test_that("input_cid_reformating processes valid and invalid CIDs", {
  input_df <- data.frame(V1 = c("2244", "999999999"), stringsAsFactors = FALSE)
  result <- input_cid_reformating(input_df)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)  # only valid retained
  expect_true("aspirin" %in% tolower(result$TITLE))
})

test_that("input_cid_reformating handles empty input", {
  input_df <- data.frame(V1 = character(0), stringsAsFactors = FALSE)
  result <- input_cid_reformating(input_df)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

test_that("input_cid_reformating handles malformed CIDs", {
  input_df <- data.frame(V1 = c("ABC123", "12 34", "-123"), stringsAsFactors = FALSE)
  result <- input_cid_reformating(input_df)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)  # all should fail
})
