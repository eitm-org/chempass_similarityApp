# testing dtxsid_pubchem

test_that("dtxsid_pubchem returns valid CID for known DTXSID", {
  input_df <- data.frame(V1 = c("DTXSID2020001"), stringsAsFactors = FALSE)
  result <- dtxsid_pubchem(input_df)
  
  expect_s3_class(result, "data.frame")  # Converted from pandas
  expect_equal(nrow(result), 1)
  expect_true("CID" %in% colnames(result))
  expect_false(result$CID[1] == "error")
})

test_that("dtxsid_pubchem handles invalid DTXSID", {
  input_df <- data.frame(V1 = c("DTXSID9999999"), stringsAsFactors = FALSE)
  result <- dtxsid_pubchem(input_df)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)  # Gets filtered by group_by(CID) %>% filter(n() == 1)
})

test_that("dtxsid_pubchem processes mixed valid and invalid DTXSIDs", {
  input_df <- data.frame(V1 = c("DTXSID2020001", "DTXSID9999999"), stringsAsFactors = FALSE)
  result <- dtxsid_pubchem(input_df)
  
  expect_s3_class(result, "data.frame")
  expect_true(all(c("DTXSID", "CID", "TITLE") %in% colnames(result)))
  expect_equal(nrow(result), 1)  # Only valid one remains after filtering
})

test_that("dtxsid_pubchem handles empty input dataframe", {
  input_df <- data.frame(V1 = character(0), stringsAsFactors = FALSE)
  result <- dtxsid_pubchem(input_df)
  
  expect_error(result, NA)  # Should not throw error
  expect_equal(nrow(result), 0)
})

