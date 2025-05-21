# testing smiles_pubchem

test_that("smiles_pubchem returns valid CID for known SMILES", {
  # Aspirin: CC(=O)OC1=CC=CC=C1C(=O)O → CID 2244
  input_df <- data.frame(V1 = c("CC(=O)OC1=CC=CC=C1C(=O)O"), stringsAsFactors = FALSE)
  result <- smiles_pubchem(input_df)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true("CID" %in% colnames(result))
  expect_false(result$CID[1] == "error")
})


test_that("smiles_pubchem returns 'error' for invalid SMILES", {
  input_df <- data.frame(V1 = c("C1CC1@@@INVALID"), stringsAsFactors = FALSE)
  result <- smiles_pubchem(input_df)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)  # filtered after grouping
})

test_that("smiles_pubchem processes mixed SMILES", {
  input_df <- data.frame(V1 = c("CC(=O)OC1=CC=CC=C1C(=O)O", "%%%", ""), stringsAsFactors = FALSE)
  result <- smiles_pubchem(input_df)
  
  expect_s3_class(result, "data.frame")
  expect_true(all(c("input", "DTXSID", "CID", "TITLE") %in% colnames(result)))
  expect_equal(nrow(result), 1)  # Only valid one remains
})

test_that("smiles_pubchem handles empty dataframe", {
  input_df <- data.frame(V1 = character(0), stringsAsFactors = FALSE)
  result <- smiles_pubchem(input_df)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

test_that("smiles_pubchem handles SMILES with spaces", {
  input_df <- data.frame(V1 = c("CC(=O) OC1=CC=CC=C1C(=O)O"), stringsAsFactors = FALSE)
  result <- smiles_pubchem(input_df)
  expect_equal(nrow(result), 0)  # gets removed after group_by + filtering
})
